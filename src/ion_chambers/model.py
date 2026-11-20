"""The ion chamber model: charge carriers drifting to the electrodes."""

import festim as F
import numpy as np
from dolfinx import fem

from ion_chambers.cellwise_diffusion_material import CellwiseDiffusionMaterial
from ion_chambers.constants import TOWNSEND
from ion_chambers.field import LaplaceField, PotentialSolver, SpaceChargePerturbation
from ion_chambers.gas import Coefficient, Gas
from ion_chambers.geometry import Geometry
from ion_chambers.mobility_electromigration import MobilityElectromigration
from ion_chambers.readout import Reading, compute_reading
from ion_chambers.source import IonisationSource
from ion_chambers.stabilisation import DiffusionFloor, cell_field_length

ANODE_ID = 1
CATHODE_ID = 2


class IonChamberModel:
    r"""Steady-state charge transport in a gas-filled ionisation chamber.

    Each charge carrier of the gas follows the Nernst-Planck flux

    .. math::

        J = -D \nabla c + c\, v, \qquad v = -z \mu \nabla \varphi

    in a given electric potential :math:`\varphi`, is created by the ionisation
    sources and reacts with the other carriers. The model is assembled and solved with
    FESTIM.

    The inputs are plain objects. Nothing is meshed or solved until
    :meth:`initialise`; :meth:`run` then solves the steady state and returns the
    :class:`Reading`. Calling :meth:`run` again after changing the activity of a
    source reuses the mesh, the potential and the wall-loss calculation.

    Args:
        geometry: the chamber geometry
        gas: the gas filling the chamber
        sources: the ionisation source, or a list of them
        bias: the potential of the anode relative to the cathode (V). It can be
            negative.
        resolution: the mesh resolution, passed to the geometry
        field: the potential solver. Defaults to :class:`LaplaceField`.
        stabilisation: the stabilisation of the drift. Defaults to
            :class:`DiffusionFloor`.

    Attributes:
        geometry: the chamber geometry
        gas: the gas filling the chamber
        sources: the ionisation sources
        bias: the potential of the anode relative to the cathode (V)
        resolution: the mesh resolution
        field: the potential solver
        stabilisation: the stabilisation of the drift
        mesh: the mesh, after ``initialise``
        potential: the electric potential (V), after ``initialise``
        reduced_field: E/N in each cell (Td), after ``initialise``
        species: the FESTIM species of each carrier, keyed by name
        mobilities: the mobility of each carrier (m2/V/s), a float or a cell-wise field
        diffusivities: the stabilised diffusion coefficient of each carrier (m2/s), a
            cell-wise field
        rates: the rate coefficient of each reaction, in the order of
            ``gas.reactions``
        pair_production: the total production rate of ion pairs (1/m3/s)
        festim_problem: the underlying :class:`festim.HydrogenTransportProblem`

    Examples:
        .. code-block:: python

            import ion_chambers as ic

            model = ic.IonChamberModel(
                geometry=ic.CoaxialChamber.from_volume(
                    volume=1e-3, rod_radius=1.5e-3, can_radius=0.05
                ),
                gas=ic.library.air(),
                sources=ic.TritiumInGas(activity=3.7e10),
                bias=500.0,
            )
            model.initialise()
            reading = model.run()
            print(reading.current)
    """

    def __init__(
        self,
        geometry: Geometry,
        gas: Gas,
        sources: IonisationSource | list[IonisationSource],
        bias: float = 500.0,
        resolution: int = 120,
        field: PotentialSolver | None = None,
        stabilisation: DiffusionFloor | None = None,
    ) -> None:
        self.geometry = geometry
        self.gas = gas
        self.sources = sources if isinstance(sources, list) else [sources]
        self.bias = bias
        self.resolution = resolution
        self.field = field or LaplaceField()
        self.stabilisation = stabilisation or DiffusionFloor()

        self.mesh = None
        self.potential = None
        self.festim_problem = None

    @property
    def bias(self) -> float:
        return self._bias

    @bias.setter
    def bias(self, value: float) -> None:
        if isinstance(value, bool) or not isinstance(value, int | float):
            raise TypeError(f"bias must be a float, not {type(value).__name__}")
        if value == 0:
            raise ValueError("bias must be non-zero: without a field nothing drifts")
        self._bias = float(value)

    def initialise(self) -> None:
        """Builds the mesh, the potential, the coefficients and the FESTIM problem."""
        self.mesh = self.geometry.build_mesh(self.resolution)
        self.potential = self.field.solve(self.mesh, self.geometry, self.bias)
        self.space_charge = SpaceChargePerturbation(self.mesh, self.geometry)
        self._create_coefficients()
        self.pair_production = sum(
            source.build(self.mesh, self.geometry, self.gas) for source in self.sources
        )
        self._create_festim_problem()

    def run(self) -> Reading:
        """Solves the steady state and returns the reading."""
        if self.festim_problem is None:
            raise RuntimeError("call initialise() before run()")
        self.festim_problem.run()
        return compute_reading(self)

    def concentrations(self) -> dict[str, fem.Function]:
        """The number density of each carrier (1/m3), keyed by name."""
        return {
            name: species.post_processing_solution
            for name, species in self.species.items()
        }

    def _create_coefficients(self) -> None:
        """Evaluates mobilities, diffusivities and rates in the local field."""
        field_strength, cell_length = cell_field_length(self.mesh, self.potential)
        self.reduced_field = field_strength / self.gas.number_density / TOWNSEND
        V0 = fem.functionspace(self.mesh, ("DG", 0))

        self.mobilities, self.diffusivities = {}, {}
        for carrier in self.gas.carriers:
            mobility = self._in_cells(carrier.mobility)
            if carrier.diffusivity is None:
                physical = mobility * self.gas.thermal_voltage
            else:
                physical = self._in_cells(carrier.diffusivity)
            D = fem.Function(V0, name=f"D_{carrier.name}")
            D.x.array[:] = self.stabilisation.apply(
                physical, mobility, field_strength, cell_length
            )
            self.diffusivities[carrier.name] = D
            self.mobilities[carrier.name] = self._as_coefficient(
                carrier.mobility, V0, f"mu_{carrier.name}"
            )

        self.rates = [
            self._as_coefficient(reaction.rate, V0, f"rate_{i}")
            for i, reaction in enumerate(self.gas.reactions)
        ]

    def _in_cells(self, coefficient: Coefficient) -> np.ndarray:
        """The value of a coefficient in every cell."""
        if callable(coefficient):
            return np.asarray(coefficient(self.reduced_field), dtype=float)
        return np.full_like(self.reduced_field, coefficient)

    def _as_coefficient(
        self, coefficient: Coefficient, V0: fem.FunctionSpace, name: str
    ) -> float | fem.Function:
        """A constant stays a float; a callable becomes a cell-wise field."""
        if not callable(coefficient):
            return float(coefficient)
        function = fem.Function(V0, name=name)
        function.x.array[:] = self._in_cells(coefficient)
        return function

    def _create_festim_problem(self) -> None:
        gas, geometry = self.gas, self.geometry
        problem = F.HydrogenTransportProblem()
        problem.mesh = F.Mesh(self.mesh, coordinate_system=geometry.coordinate_system)

        self.species = {
            carrier.name: F.Species(carrier.name) for carrier in gas.carriers
        }
        problem.species = list(self.species.values())

        material = CellwiseDiffusionMaterial(self.diffusivities, name=gas.name)
        volume = F.VolumeSubdomain(id=1, material=material)
        anode = F.SurfaceSubdomain(id=ANODE_ID, locator=geometry.is_anode)
        cathode = F.SurfaceSubdomain(id=CATHODE_ID, locator=geometry.is_cathode)
        problem.subdomains = [volume, anode, cathode]
        problem.temperature = gas.temperature

        problem.sources = [
            F.ParticleSource(
                value=self.pair_production, volume=volume, species=self.species[name]
            )
            for name in gas.ionisation_products
        ]
        problem.reactions = [
            F.GenericReaction(
                reactant=[self.species[name] for name in reaction.reactants],
                product=[self.species[name] for name in reaction.products] or None,
                forward_rate=rate,
                volume=volume,
            )
            for reaction, rate in zip(gas.reactions, self.rates)
        ]
        problem.drift_terms = [
            MobilityElectromigration(
                species=self.species[carrier.name],
                charge=carrier.charge,
                mobility=self.mobilities[carrier.name],
                potential=self.potential,
                subdomain=volume,
            )
            for carrier in gas.carriers
        ]

        # An electrode absorbs the carriers the field drives into it (outflow) and
        # emits none (zero density for the carriers the field pushes away).
        problem.boundary_conditions = []
        for carrier in gas.carriers:
            species = self.species[carrier.name]
            drawn_to_anode = carrier.charge * self.bias < 0
            collector, repeller = (
                (anode, cathode) if drawn_to_anode else (cathode, anode)
            )
            problem.boundary_conditions += [
                F.OutflowBC(subdomain=collector, species=species),
                F.FixedConcentrationBC(subdomain=repeller, value=0.0, species=species),
            ]

        problem.settings = F.Settings(
            atol=1e-30, rtol=1e-10, max_iterations=40, transient=False
        )
        problem.initialise()
        self.festim_problem = problem
