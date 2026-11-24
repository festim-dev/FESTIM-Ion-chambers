"""Transport of the charge carriers in the gas, in a given electric potential."""

from __future__ import annotations

import festim as F
import numpy as np
import ufl
from dolfinx import fem
from scifem import assemble_scalar

from ion_chambers.constants import TOWNSEND
from ion_chambers.electrostatic_problem import ElectrostaticProblem
from ion_chambers.exports import ChamberQuantity
from ion_chambers.gas import Coefficient, Gas
from ion_chambers.source import IonisationSource


class MobilityElectromigration(F.ElectromigrationTerm):
    """Electromigration with the drift velocity ``v = -z mu grad(phi)``.

    :class:`festim.ElectromigrationTerm` derives the mobility from the diffusion
    coefficient with the Einstein relation. That is wrong here for two reasons: the
    diffusion coefficient carries a numerical floor (see :class:`DiffusionFloor`), and
    electrons heated by the field do not follow the Einstein relation at the gas
    temperature.

    Args:
        species: the species the drift acts on
        charge: the charge number of the species
        mobility: the mobility (m2/V/s), a constant or a field
        potential: the electric potential (V)
        subdomain: the volume subdomain where the drift applies
    """

    def __init__(
        self,
        species: F.Species,
        charge: int,
        mobility: float | fem.Function,
        potential: fem.Function,
        subdomain: F.VolumeSubdomain,
    ) -> None:
        super().__init__(
            species=species, charge=charge, potential=potential, subdomain=subdomain
        )
        self.mobility = mobility

    def drift_velocity(self, D, temperature):
        return -self.charge * self.mobility * ufl.grad(self.potential.fenics_object)


class CellwiseDiffusionMaterial(F.Material):
    """Material whose diffusion coefficient is a separate field for each species.

    :class:`festim.Material` accepts a field through ``D``, but uses that one field for
    every species. Here each charge carrier needs its own.

    Args:
        D_fields: the diffusion coefficient of each species (m2/s), keyed by species
            name
        name: the name of the material
    """

    def __init__(self, D_fields: dict[str, fem.Function], name: str | None = None):
        # The parent class requires D_0 and E_D for every species. They are not used:
        # get_diffusion_coefficient below returns the fields instead.
        super().__init__(
            D_0=dict.fromkeys(D_fields, 1.0),
            E_D=dict.fromkeys(D_fields, 0.0),
            name=name,
        )
        self.D_fields = D_fields

    def get_diffusion_coefficient(self, mesh=None, temperature=None, species=None):
        return self.D_fields[species.name]


class DiffusionFloor:
    r"""A cell-local lower bound on the diffusion coefficient.

    FESTIM does not stabilise drift-dominated transport, and the cell Péclet number of
    the ions reaches about 100. Each carrier therefore gets

    .. math::

        D = \max\left(D_\mathrm{phys},\ \tfrac{1}{2}\, \mathrm{factor}\, \mu |E| h_E
        \right)

    with :math:`h_E` the cell length along the field. This is first-order upwinding
    written as diffusion. The floor vanishes where the field is weak, so it cannot
    fake losses in low-field corners, and the error it leaves is first order in the
    cell size. Extrapolate from two meshes for numbers you intend to quote.

    Args:
        factor: the multiplier of the floor. 1 gives a cell Péclet number of 2; 0
            switches the stabilisation off.
    """

    def __init__(self, factor: float = 1.0) -> None:
        if factor < 0:
            raise ValueError(f"factor must not be negative, not {factor}")
        self.factor = factor

    def apply(
        self,
        diffusivity: np.ndarray,
        mobility: np.ndarray,
        field_strength: np.ndarray,
        cell_length: np.ndarray,
    ) -> np.ndarray:
        """Returns the stabilised diffusion coefficient in each cell.

        Args:
            diffusivity: the physical diffusion coefficient in each cell (m2/s)
            mobility: the mobility in each cell (m2/V/s)
            field_strength: ``|E|`` in each cell (V/m)
            cell_length: the cell length along the field (m)
        """
        floor = self.factor * 0.5 * mobility * field_strength * cell_length
        return np.maximum(diffusivity, floor)


def cell_field_length(potential: fem.Function) -> tuple[np.ndarray, np.ndarray]:
    """The field strength in each cell and the length of the cell along the field.

    Only for triangular meshes.

    Args:
        potential: the electric potential (V)

    Returns:
        ``|E|`` at the cell centroids (V/m) and the cell length along E (m)
    """
    mesh = potential.function_space.mesh
    vertices = mesh.geometry.x[:, :2][mesh.geometry.dofmaps[0]]  # (cells, 3, 2)
    centroids = vertices.mean(axis=1)
    n_cells = len(centroids)

    # the electric field E = -grad(potential), evaluated at the cell centroids
    V = fem.functionspace(mesh, ("DG", 1, (2,)))
    field = fem.Function(V)
    field.interpolate(
        fem.Expression(-ufl.grad(potential), V.element.interpolation_points)
    )
    points = np.column_stack([centroids, np.zeros(n_cells)])
    cells = np.arange(n_cells, dtype=np.int32)
    field_at_centroids = field.eval(points, cells)

    strength = np.linalg.norm(field_at_centroids, axis=1)
    direction = field_at_centroids / np.maximum(strength[:, None], 1e-300)

    # the length of a cell along the field: its longest edge projected on the field
    edges = np.stack(
        [
            vertices[:, 1] - vertices[:, 0],
            vertices[:, 2] - vertices[:, 1],
            vertices[:, 0] - vertices[:, 2],
        ],
        axis=1,
    )
    projections = np.abs(np.einsum("cej,cj->ce", edges, direction))
    length = projections.max(axis=1)
    return strength, length


class ChargeTransportProblem(F.HydrogenTransportProblem):
    r"""Steady-state transport of the charge carriers of a gas.

    Each charge carrier follows the Nernst-Planck flux

    .. math::

        J = -D \nabla c + c\, v, \qquad v = -z \mu \nabla \varphi

    in the electric potential :math:`\varphi` of an :class:`ElectrostaticProblem`, is
    created by the ionisation sources and reacts with the other carriers.

    It is a :class:`festim.HydrogenTransportProblem` whose species, reactions, drift
    terms and boundary conditions are set up from the gas and the electrostatic
    problem during :meth:`initialise`. The electrostatic problem must have been run
    before. The problem is solved on the mesh of the electrostatic problem, and its
    electrodes are the surfaces that carry a potential there.

    ``sources``, ``drift_terms`` and ``boundary_conditions`` can be set as on any
    FESTIM problem: :meth:`initialise` adds its own to them. A boundary condition set
    for a carrier on an electrode replaces the automatic one. Sources and drift terms
    are given on the volume subdomain of the electrostatic problem. ``mesh``,
    ``subdomains``, ``species``, ``reactions`` and ``temperature`` cannot be set:
    they come from the electrostatic problem and the gas.

    Calling ``run()`` again after changing the activity of a source reuses everything.
    Any other change needs a new problem.

    Args:
        electrostatic_problem: the problem that provides the electric potential, the
            mesh and the electrodes
        gas: the gas filling the chamber
        ionisation_sources: the ionisation sources
        exports: the quantities to compute after each run: the classes of
            :mod:`ion_chambers.exports` and FESTIM's own exports
        stabilisation: the stabilisation of the drift. Defaults to
            :class:`DiffusionFloor`.
        settings: the solver settings. The default suits this problem.

    Attributes:
        electrostatic_problem: the problem that provides the electric potential
        gas: the gas filling the chamber
        ionisation_sources: the ionisation sources
        stabilisation: the stabilisation of the drift
        gas_volume: the volume subdomain of the gas, after ``initialise``. A
            :class:`festim.VolumeSubdomain` with the id of the volume subdomain of
            the electrostatic problem.
        potential: the electric potential (V), after ``initialise``
        reduced_field: E/N in each cell (Td), after ``initialise``
        mobilities: the mobility of each carrier (m2/V/s), keyed by name. A float, or
            a cell-wise field when it depends on E/N.
        diffusivities: the stabilised diffusion coefficient of each carrier (m2/s),
            keyed by name, a cell-wise field
        rates: the rate coefficient of each reaction, in the order of
            ``gas.reactions``
        pair_production: the total production rate of ion pairs (1/m3/s)
        volume_of_gas: the volume of gas in the chamber (m3)

    Examples:
        .. code-block:: python

            current = ic.ChamberCurrent(surface=chamber.anode)

            transport = ic.ChargeTransportProblem()
            transport.electrostatic_problem = electrostatics  # already run
            transport.gas = gas
            transport.ionisation_sources = [ic.TritiumInGas(activity=3.7e10)]
            transport.exports = [current]
            transport.initialise()
            transport.run()
            print(current.value)
    """

    def __init__(
        self,
        electrostatic_problem: ElectrostaticProblem | None = None,
        gas: Gas | None = None,
        ionisation_sources: list[IonisationSource] | None = None,
        exports: list | None = None,
        stabilisation: DiffusionFloor | None = None,
        settings: F.Settings | None = None,
    ) -> None:
        if settings is None:
            settings = F.Settings(
                atol=1e-30, rtol=1e-10, max_iterations=40, transient=False
            )
        super().__init__(exports=exports, settings=settings)

        self.electrostatic_problem = electrostatic_problem
        self.gas = gas
        self.ionisation_sources = ionisation_sources or []
        self.stabilisation = stabilisation or DiffusionFloor()

        # Created during initialise, once the electric field and therefore the
        # diffusion coefficients are known.
        self.gas_volume = None
        self.potential = None

    def initialise(self) -> None:
        """Sets up the FESTIM problem from the gas and the electrostatic problem."""
        if self.electrostatic_problem is None or self.electrostatic_problem.u is None:
            raise RuntimeError(
                "the electrostatic problem must be initialised and run before the "
                "transport problem is initialised"
            )
        self._check_inputs_taken_from_elsewhere()
        self.mesh = self.electrostatic_problem.mesh
        self.potential = self.electrostatic_problem.u

        cell_name = self.mesh.mesh.topology.cell_name()
        if cell_name != "triangle":
            raise ValueError(f"the mesh must be made of triangles, not {cell_name}s")

        # The same cell and facet tags as the electrostatic problem, whether it made
        # them from the locators of its subdomains or was given them with the mesh.
        self.volume_meshtags = self.electrostatic_problem.volume_meshtags
        self.facet_meshtags = self.electrostatic_problem.facet_meshtags

        self._create_coefficients()
        self._create_subdomains()
        self._create_species_and_reactions()
        self._create_sources()
        self._create_drift_terms()
        self._create_boundary_conditions()
        self.temperature = self.gas.temperature

        super().initialise()

        weight = self.mesh.coordinate_system.integration_weight(self.mesh.mesh)
        self.volume_of_gas = assemble_scalar(weight * self.dx(self.gas_volume.id))

    def _check_inputs_taken_from_elsewhere(self) -> None:
        """Raises if an attribute that ``initialise`` fills was set by the user.

        The mesh and the subdomains come from the electrostatic problem, and the
        species, the reactions and the temperature from the gas. A value set on the
        transport problem would be replaced, so it is refused instead.
        """
        taken_from = {
            "mesh": "electrostatic_problem",
            "subdomains": "electrostatic_problem",
            "species": "the carriers of gas",
            "reactions": "gas",
            "temperature": "gas",
        }
        for name, origin in taken_from.items():
            value = getattr(self, name)
            is_unset = value is None or (isinstance(value, list) and len(value) == 0)
            if not is_unset:
                raise ValueError(
                    f"{name} must not be set on a ChargeTransportProblem: it is "
                    f"taken from {origin}"
                )

    def _create_coefficients(self) -> None:
        """Evaluates the mobilities, diffusivities and rates in the local field."""
        field_strength, cell_length = cell_field_length(self.potential)
        self.reduced_field = field_strength / self.gas.number_density / TOWNSEND
        V0 = fem.functionspace(self.mesh.mesh, ("DG", 0))

        self.mobilities = {}
        self.diffusivities = {}
        for carrier in self.gas.carriers:
            mobility_in_cells = self._value_in_cells(carrier.mobility)
            if carrier.diffusivity is None:
                # Einstein relation
                diffusivity_in_cells = mobility_in_cells * self.gas.thermal_voltage
            else:
                diffusivity_in_cells = self._value_in_cells(carrier.diffusivity)

            diffusivity = fem.Function(V0, name=f"D_{carrier.name}")
            diffusivity.x.array[:] = self.stabilisation.apply(
                diffusivity_in_cells, mobility_in_cells, field_strength, cell_length
            )
            self.diffusivities[carrier.name] = diffusivity
            self.mobilities[carrier.name] = self._as_fenics_coefficient(
                carrier.mobility, V0, name=f"mobility_{carrier.name}"
            )

        self.rates = []
        for i, reaction in enumerate(self.gas.reactions):
            rate = self._as_fenics_coefficient(reaction.rate, V0, name=f"rate_{i}")
            self.rates.append(rate)

    def _value_in_cells(self, coefficient: Coefficient) -> np.ndarray:
        """The value of a coefficient in every cell, as an array."""
        if callable(coefficient):
            return np.asarray(coefficient(self.reduced_field), dtype=float)
        return np.full_like(self.reduced_field, coefficient)

    def _as_fenics_coefficient(
        self, coefficient: Coefficient, V0: fem.FunctionSpace, name: str
    ) -> float | fem.Function:
        """A constant stays a float; a callable of E/N becomes a cell-wise field."""
        if not callable(coefficient):
            return float(coefficient)
        function = fem.Function(V0, name=name)
        function.x.array[:] = self._value_in_cells(coefficient)
        return function

    def _create_subdomains(self) -> None:
        """The gas volume, plus the surfaces of the electrostatic problem.

        The gas volume is a new subdomain with the id of the one of the electrostatic
        problem, because it carries the material of this problem.
        """
        volumes = self.electrostatic_problem.volume_subdomains
        if len(volumes) != 1:
            raise ValueError(
                f"the electrostatic problem must have one volume subdomain, the gas, "
                f"not {len(volumes)}"
            )
        material = CellwiseDiffusionMaterial(
            D_fields=self.diffusivities, name=self.gas.name
        )
        self.gas_volume = F.VolumeSubdomain(id=volumes[0].id, material=material)
        self.subdomains = [self.gas_volume]
        for surface in self.electrostatic_problem.surface_subdomains:
            self.subdomains.append(surface)

    def _create_species_and_reactions(self) -> None:
        # A ChargeCarrier is a festim.Species and a CarrierReaction is a
        # festim.GenericReaction, so FESTIM takes them as they are. The reactions
        # were created without a volume or a final rate: both are known now.
        self.species = self.gas.carriers
        for i, reaction in enumerate(self.gas.reactions):
            reaction.volume = self.gas_volume
            reaction.forward_rate = self.rates[i]
        self.reactions = self.gas.reactions

    def _create_sources(self) -> None:
        """Every ionisation creates one of each of the two ionisation products.

        The sources are added to those the user set in ``sources``. A source the
        user gave on the volume of the electrostatic problem is moved to the gas
        volume of this problem, which has the same id and carries the material.
        """
        electrostatic_volume = self.electrostatic_problem.volume_subdomains[0]
        for source in self.sources:
            if source.volume is electrostatic_volume:
                source.volume = self.gas_volume

        self.pair_production = 0
        for source in self.ionisation_sources:
            self.pair_production += source.create_pair_production(self.mesh, self.gas)

        for carrier in self.gas.ionisation_products:
            particle_source = F.ParticleSource(
                value=self.pair_production, volume=self.gas_volume, species=carrier
            )
            self.sources.append(particle_source)

    def _create_drift_terms(self) -> None:
        """The drift of every carrier in the electric field.

        The terms are added to those the user set in ``drift_terms``. A term the
        user gave on the volume of the electrostatic problem is moved to the gas
        volume of this problem, which has the same id and carries the material.
        """
        electrostatic_volume = self.electrostatic_problem.volume_subdomains[0]
        for drift_term in self.drift_terms:
            if drift_term.subdomain is electrostatic_volume:
                drift_term.subdomain = self.gas_volume

        for carrier in self.gas.carriers:
            drift = MobilityElectromigration(
                species=carrier,
                charge=carrier.charge,
                mobility=self.mobilities[carrier.name],
                potential=self.potential,
                subdomain=self.gas_volume,
            )
            self.drift_terms.append(drift)

    def _create_boundary_conditions(self) -> None:
        """An electrode absorbs the carriers the field drives into it (outflow) and
        emits none (zero density for the carriers the field pushes away).

        Negative carriers are driven to the electrodes above the mean electrode
        potential, positive carriers to those below it.

        The conditions are added to those the user set in ``boundary_conditions``.
        Where the user gave a condition for a carrier on an electrode, that one is
        kept and none is created.
        """
        set_by_user = []
        for bc in self.boundary_conditions:
            set_by_user.append((bc.subdomain.id, getattr(bc, "species", None)))

        electrode_bcs = self.electrostatic_problem.boundary_conditions
        total = 0.0
        for bc in electrode_bcs:
            total += bc.value
        mean_potential = total / len(electrode_bcs)

        for bc in electrode_bcs:
            electrode = bc.subdomain
            electrode_is_positive = bc.value > mean_potential
            for carrier in self.gas.carriers:
                if (electrode.id, carrier) in set_by_user:
                    continue
                carrier_is_negative = carrier.charge < 0
                attracted = carrier_is_negative == electrode_is_positive
                if attracted:
                    carrier_bc = F.OutflowBC(subdomain=electrode, species=carrier)
                else:
                    carrier_bc = F.FixedConcentrationBC(
                        subdomain=electrode, value=0.0, species=carrier
                    )
                self.boundary_conditions.append(carrier_bc)

    def initialise_exports(self) -> None:
        """Initialises FESTIM's exports, leaving the chamber quantities aside.

        FESTIM refuses its derived quantities on a cylindrical mesh, because most of
        them do not integrate with the factor 2 pi r. The chamber quantities do, so
        they are hidden from FESTIM's check here and computed in post_processing.
        """
        all_exports = self.exports

        festim_exports = []
        for export in all_exports:
            if isinstance(export, ChamberQuantity):
                # start from empty data, as FESTIM does for its own quantities
                export.t = []
                export.data = []
            else:
                festim_exports.append(export)

        self.exports = festim_exports
        super().initialise_exports()
        self.exports = all_exports

    def post_processing(self) -> None:
        """Computes FESTIM's exports, then the chamber quantities."""
        super().post_processing()
        for export in self.exports:
            if isinstance(export, ChamberQuantity):
                export.compute(self)
                export.t.append(float(self.t))
                export.write(t=float(self.t))
