"""Quantities computed from a solved transport problem.

They are used like FESTIM's derived quantities: create them, add them to the
``exports`` of a :class:`ion_chambers.ChargeTransportProblem`, and read ``value``
after ``run()``. ``data`` keeps the value of every run, and a ``filename`` writes
them to a csv file.

FESTIM's field exports work too, for instance ``festim.SpeciesExport`` to write the
density of a carrier to a file. Of FESTIM's derived quantities, only
``festim.SurfaceFlux`` supports the cylindrical meshes used here.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import festim as F
import numpy as np
import ufl
from dolfinx import fem
from scifem import assemble_scalar

from ion_chambers.constants import E_CHARGE
from ion_chambers.electrostatic_problem import ElectrostaticProblem, FixedPotentialBC

if TYPE_CHECKING:
    from ion_chambers.charge_transport_problem import ChargeTransportProblem
    from ion_chambers.gas import ChargeCarrier


def integration_weight(problem: ChargeTransportProblem):
    """The factor 2 pi r of a cylindrical mesh, or 1 on a cartesian mesh."""
    return problem.mesh.coordinate_system.integration_weight(problem.mesh.mesh)


def carrier_flux(problem: ChargeTransportProblem, carrier: ChargeCarrier):
    """The Nernst-Planck number flux of a carrier (1/m2/s), a ufl expression."""
    density = carrier.post_processing_solution
    mobility = problem.mobilities[carrier.name]
    diffusivity = problem.diffusivities[carrier.name]
    velocity = -carrier.charge * mobility * ufl.grad(problem.potential)
    return -diffusivity * ufl.grad(density) + density * velocity


def charge_flux(problem: ChargeTransportProblem):
    """The net flux of charge, in units of the elementary charge (1/m2/s)."""
    total = 0
    for carrier in problem.gas.carriers:
        total = total + carrier.charge * carrier_flux(problem, carrier)
    return total


def net_charge_density(problem: ChargeTransportProblem):
    """The net charge number density (1/m3), a ufl expression."""
    total = 0
    for carrier in problem.gas.carriers:
        total = total + carrier.charge * carrier.post_processing_solution
    return total


def generation_rate(problem: ChargeTransportProblem) -> float:
    """The production rate of ion pairs in the whole gas volume (1/s)."""
    weight = integration_weight(problem)
    dx = problem.dx(problem.gas_volume.id)
    return assemble_scalar(weight * problem.pair_production * dx)


def recombination_rate(problem: ChargeTransportProblem) -> float:
    """The loss rate of ion pairs to recombination in the whole gas volume (1/s)."""
    weight = integration_weight(problem)
    dx = problem.dx(problem.gas_volume.id)
    total = 0.0
    for i, reaction in enumerate(problem.gas.reactions):
        if not reaction.is_recombination:
            continue
        # mass action: the rate coefficient times the density of every reactant
        local_rate = problem.rates[i]
        for carrier in reaction.reactant:
            local_rate = local_rate * carrier.post_processing_solution
        total += assemble_scalar(weight * local_rate * dx)
    return total


def ideal_pair_rate(problem: ChargeTransportProblem) -> float:
    """The pair production rate (1/s) if no beta energy were lost to the walls."""
    total = 0.0
    for source in problem.ionisation_sources:
        total += source.ideal_pair_rate(problem.gas, problem.volume_of_gas)
    return total


class ChamberQuantity(F.exports.DerivedQuantity):
    """Base class of the quantities computed from a transport problem.

    A :class:`festim.exports.DerivedQuantity` whose ``compute`` takes the whole
    problem, because these quantities involve every carrier and the potential.

    Args:
        filename: the name of a csv or txt file to write the values to, optional

    Attributes:
        value: the value after the last run
        data: the values of all the runs
        t: the time of each run (always 0 for a steady-state problem)
    """

    def compute(self, problem: ChargeTransportProblem) -> None:
        """Computes the quantity and stores it in ``value`` and ``data``."""
        self.value = float(self.evaluate(problem))
        self.data.append(self.value)

    def evaluate(self, problem: ChargeTransportProblem) -> float:
        """Returns the quantity for a solved problem. Implemented by each subclass."""
        raise NotImplementedError


class ChamberCurrent(ChamberQuantity):
    r"""The current the chamber reads (A), from the Shockley-Ramo theorem.

    .. math::

        I = e \int \sum_s z_s J_s \cdot (-\nabla \varphi_w)\, dV

    where the weighting potential :math:`\varphi_w` is 1 on the readout electrode and
    0 on the other one. In steady state this equals the current through either
    electrode. Being a volume integral, it does not suffer from the error of the
    gradient at the boundary, which is why it is the recommended reading.

    The sign follows the potential of the readout electrode: positive when it is the
    more positive electrode.

    Args:
        surface: the surface of the readout electrode
        filename: the name of a csv or txt file to write the values to, optional
    """

    def __init__(
        self, surface: F.SurfaceSubdomain, filename: str | None = None
    ) -> None:
        super().__init__(filename=filename)
        self.surface = surface

    @property
    def title(self) -> str:
        return "Chamber current (A)"

    def evaluate(self, problem: ChargeTransportProblem) -> float:
        # The potential difference between the readout electrode and the other one.
        # Dividing the potential by it gives the weighting potential.
        readout_potential = None
        other_potential = None
        for bc in problem.electrostatic_problem.boundary_conditions:
            if bc.subdomain is self.surface:
                readout_potential = bc.value
            else:
                other_potential = bc.value
        if readout_potential is None:
            raise ValueError(
                "the surface of ChamberCurrent is not an electrode of the "
                "electrostatic problem"
            )
        bias = readout_potential - other_potential

        weighting_field = -ufl.grad(problem.potential) / bias
        integrand = ufl.dot(charge_flux(problem), weighting_field)
        weight = integration_weight(problem)
        dx = problem.dx(problem.gas_volume.id)
        return E_CHARGE * assemble_scalar(weight * integrand * dx)


class ElectrodeCurrent(ChamberQuantity):
    """The charge flowing out of the gas through an electrode per unit time (A).

    Computed from the carrier fluxes at the electrode. It is less accurate than
    :class:`ChamberCurrent` on a given mesh, and comparing the two is a check of the
    discretisation.

    It is negative on the electrode that collects the negative carriers.

    Args:
        surface: the surface of the electrode
        filename: the name of a csv or txt file to write the values to, optional
    """

    def __init__(
        self, surface: F.SurfaceSubdomain, filename: str | None = None
    ) -> None:
        super().__init__(filename=filename)
        self.surface = surface

    @property
    def title(self) -> str:
        return f"Current through surface {self.surface.id} (A)"

    def evaluate(self, problem: ChargeTransportProblem) -> float:
        normal = ufl.FacetNormal(problem.mesh.mesh)
        weight = integration_weight(problem)
        ds = problem.ds(self.surface.id)
        outward_flux = ufl.dot(charge_flux(problem), normal)
        return E_CHARGE * assemble_scalar(weight * outward_flux * ds)


class IdealCurrent(ChamberQuantity):
    """The current if every beta stopped in the gas and every charge were collected
    (A)."""

    @property
    def title(self) -> str:
        return "Ideal current (A)"

    def evaluate(self, problem: ChargeTransportProblem) -> float:
        return E_CHARGE * ideal_pair_rate(problem)


class PairGeneration(ChamberQuantity):
    """The production rate of ion pairs in the gas (1/s)."""

    @property
    def title(self) -> str:
        return "Pair generation (1/s)"

    def evaluate(self, problem: ChargeTransportProblem) -> float:
        return generation_rate(problem)


class Recombination(ChamberQuantity):
    """The loss rate of ion pairs to recombination in the gas (1/s)."""

    @property
    def title(self) -> str:
        return "Recombination (1/s)"

    def evaluate(self, problem: ChargeTransportProblem) -> float:
        return recombination_rate(problem)


class WallEfficiency(ChamberQuantity):
    """The share of the beta energy deposited in the gas rather than in the walls."""

    @property
    def title(self) -> str:
        return "Wall efficiency"

    def evaluate(self, problem: ChargeTransportProblem) -> float:
        ideal = ideal_pair_rate(problem)
        if ideal == 0:
            return float("nan")
        return generation_rate(problem) / ideal


class CollectionEfficiency(ChamberQuantity):
    """The share of the ion pairs that escape recombination."""

    @property
    def title(self) -> str:
        return "Collection efficiency"

    def evaluate(self, problem: ChargeTransportProblem) -> float:
        generation = generation_rate(problem)
        if generation == 0:
            return float("nan")
        return (generation - recombination_rate(problem)) / generation


class SpaceChargeField(ChamberQuantity):
    r"""The field the charge in the gas would create, relative to the applied field.

    The transport problem uses a fixed potential and ignores the charge in the gas.
    This quantity tells how far that holds. It solves a second
    :class:`ion_chambers.ElectrostaticProblem` with the converged charge density and
    grounded electrodes, and returns

    .. math::

        \frac{\lVert \nabla \varphi_1 \rVert}{\lVert \nabla \varphi \rVert}

    with :math:`\varphi_1` the potential of the space charge and :math:`\varphi` the
    applied potential, in the L2 norm. Treat results as indicative once it exceeds
    about 0.1.
    """

    def __init__(self, filename: str | None = None) -> None:
        super().__init__(filename=filename)
        self._space_charge_problem = None
        self._charge_density = None

    @property
    def title(self) -> str:
        return "Space-charge field / applied field"

    def _create_space_charge_problem(self, problem: ChargeTransportProblem) -> None:
        """The electrostatic problem of the space charge alone: the same electrodes,
        all grounded, with the charge in the gas as source."""
        applied = problem.electrostatic_problem
        V1 = fem.functionspace(problem.mesh.mesh, ("Lagrange", 1))
        self._charge_density = fem.Function(V1, name="charge_density")

        grounded_electrodes = []
        for bc in applied.boundary_conditions:
            grounded_electrodes.append(
                FixedPotentialBC(subdomain=bc.subdomain, value=0.0)
            )

        self._space_charge_problem = ElectrostaticProblem(
            mesh=applied.mesh,
            subdomains=applied.subdomains,
            boundary_conditions=grounded_electrodes,
            charge_density=self._charge_density,
            element_degree=applied.element_degree,
        )
        self._space_charge_problem.initialise()

    def evaluate(self, problem: ChargeTransportProblem) -> float:
        if self._space_charge_problem is None:
            self._create_space_charge_problem(problem)

        # update the charge density with the latest solution, then solve
        expression = fem.Expression(
            E_CHARGE * net_charge_density(problem),
            self._charge_density.function_space.element.interpolation_points,
        )
        self._charge_density.interpolate(expression)
        self._space_charge_problem.run()

        induced_field = ufl.grad(self._space_charge_problem.u)
        applied_field = ufl.grad(problem.potential)
        weight = integration_weight(problem)
        dx = problem.dx(problem.gas_volume.id)
        induced_norm_squared = assemble_scalar(
            weight * ufl.dot(induced_field, induced_field) * dx
        )
        applied_norm_squared = assemble_scalar(
            weight * ufl.dot(applied_field, applied_field) * dx
        )
        return np.sqrt(induced_norm_squared / applied_norm_squared)
