"""Quantities computed from a solved transport problem.

They are used like FESTIM's derived quantities: create them with the surface or the
volume they are computed on, add them to the ``exports`` of a
:class:`ion_chambers.ChargeTransportProblem`, and read ``value`` after ``run()``.
``data`` keeps the value of every run, and a ``filename`` writes them to a csv file.

FESTIM's quantities are computed from the field of one species. The quantities here
involve every charge carrier and the electric potential, so they take no field and
their ``compute`` takes the whole problem.

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


def _charge_flux(problem: ChargeTransportProblem):
    """The net flux of charge, in units of the elementary charge (1/m2/s).

    The sum over the carriers of the charge number times the Nernst-Planck flux
    ``-D grad(c) + c v``, with the drift velocity ``v = -z mu grad(potential)``.
    """
    total = 0
    for carrier in problem.gas.carriers:
        density = carrier.post_processing_solution
        mobility = problem.mobilities[carrier.name]
        diffusivity = problem.diffusivities[carrier.name]
        velocity = -carrier.charge * mobility * ufl.grad(problem.potential)
        flux = -diffusivity * ufl.grad(density) + density * velocity
        total = total + carrier.charge * flux
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
        """Computes the quantity for a solved problem and stores it in ``value`` and
        ``data``. Implemented by each subclass.

        Args:
            problem: the transport problem, after a solve
        """
        raise NotImplementedError


class ChamberSurfaceQuantity(ChamberQuantity):
    """Base class of the quantities computed on a surface, the counterpart of
    :class:`festim.SurfaceQuantity`.

    Args:
        surface: the surface subdomain
        filename: the name of a csv or txt file to write the values to, optional

    Attributes:
        surface: the surface subdomain
    """

    def __init__(
        self, surface: F.SurfaceSubdomain, filename: str | None = None
    ) -> None:
        super().__init__(filename=filename)
        if not isinstance(surface, F.SurfaceSubdomain):
            raise TypeError(
                f"surface must be a festim.SurfaceSubdomain, "
                f"not {type(surface).__name__}"
            )
        self.surface = surface


class ChamberVolumeQuantity(ChamberQuantity):
    """Base class of the quantities computed in a volume, the counterpart of
    :class:`festim.VolumeQuantity`.

    Args:
        volume: the volume subdomain
        filename: the name of a csv or txt file to write the values to, optional

    Attributes:
        volume: the volume subdomain
    """

    def __init__(self, volume: F.VolumeSubdomain, filename: str | None = None) -> None:
        super().__init__(filename=filename)
        if not isinstance(volume, F.VolumeSubdomain):
            raise TypeError(
                f"volume must be a festim.VolumeSubdomain, not {type(volume).__name__}"
            )
        self.volume = volume


class ChamberCurrent(ChamberSurfaceQuantity):
    r"""The current the chamber reads (A), from the Shockley-Ramo theorem.

    .. math::

        I = e \int \sum_s z_s J_s \cdot (-\nabla \varphi_w)\, dV

    where the weighting potential :math:`\varphi_w` is 1 on the readout electrode and
    0 on the other one. In steady state this equals the current through either
    electrode. Being an integral over the gas, it does not suffer from the error of
    the gradient at the boundary, which is why it is the recommended reading.

    The sign follows the potential of the readout electrode: positive when it is the
    more positive electrode.

    Args:
        surface: the surface of the readout electrode
        filename: the name of a csv or txt file to write the values to, optional
    """

    @property
    def title(self) -> str:
        return "Chamber current (A)"

    def compute(self, problem: ChargeTransportProblem) -> None:
        # The potential difference between the readout electrode and the other one.
        # Dividing the potential by it gives the weighting potential.
        electrode_bcs = problem.electrostatic_problem.boundary_conditions
        if len(electrode_bcs) != 2:
            raise ValueError(
                f"ChamberCurrent needs two electrodes, not {len(electrode_bcs)}: "
                f"use ElectrodeCurrent for the current through one electrode"
            )
        readout_potential = None
        other_potential = None
        for bc in electrode_bcs:
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

        weight = problem.mesh.coordinate_system.integration_weight(problem.mesh.mesh)
        dx = problem.dx(problem.gas_volume.id)
        integrand = ufl.dot(_charge_flux(problem), weighting_field)
        self.value = E_CHARGE * assemble_scalar(weight * integrand * dx)
        self.data.append(self.value)


class ElectrodeCurrent(ChamberSurfaceQuantity):
    """The charge flowing out of the gas through an electrode per unit time (A).

    Computed from the carrier fluxes at the electrode. It is less accurate than
    :class:`ChamberCurrent` on a given mesh, and comparing the two is a check of the
    discretisation. For the flux of a single carrier, use
    ``festim.SurfaceFlux(field=carrier, surface=surface)``.

    It is negative on the electrode that collects the negative carriers.

    Args:
        surface: the surface of the electrode
        filename: the name of a csv or txt file to write the values to, optional
    """

    @property
    def title(self) -> str:
        return f"Current through surface {self.surface.id} (A)"

    def compute(self, problem: ChargeTransportProblem) -> None:
        weight = problem.mesh.coordinate_system.integration_weight(problem.mesh.mesh)
        ds = problem.ds(self.surface.id)
        normal = ufl.FacetNormal(problem.mesh.mesh)
        outward_flux = ufl.dot(_charge_flux(problem), normal)
        self.value = E_CHARGE * assemble_scalar(weight * outward_flux * ds)
        self.data.append(self.value)


class PairGeneration(ChamberVolumeQuantity):
    """The production rate of ion pairs in a volume (1/s).

    Args:
        volume: the volume subdomain of the gas
        filename: the name of a csv or txt file to write the values to, optional
    """

    @property
    def title(self) -> str:
        return f"Pair generation volume {self.volume.id} (1/s)"

    def compute(self, problem: ChargeTransportProblem) -> None:
        weight = problem.mesh.coordinate_system.integration_weight(problem.mesh.mesh)
        dx = problem.dx(self.volume.id)
        self.value = assemble_scalar(weight * problem.pair_production * dx)
        self.data.append(self.value)


class Recombination(ChamberVolumeQuantity):
    """The loss rate of ion pairs to recombination in a volume (1/s).

    Args:
        volume: the volume subdomain of the gas
        filename: the name of a csv or txt file to write the values to, optional
    """

    @property
    def title(self) -> str:
        return f"Recombination volume {self.volume.id} (1/s)"

    def compute(self, problem: ChargeTransportProblem) -> None:
        weight = problem.mesh.coordinate_system.integration_weight(problem.mesh.mesh)
        dx = problem.dx(self.volume.id)
        total = 0.0
        for i, reaction in enumerate(problem.gas.reactions):
            if not reaction.is_recombination:
                continue
            # mass action: the rate coefficient times the density of every reactant
            local_rate = problem.rates[i]
            for carrier in reaction.reactant:
                local_rate = local_rate * carrier.post_processing_solution
            total += assemble_scalar(weight * local_rate * dx)
        self.value = total
        self.data.append(self.value)


class IdealCurrent(ChamberVolumeQuantity):
    """The current if every beta stopped in the gas and every charge were collected
    (A).

    Args:
        volume: the volume subdomain of the gas
        filename: the name of a csv or txt file to write the values to, optional
    """

    @property
    def title(self) -> str:
        return f"Ideal current volume {self.volume.id} (A)"

    def compute(self, problem: ChargeTransportProblem) -> None:
        weight = problem.mesh.coordinate_system.integration_weight(problem.mesh.mesh)
        dx = problem.dx(self.volume.id)
        size_of_volume = assemble_scalar(weight * dx)

        ideal_pair_rate = 0.0
        for source in problem.ionisation_sources:
            ideal_pair_rate += source.ideal_pair_rate(problem.gas, size_of_volume)
        self.value = E_CHARGE * ideal_pair_rate
        self.data.append(self.value)


class WallEfficiency(ChamberVolumeQuantity):
    """The share of the beta energy deposited in the gas rather than in the walls.

    The current carried by the pairs that are generated, divided by the ideal current.

    Args:
        volume: the volume subdomain of the gas
        filename: the name of a csv or txt file to write the values to, optional
    """

    @property
    def title(self) -> str:
        return f"Wall efficiency volume {self.volume.id}"

    def compute(self, problem: ChargeTransportProblem) -> None:
        generation = PairGeneration(volume=self.volume)
        generation.compute(problem)
        ideal_current = IdealCurrent(volume=self.volume)
        ideal_current.compute(problem)

        if ideal_current.value == 0:
            self.value = float("nan")
        else:
            self.value = E_CHARGE * generation.value / ideal_current.value
        self.data.append(self.value)


class CollectionEfficiency(ChamberVolumeQuantity):
    """The share of the ion pairs that escape recombination.

    Args:
        volume: the volume subdomain of the gas
        filename: the name of a csv or txt file to write the values to, optional
    """

    @property
    def title(self) -> str:
        return f"Collection efficiency volume {self.volume.id}"

    def compute(self, problem: ChargeTransportProblem) -> None:
        generation = PairGeneration(volume=self.volume)
        generation.compute(problem)
        recombination = Recombination(volume=self.volume)
        recombination.compute(problem)

        if generation.value == 0:
            self.value = float("nan")
        else:
            self.value = (generation.value - recombination.value) / generation.value
        self.data.append(self.value)


class SpaceChargeField(ChamberVolumeQuantity):
    r"""The field the charge in the gas would create, relative to the applied field.

    The transport problem uses a fixed potential and ignores the charge in the gas.
    This quantity tells how far that holds. It solves a second
    :class:`ion_chambers.ElectrostaticProblem` with the converged charge density and
    grounded electrodes, and returns

    .. math::

        \frac{\lVert \nabla \varphi_1 \rVert}{\lVert \nabla \varphi \rVert}

    with :math:`\varphi_1` the potential of the space charge and :math:`\varphi` the
    applied potential, in the L2 norm over the volume. Treat results as indicative
    once it exceeds about 0.1.

    Args:
        volume: the volume subdomain of the gas
        filename: the name of a csv or txt file to write the values to, optional
    """

    def __init__(self, volume: F.VolumeSubdomain, filename: str | None = None) -> None:
        super().__init__(volume=volume, filename=filename)
        self._space_charge_problem = None
        self._charge_density = None

    @property
    def title(self) -> str:
        return f"Space-charge field / applied field volume {self.volume.id}"

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
        # the same cell and facet tags as the applied problem
        self._space_charge_problem.volume_meshtags = applied.volume_meshtags
        self._space_charge_problem.facet_meshtags = applied.facet_meshtags
        self._space_charge_problem.initialise()

    def compute(self, problem: ChargeTransportProblem) -> None:
        if self._space_charge_problem is None:
            self._create_space_charge_problem(problem)

        # update the charge density with the latest solution, then solve
        net_charge_density = 0
        for carrier in problem.gas.carriers:
            net_charge_density += carrier.charge * carrier.post_processing_solution
        expression = fem.Expression(
            E_CHARGE * net_charge_density,
            self._charge_density.function_space.element.interpolation_points,
        )
        self._charge_density.interpolate(expression)
        self._space_charge_problem.run()

        induced_field = ufl.grad(self._space_charge_problem.u)
        applied_field = ufl.grad(problem.potential)
        weight = problem.mesh.coordinate_system.integration_weight(problem.mesh.mesh)
        dx = problem.dx(self.volume.id)
        induced_norm_squared = assemble_scalar(
            weight * ufl.dot(induced_field, induced_field) * dx
        )
        applied_norm_squared = assemble_scalar(
            weight * ufl.dot(applied_field, applied_field) * dx
        )
        self.value = float(np.sqrt(induced_norm_squared / applied_norm_squared))
        self.data.append(self.value)
