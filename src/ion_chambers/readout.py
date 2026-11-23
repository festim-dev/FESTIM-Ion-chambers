"""What the chamber reads, and the checks that say how far to trust it."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import TYPE_CHECKING

from mpi4py import MPI

import numpy as np
import ufl

from ion_chambers.constants import E_CHARGE
from ion_chambers.integration import assemble, volume_weight

if TYPE_CHECKING:
    from ion_chambers.charge_carrier import ChargeCarrier
    from ion_chambers.model import IonChamberModel


@dataclass(frozen=True)
class Reading:
    """The result of one steady-state solve.

    ``current`` is the reading. ``anode_current``, ``cathode_current`` and
    ``balance_current`` are independent estimates of the same quantity; how closely
    they agree with it measures the discretisation error.

    Attributes:
        current: the chamber current from the Shockley-Ramo theorem (A)
        anode_current: the current from the carrier fluxes into the anode (A)
        cathode_current: the current from the carrier fluxes into the cathode (A)
        balance_current: the current from the volume balance, generation minus
            recombination (A)
        ideal_current: the current if every beta stopped in the gas and every charge
            were collected (A)
        generation: the production rate of ion pairs in the gas (1/s)
        recombination: the loss rate of ion pairs to recombination (1/s)
        wall_efficiency: the share of the beta energy deposited in the gas
        collection_efficiency: the share of the ion pairs that escape recombination
        space_charge_field: the field of the space charge relative to the applied
            field (ratio of the L2 norms). The fixed-field assumption holds while
            this is small.
        space_charge_potential: the largest potential of the space charge relative to
            the bias
        mean_net_charge_density: the volume average of the absolute net charge
            number density (1/m3)
    """

    current: float
    anode_current: float
    cathode_current: float
    balance_current: float
    ideal_current: float
    generation: float
    recombination: float
    wall_efficiency: float
    collection_efficiency: float
    space_charge_field: float
    space_charge_potential: float
    mean_net_charge_density: float

    def as_dict(self) -> dict[str, float]:
        """The reading as a dictionary, for saving."""
        return asdict(self)


def carrier_flux(model: IonChamberModel, carrier: ChargeCarrier):
    """The Nernst-Planck number flux of a carrier (1/m2/s), a ufl expression."""
    density = carrier.post_processing_solution
    mobility = model.mobilities[carrier.name]
    diffusivity = model.diffusivities[carrier.name]
    velocity = -carrier.charge * mobility * ufl.grad(model.potential)
    return -diffusivity * ufl.grad(density) + density * velocity


def charge_flux(model: IonChamberModel):
    """The net flux of charge, in units of the elementary charge (1/m2/s)."""
    total = 0
    for carrier in model.gas.carriers:
        total = total + carrier.charge * carrier_flux(model, carrier)
    return total


def net_charge_density(model: IonChamberModel):
    """The net charge number density (1/m3), a ufl expression."""
    total = 0
    for carrier in model.gas.carriers:
        total = total + carrier.charge * carrier.post_processing_solution
    return total


def _weight(model: IonChamberModel):
    return volume_weight(model.mesh, model.geometry.coordinate_system)


def ramo_current(model: IonChamberModel) -> float:
    r"""The chamber current from the Shockley-Ramo theorem (A).

    .. math::

        I = e \int \sum_s q_s J_s \cdot (-\nabla \varphi_w)\, dV,
        \qquad \varphi_w = \varphi / V

    In steady state this equals the current through either electrode. Being a volume
    integral, it does not suffer from the error of the gradient at the boundary, which
    is why it is used as the reading.
    """
    weighting_field = -ufl.grad(model.potential) / model.bias
    integrand = ufl.dot(charge_flux(model), weighting_field)
    return E_CHARGE * assemble(_weight(model) * integrand * ufl.dx(domain=model.mesh))


def electrode_current(model: IonChamberModel, surface_id: int) -> float:
    """The charge flowing out of the gas through an electrode per unit time (A).

    Args:
        model: the model, after a run
        surface_id: the id of the electrode surface
    """
    n = ufl.FacetNormal(model.mesh)
    ds = model.festim_problem.ds(surface_id)
    return E_CHARGE * assemble(_weight(model) * ufl.dot(charge_flux(model), n) * ds)


def generation_rate(model: IonChamberModel) -> float:
    """The production rate of ion pairs in the whole gas volume (1/s)."""
    return assemble(_weight(model) * model.pair_production * ufl.dx(domain=model.mesh))


def recombination_rate(model: IonChamberModel) -> float:
    """The loss rate of ion pairs to recombination in the whole gas volume (1/s)."""
    total = 0.0
    for i, reaction in enumerate(model.gas.reactions):
        if not reaction.is_recombination:
            continue
        # mass action: the rate coefficient times the density of every reactant
        local_rate = model.rates[i]
        for carrier in reaction.reactant:
            local_rate = local_rate * carrier.post_processing_solution
        total += assemble(_weight(model) * local_rate * ufl.dx(domain=model.mesh))
    return total


def space_charge_distortion(model: IonChamberModel) -> tuple[float, float]:
    """How strongly the charge in the gas would distort the applied field.

    Solves Poisson's equation for the converged charge density with grounded
    electrodes, and compares the result with the applied potential.

    Returns:
        the ratio of the L2 norms of the space-charge field and the applied field,
        and the largest space-charge potential relative to the bias
    """
    perturbation = model.space_charge.solve(E_CHARGE * net_charge_density(model))
    weight, dx = _weight(model), ufl.dx(domain=model.mesh)
    applied = ufl.grad(model.potential)
    induced = ufl.grad(perturbation)
    field_ratio = np.sqrt(
        assemble(weight * ufl.dot(induced, induced) * dx)
        / assemble(weight * ufl.dot(applied, applied) * dx)
    )
    largest = model.mesh.comm.allreduce(
        np.max(np.abs(perturbation.x.array)), op=MPI.MAX
    )
    return field_ratio, largest / abs(model.bias)


def compute_reading(model: IonChamberModel) -> Reading:
    """Evaluates every quantity of the :class:`Reading` on a solved model."""
    from ion_chambers.model import ANODE_ID, CATHODE_ID

    generation = generation_rate(model)
    recombination = recombination_rate(model)
    ideal_pairs = sum(
        source.ideal_pair_rate(model.geometry, model.gas) for source in model.sources
    )
    field_ratio, potential_ratio = space_charge_distortion(model)
    mean_net_charge = (
        assemble(
            _weight(model) * abs(net_charge_density(model)) * ufl.dx(domain=model.mesh)
        )
        / model.geometry.volume
    )
    nan = float("nan")
    return Reading(
        current=ramo_current(model),
        # negative charge leaves the gas through the anode when the bias is positive
        anode_current=-electrode_current(model, ANODE_ID),
        cathode_current=electrode_current(model, CATHODE_ID),
        balance_current=E_CHARGE * (generation - recombination),
        ideal_current=E_CHARGE * ideal_pairs,
        generation=generation,
        recombination=recombination,
        wall_efficiency=generation / ideal_pairs if ideal_pairs > 0 else nan,
        collection_efficiency=(
            (generation - recombination) / generation if generation > 0 else nan
        ),
        space_charge_field=float(field_ratio),
        space_charge_potential=float(potential_ratio),
        mean_net_charge_density=mean_net_charge,
    )
