"""Shared fixtures."""

import pytest
from prototype_case import air, reference_chamber, solve_potential

import ion_chambers as ic

ACTIVITIES = [3.7e4, 3.7e5, 3.7e10, 3.7e13]  # Bq/m3


@pytest.fixture
def chamber() -> ic.CoaxialChamber:
    """The reference 1 L chamber."""
    return reference_chamber()


@pytest.fixture
def coarse_electrostatics(chamber) -> ic.ElectrostaticProblem:
    """The potential of the reference chamber at 500 V on a coarse mesh, solved."""
    return solve_potential(chamber, resolution=30)


@pytest.fixture
def coarse_transport(coarse_electrostatics) -> ic.ChargeTransportProblem:
    """A cheap transport problem, not initialised: coarse mesh, no wall losses."""
    return ic.ChargeTransportProblem(
        electrostatic_problem=coarse_electrostatics,
        gas=air(),
        ionisation_sources=[ic.TritiumInGas(activity=3.7e4)],
    )


@pytest.fixture(scope="session")
def reference_results() -> dict[str, ic.ChamberQuantity]:
    """The quantities of the reference chamber in air at 500 V, 60 radial cells.

    The problem is run once for each of ACTIVITIES, so ``data`` of every quantity
    holds one value per activity, in that order.
    """
    chamber = reference_chamber()
    quantities = {
        "current": ic.ChamberCurrent(surface=chamber.anode),
        "anode_current": ic.ElectrodeCurrent(surface=chamber.anode),
        "cathode_current": ic.ElectrodeCurrent(surface=chamber.cathode),
        "ideal_current": ic.IdealCurrent(),
        "generation": ic.PairGeneration(),
        "recombination": ic.Recombination(),
        "wall_efficiency": ic.WallEfficiency(),
        "collection_efficiency": ic.CollectionEfficiency(),
        "space_charge_field": ic.SpaceChargeField(),
    }
    tritium = ic.TritiumInGas(activity=0.0, geometry=chamber)
    transport = ic.ChargeTransportProblem(
        electrostatic_problem=solve_potential(chamber, resolution=60),
        gas=air(),
        ionisation_sources=[tritium],
        exports=list(quantities.values()),
    )
    transport.initialise()
    for activity in ACTIVITIES:
        tritium.activity = activity
        transport.run()
    return quantities
