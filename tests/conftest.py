"""Shared fixtures."""

import pytest
from prototype_case import air, reference_chamber

import ion_chambers as ic

ACTIVITIES = [3.7e4, 3.7e5, 3.7e10, 3.7e13]  # Bq/m3


@pytest.fixture
def chamber() -> ic.CoaxialChamber:
    """The reference 1 L chamber."""
    return reference_chamber()


@pytest.fixture
def coarse_model() -> ic.IonChamberModel:
    """A cheap model, not initialised: coarse mesh, no wall-loss Monte Carlo."""
    return ic.IonChamberModel(
        geometry=reference_chamber(),
        gas=air(),
        sources=ic.TritiumInGas(activity=3.7e4, wall_loss=False),
        resolution=30,
    )


@pytest.fixture(scope="session")
def reference_model() -> ic.IonChamberModel:
    """The reference chamber in air at 500 V, initialised, 60 radial cells."""
    model = ic.IonChamberModel(
        geometry=reference_chamber(),
        gas=air(),
        sources=ic.TritiumInGas(activity=0.0),
        bias=500.0,
        resolution=60,
    )
    model.initialise()
    return model


@pytest.fixture(scope="session")
def reference_readings(reference_model) -> dict[float, ic.Reading]:
    """Readings of the reference model, keyed by activity (Bq/m3)."""
    readings = {}
    for activity in ACTIVITIES:
        reference_model.sources[0].activity = activity
        readings[activity] = reference_model.run()
    return readings
