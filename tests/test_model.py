"""Integration tests: full solves of small models."""

import numpy as np
import pytest

import ion_chambers as ic

pytestmark = pytest.mark.slow


def electron_ion_gas(mobility=0.04) -> ic.Gas:
    """A non-attaching gas: electrons and positive ions only."""
    return ic.Gas(
        name="non-attaching",
        w_value=33.97,
        density=1.205,
        carriers=[
            ic.ChargeCarrier(name="e", charge=-1, mobility=mobility, diffusivity=0.04),
            ic.ChargeCarrier(name="P", charge=1, mobility=2.0e-4),
        ],
        reactions=[ic.Reaction(reactants=["P", "e"], rate=2.0e-13)],
    )


# ---- input checks ----------------------------------------------------------------


def test_bias_must_be_non_zero(coarse_model):
    with pytest.raises(ValueError, match="bias must be non-zero"):
        coarse_model.bias = 0.0


def test_bias_rejects_wrong_type(coarse_model):
    with pytest.raises(TypeError, match="bias must be a float, not str"):
        coarse_model.bias = "500 V"


def test_run_before_initialise_fails(coarse_model):
    with pytest.raises(RuntimeError, match=r"call initialise\(\) before run\(\)"):
        coarse_model.run()


def test_a_single_source_is_wrapped_in_a_list(coarse_model):
    assert len(coarse_model.sources) == 1


# ---- regression against the prototype (chamber_model.py, N_r = 60) ------------------


@pytest.mark.parametrize(
    "activity, expected",
    [
        (3.7e4, 9.40662259521735e-16),
        (3.7e10, 9.33952139344237e-10),
        (3.7e13, 3.2356679629878293e-07),
    ],
)
def test_current_matches_the_prototype(reference_readings, activity, expected):
    assert np.isclose(reference_readings[activity].current, expected, rtol=1e-8)


def test_wall_efficiency_matches_the_prototype(reference_readings):
    reading = reference_readings[3.7e4]
    assert np.isclose(reading.wall_efficiency, 0.9895863173320661, rtol=1e-8)


def test_collection_efficiency_matches_the_prototype(reference_readings):
    reading = reference_readings[3.7e13]
    assert np.isclose(reading.collection_efficiency, 0.36235411253601335, rtol=1e-8)


def test_space_charge_field_matches_the_prototype(reference_readings):
    reading = reference_readings[3.7e13]
    assert np.isclose(reading.space_charge_field, 0.6501553934189711, rtol=1e-8)


# ---- physics checks ----------------------------------------------------------------


@pytest.mark.parametrize("activity", [3.7e4, 3.7e10, 3.7e13])
def test_ramo_current_agrees_with_the_anode_flux(reference_readings, activity):
    reading = reference_readings[activity]
    assert np.isclose(reading.current, reading.anode_current, rtol=1e-3)


@pytest.mark.parametrize("activity", [3.7e4, 3.7e10])
def test_both_electrodes_carry_the_same_current(reference_readings, activity):
    reading = reference_readings[activity]
    assert np.isclose(reading.anode_current, reading.cathode_current, rtol=5e-3)


def test_response_is_linear_at_low_activity(reference_readings):
    low, ten_times = reference_readings[3.7e4], reference_readings[3.7e5]
    assert np.isclose(ten_times.current, 10 * low.current, rtol=1e-5)


def test_recombination_is_negligible_at_low_activity(reference_readings):
    assert reference_readings[3.7e4].collection_efficiency > 1 - 1e-6


def test_recombination_dominates_at_high_activity(reference_readings):
    assert reference_readings[3.7e13].collection_efficiency < 0.5


def test_space_charge_grows_with_activity(reference_readings):
    fields = [r.space_charge_field for r in reference_readings.values()]
    assert fields == sorted(fields)


def test_current_never_exceeds_the_ideal_current(reference_readings):
    assert all(r.current < r.ideal_current for r in reference_readings.values())


def test_no_wall_loss_generates_every_pair(coarse_model):
    coarse_model.initialise()
    assert np.isclose(coarse_model.run().wall_efficiency, 1.0)


def test_reversing_the_bias_reverses_the_current(coarse_model):
    coarse_model.initialise()
    forward = coarse_model.run().current
    coarse_model.bias = -coarse_model.bias
    coarse_model.initialise()
    assert np.isclose(coarse_model.run().current, -forward, rtol=1e-6)


def test_higher_bias_collects_more(coarse_model):
    coarse_model.sources[0].activity = 3.7e12
    efficiencies = []
    for bias in (100.0, 1000.0):
        coarse_model.bias = bias
        coarse_model.initialise()
        efficiencies.append(coarse_model.run().collection_efficiency)
    assert efficiencies[1] > efficiencies[0]


def test_lower_pressure_loses_more_to_the_walls(coarse_model):
    efficiencies = []
    for pressure in (101325.0, 20265.0):
        coarse_model.gas = ic.library.air(pressure=pressure)
        coarse_model.sources = [ic.TritiumInGas(activity=3.7e4, n_rays=300)]
        coarse_model.initialise()
        efficiencies.append(coarse_model.run().wall_efficiency)
    assert efficiencies[1] < efficiencies[0]


# ---- flexibility of the gas description ---------------------------------------------


def test_gas_without_negative_ions_runs(coarse_model):
    coarse_model.gas = electron_ion_gas()
    coarse_model.initialise()
    reading = coarse_model.run()
    assert np.isclose(reading.current, reading.ideal_current, rtol=0.1)


def test_gas_without_reactions_has_no_recombination(coarse_model):
    gas = electron_ion_gas()
    gas.reactions = []
    coarse_model.gas = gas
    coarse_model.initialise()
    assert coarse_model.run().recombination == 0.0


def test_callable_mobility_matches_the_equivalent_constant(coarse_model):
    currents = []
    for mobility in (0.04, lambda reduced_field: np.full_like(reduced_field, 0.04)):
        coarse_model.gas = electron_ion_gas(mobility=mobility)
        coarse_model.initialise()
        currents.append(coarse_model.run().current)
    assert np.isclose(currents[0], currents[1], rtol=1e-9)


def test_callable_coefficients_see_the_reduced_field(coarse_model):
    """Near the rod of the reference chamber at 500 V, E/N is a few townsend."""
    coarse_model.initialise()
    assert 1.0 < coarse_model.reduced_field.max() < 10.0


def test_two_sources_add_up(coarse_model):
    coarse_model.sources = [
        ic.TritiumInGas(activity=1e4, wall_loss=False),
        ic.TritiumInGas(activity=3e4, wall_loss=False),
    ]
    coarse_model.initialise()
    reading = coarse_model.run()
    one_source = ic.TritiumInGas(activity=4e4).ideal_pair_rate(
        coarse_model.geometry, coarse_model.gas
    )
    assert np.isclose(reading.generation, one_source)


def test_response_curve_extrapolates_from_two_meshes():
    def build_model(resolution: int) -> ic.IonChamberModel:
        return ic.IonChamberModel(
            geometry=ic.library.reference_chamber(),
            gas=ic.library.air(),
            sources=ic.TritiumInGas(activity=0.0, wall_loss=False),
            resolution=resolution,
        )

    (row,) = ic.response_curve(build_model, [3.7e4], resolutions=(20, 40))
    expected = 2 * row["current_fine"] - row["current_coarse"]
    assert np.isclose(row["current_extrapolated"], expected)
