import festim as F
import numpy as np
import pytest

import ion_chambers as ic
from ion_chambers.constants import K_B_EV


def make_gas(**kwargs) -> ic.Gas:
    arguments = {
        "name": "test gas",
        "w_value": 30.0,
        "density": 1.0,
        "carriers": [
            ic.ChargeCarrier(name="e", charge=-1, mobility=0.1),
            ic.ChargeCarrier(name="P", charge=1, mobility=1e-4),
        ],
    }
    arguments.update(kwargs)
    return ic.Gas(**arguments)


def test_boltzmann_constant_matches_festim():
    assert K_B_EV == F.k_B


@pytest.mark.parametrize("mobility", ["fast", None, True])
def test_carrier_rejects_wrong_mobility_type(mobility):
    with pytest.raises(TypeError, match="mobility must be a float or a callable"):
        ic.ChargeCarrier(name="e", charge=-1, mobility=mobility)


@pytest.mark.parametrize("charge", [1.0, "+", True])
def test_carrier_rejects_wrong_charge_type(charge):
    with pytest.raises(TypeError, match="charge must be an int"):
        ic.ChargeCarrier(name="P", charge=charge, mobility=1e-4)


def test_carrier_rejects_zero_charge():
    with pytest.raises(ValueError, match="charge must be non-zero"):
        ic.ChargeCarrier(name="T2", charge=0, mobility=1e-4)


def test_carrier_accepts_callable_mobility():
    carrier = ic.ChargeCarrier(name="e", charge=-1, mobility=lambda E_N: 0.1 / E_N)
    assert callable(carrier.mobility)


def test_reaction_rejects_reactants_not_in_a_list():
    with pytest.raises(TypeError, match="reactants must be a list of str"):
        ic.Reaction(reactants="e", rate=1.0)


def test_reaction_needs_a_reactant():
    with pytest.raises(ValueError, match="at least one carrier"):
        ic.Reaction(reactants=[], rate=1.0)


@pytest.mark.parametrize("products, expected", [([], True), (["N"], False)])
def test_reaction_is_recombination_when_it_has_no_product(products, expected):
    reaction = ic.Reaction(reactants=["e"], rate=1.0, products=products)
    assert reaction.is_recombination is expected


@pytest.mark.parametrize("attribute", ["w_value", "density", "pressure", "temperature"])
def test_gas_rejects_non_positive_values(attribute):
    with pytest.raises(ValueError, match=f"{attribute} must be positive"):
        make_gas(**{attribute: 0.0})


@pytest.mark.parametrize("attribute", ["w_value", "density", "pressure", "temperature"])
def test_gas_rejects_wrong_types(attribute):
    with pytest.raises(TypeError, match=f"{attribute} must be a float, not str"):
        make_gas(**{attribute: "high"})


def test_gas_rejects_reaction_with_unknown_carrier():
    with pytest.raises(ValueError, match="'N', which is not one of the carriers"):
        make_gas(reactions=[ic.Reaction(reactants=["P", "N"], rate=1e-12)])


def test_gas_rejects_unknown_ionisation_product():
    with pytest.raises(ValueError, match="ionisation_products refers to 'X'"):
        make_gas(ionisation_products=("e", "X"))


def test_gas_rejects_ionisation_products_of_equal_charge():
    with pytest.raises(ValueError, match="opposite charges"):
        make_gas(ionisation_products=("e", "e"))


def test_gas_rejects_duplicate_carrier_names():
    carriers = [ic.ChargeCarrier(name="e", charge=-1, mobility=0.1)] * 2
    with pytest.raises(ValueError, match="unique names"):
        make_gas(carriers=carriers)


def test_number_density_at_room_conditions():
    assert np.isclose(make_gas().number_density, 2.5035e25, rtol=1e-4)


def test_thermal_voltage_at_room_temperature():
    assert np.isclose(make_gas().thermal_voltage, 0.02526, rtol=1e-3)


def test_carrier_lookup_by_name():
    assert make_gas().carrier("P").charge == 1


def test_carrier_lookup_fails_for_unknown_name():
    with pytest.raises(ValueError, match="'N' is not one of the carriers"):
        make_gas().carrier("N")


def test_air_has_three_carriers():
    assert [c.name for c in ic.library.air().carriers] == ["e", "P", "N"]


def test_air_density_follows_pressure():
    assert np.isclose(ic.library.air(pressure=50662.5).density, 1.205 / 2)


def test_air_ion_mobility_is_inversely_proportional_to_pressure():
    low = ic.library.air(pressure=50662.5).carrier("P").mobility
    assert np.isclose(low, 2 * ic.library.air().carrier("P").mobility)


def test_air_attachment_is_three_body():
    """Halving the pressure divides the attachment frequency by four."""
    low = ic.library.air(pressure=50662.5).reactions[0].rate
    assert np.isclose(low, ic.library.air().reactions[0].rate / 4)


def test_air_density_follows_temperature():
    assert np.isclose(ic.library.air(temperature=586.3).density, 1.205 / 2)
