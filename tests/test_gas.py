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
    """The package's Boltzmann constant in eV/K is exactly FESTIM's."""
    assert K_B_EV == F.k_B


@pytest.mark.parametrize("mobility", ["fast", None, True])
def test_carrier_rejects_wrong_mobility_type(mobility):
    """A mobility that is neither a number nor a callable raises a TypeError."""
    with pytest.raises(TypeError, match="mobility must be a float or a callable"):
        ic.ChargeCarrier(name="e", charge=-1, mobility=mobility)


@pytest.mark.parametrize("charge", [1.0, "+", True])
def test_carrier_rejects_wrong_charge_type(charge):
    """A charge that is not an int raises a TypeError."""
    with pytest.raises(TypeError, match="charge must be an int"):
        ic.ChargeCarrier(name="P", charge=charge, mobility=1e-4)


def test_carrier_rejects_zero_charge():
    """A neutral carrier raises a ValueError, since it would not drift."""
    with pytest.raises(ValueError, match="charge must be non-zero"):
        ic.ChargeCarrier(name="T2", charge=0, mobility=1e-4)


def test_carrier_accepts_callable_mobility():
    """A function of the reduced field is accepted as a mobility and stored as given."""
    carrier = ic.ChargeCarrier(name="e", charge=-1, mobility=lambda E_N: 0.1 / E_N)
    assert callable(carrier.mobility)


def test_reaction_rejects_reactants_not_in_a_list():
    """Reactants given as a bare string instead of a list raise a TypeError."""
    with pytest.raises(TypeError, match="reactants must be a list of str"):
        ic.Reaction(reactants="e", rate=1.0)


def test_reaction_needs_a_reactant():
    """A reaction with an empty list of reactants raises a ValueError."""
    with pytest.raises(ValueError, match="at least one carrier"):
        ic.Reaction(reactants=[], rate=1.0)


@pytest.mark.parametrize("products, expected", [([], True), (["N"], False)])
def test_reaction_is_recombination_when_it_has_no_product(products, expected):
    """A reaction counts as recombination only when it has no product."""
    reaction = ic.Reaction(reactants=["e"], rate=1.0, products=products)
    assert reaction.is_recombination is expected


@pytest.mark.parametrize("attribute", ["w_value", "density", "pressure", "temperature"])
def test_gas_rejects_non_positive_values(attribute):
    """A zero W value, density, pressure or temperature raises a ValueError."""
    with pytest.raises(ValueError, match=f"{attribute} must be positive"):
        make_gas(**{attribute: 0.0})


@pytest.mark.parametrize("attribute", ["w_value", "density", "pressure", "temperature"])
def test_gas_rejects_wrong_types(attribute):
    """A string for the W value, density, pressure or temperature raises a TypeError."""
    with pytest.raises(TypeError, match=f"{attribute} must be a float, not str"):
        make_gas(**{attribute: "high"})


def test_gas_rejects_reaction_with_unknown_carrier():
    """A reaction naming a carrier the gas does not have raises a ValueError."""
    with pytest.raises(ValueError, match="'N', which is not one of the carriers"):
        make_gas(reactions=[ic.Reaction(reactants=["P", "N"], rate=1e-12)])


def test_gas_rejects_unknown_ionisation_product():
    """An ionisation product the gas does not have raises a ValueError."""
    with pytest.raises(ValueError, match="ionisation_products refers to 'X'"):
        make_gas(ionisation_products=("e", "X"))


def test_gas_rejects_ionisation_products_of_equal_charge():
    """Ionisation products whose charges do not cancel raise a ValueError."""
    with pytest.raises(ValueError, match="opposite charges"):
        make_gas(ionisation_products=("e", "e"))


def test_gas_rejects_duplicate_carrier_names():
    """Two carriers with the same name raise a ValueError."""
    carriers = [ic.ChargeCarrier(name="e", charge=-1, mobility=0.1)] * 2
    with pytest.raises(ValueError, match="unique names"):
        make_gas(carriers=carriers)


def test_number_density_at_room_conditions():
    """An ideal gas at 1 atm and 293.15 K holds 2.5035e25 molecules per m3."""
    assert np.isclose(make_gas().number_density, 2.5035e25, rtol=1e-4)


def test_thermal_voltage_at_room_temperature():
    """The thermal voltage k_B T / e is 25.26 mV at 293.15 K."""
    assert np.isclose(make_gas().thermal_voltage, 0.02526, rtol=1e-3)


def test_carrier_lookup_by_name():
    """Looking up a carrier by its name returns that carrier."""
    assert make_gas().carrier("P").charge == 1


def test_carrier_lookup_fails_for_unknown_name():
    """Looking up a carrier the gas does not have raises a ValueError."""
    with pytest.raises(ValueError, match="'N' is not one of the carriers"):
        make_gas().carrier("N")
