import festim as F
import numpy as np
import pytest

import ion_chambers as ic
from ion_chambers.constants import K_B_EV


def make_carriers() -> tuple[ic.ChargeCarrier, ic.ChargeCarrier]:
    electrons = ic.ChargeCarrier(name="e", charge=-1, mobility=0.1)
    positive_ions = ic.ChargeCarrier(name="P", charge=1, mobility=1e-4)
    return electrons, positive_ions


def make_gas(**kwargs) -> ic.Gas:
    electrons, positive_ions = make_carriers()
    arguments = {
        "name": "test gas",
        "w_value": 30.0,
        "density": 1.0,
        "carriers": [electrons, positive_ions],
        "ionisation_products": (electrons, positive_ions),
    }
    arguments.update(kwargs)
    return ic.Gas(**arguments)


def test_boltzmann_constant_matches_festim():
    """The package's Boltzmann constant in eV/K is exactly FESTIM's."""
    assert K_B_EV == F.k_B


# ---- ChargeCarrier ------------------------------------------------------------------


def test_carrier_is_a_festim_species():
    """A ChargeCarrier can be used wherever FESTIM expects a Species."""
    electrons, _ = make_carriers()
    assert isinstance(electrons, F.Species)


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


def test_carrier_rejects_wrong_name_type():
    """A name that is not a string raises a TypeError."""
    with pytest.raises(TypeError, match="name must be a str, not int"):
        ic.ChargeCarrier(name=1, charge=1, mobility=1e-4)


def test_carrier_accepts_callable_mobility():
    """A function of the reduced field is accepted as a mobility and stored as
    given."""
    carrier = ic.ChargeCarrier(name="e", charge=-1, mobility=lambda E_N: 0.1 / E_N)
    assert callable(carrier.mobility)


# ---- CarrierReaction ----------------------------------------------------------------


def test_reaction_is_a_festim_generic_reaction():
    """A CarrierReaction can be used wherever FESTIM expects a GenericReaction."""
    electrons, positive_ions = make_carriers()
    reaction = ic.CarrierReaction(reactant=[positive_ions, electrons], rate=1e-13)
    assert isinstance(reaction, F.GenericReaction)


def test_reaction_needs_no_volume_when_created():
    """A CarrierReaction has no volume until a model assigns one."""
    electrons, positive_ions = make_carriers()
    reaction = ic.CarrierReaction(reactant=[positive_ions, electrons], rate=1e-13)
    assert reaction.volume is None


def test_reaction_rejects_a_volume_of_the_wrong_type():
    """Assigning something that is not a VolumeSubdomain as volume raises a
    TypeError."""
    electrons, positive_ions = make_carriers()
    reaction = ic.CarrierReaction(reactant=[positive_ions, electrons], rate=1e-13)
    with pytest.raises(TypeError, match="volume must be a festim"):
        reaction.volume = "gas"


def test_reaction_accepts_a_single_reactant():
    """A reactant given on its own is stored as a list of one, as in FESTIM."""
    electrons, _ = make_carriers()
    reaction = ic.CarrierReaction(reactant=electrons, rate=1e8)
    assert reaction.reactant == [electrons]


@pytest.mark.parametrize("rate", ["fast", None])
def test_reaction_rejects_wrong_rate_type(rate):
    """A rate that is neither a number nor a callable raises a TypeError."""
    electrons, _ = make_carriers()
    with pytest.raises(TypeError, match="rate must be a float or a callable"):
        ic.CarrierReaction(reactant=electrons, rate=rate)


def test_reaction_without_product_is_recombination():
    """A reaction with no product counts as recombination."""
    electrons, positive_ions = make_carriers()
    reaction = ic.CarrierReaction(reactant=[positive_ions, electrons], rate=1e-13)
    assert reaction.is_recombination is True


def test_reaction_with_product_is_not_recombination():
    """A reaction with a product, such as attachment, is not recombination."""
    electrons, _ = make_carriers()
    negative_ions = ic.ChargeCarrier(name="N", charge=-1, mobility=1e-4)
    reaction = ic.CarrierReaction(reactant=electrons, product=negative_ions, rate=1e8)
    assert reaction.is_recombination is False


# ---- Gas ----------------------------------------------------------------------------


@pytest.mark.parametrize("attribute", ["w_value", "density", "pressure", "temperature"])
def test_gas_rejects_non_positive_values(attribute):
    """A zero W value, density, pressure or temperature raises a ValueError."""
    with pytest.raises(ValueError, match=f"{attribute} must be positive"):
        make_gas(**{attribute: 0.0})


@pytest.mark.parametrize("attribute", ["w_value", "density", "pressure", "temperature"])
def test_gas_rejects_wrong_types(attribute):
    """A string for the W value, density, pressure or temperature raises a
    TypeError."""
    with pytest.raises(TypeError, match=f"{attribute} must be a float, not str"):
        make_gas(**{attribute: "high"})


def test_gas_rejects_a_plain_festim_species_as_carrier():
    """A carrier that is a plain festim.Species, with no charge, raises a
    TypeError."""
    with pytest.raises(TypeError, match="carriers must be ion_chambers"):
        make_gas(carriers=[F.Species("H")])


def test_gas_rejects_a_plain_festim_reaction():
    """A reaction that is not a CarrierReaction raises a TypeError."""
    with pytest.raises(TypeError, match="reactions must be ion_chambers"):
        make_gas(reactions=["P + e -> 0"])


def test_gas_rejects_reaction_with_a_carrier_it_does_not_have():
    """A reaction involving a carrier that is not in the gas raises a ValueError."""
    electrons, positive_ions = make_carriers()
    stranger = ic.ChargeCarrier(name="N", charge=-1, mobility=1e-4)
    reaction = ic.CarrierReaction(reactant=[positive_ions, stranger], rate=1e-12)
    with pytest.raises(ValueError, match="'N', which is not one of the carriers"):
        ic.Gas(
            name="test gas",
            w_value=30.0,
            density=1.0,
            carriers=[electrons, positive_ions],
            ionisation_products=(electrons, positive_ions),
            reactions=[reaction],
        )


def test_gas_rejects_ionisation_product_it_does_not_have():
    """An ionisation product that is not a carrier of the gas raises a ValueError."""
    electrons, positive_ions = make_carriers()
    stranger = ic.ChargeCarrier(name="X", charge=1, mobility=1e-4)
    with pytest.raises(ValueError, match="ionisation_products contains 'X'"):
        ic.Gas(
            name="test gas",
            w_value=30.0,
            density=1.0,
            carriers=[electrons, positive_ions],
            ionisation_products=(electrons, stranger),
        )


def test_gas_rejects_ionisation_products_of_equal_charge():
    """Ionisation products whose charges do not cancel raise a ValueError."""
    electrons, positive_ions = make_carriers()
    with pytest.raises(ValueError, match="opposite charges"):
        ic.Gas(
            name="test gas",
            w_value=30.0,
            density=1.0,
            carriers=[electrons, positive_ions],
            ionisation_products=(electrons, electrons),
        )


def test_gas_rejects_duplicate_carrier_names():
    """Two carriers with the same name raise a ValueError."""
    first = ic.ChargeCarrier(name="e", charge=-1, mobility=0.1)
    second = ic.ChargeCarrier(name="e", charge=1, mobility=0.1)
    with pytest.raises(ValueError, match="unique names"):
        ic.Gas(
            name="test gas",
            w_value=30.0,
            density=1.0,
            carriers=[first, second],
            ionisation_products=(first, second),
        )


def test_number_density_at_room_conditions():
    """An ideal gas at 1 atm and 293.15 K holds 2.5035e25 molecules per m3."""
    assert np.isclose(make_gas().number_density, 2.5035e25, rtol=1e-4)


def test_thermal_voltage_at_room_temperature():
    """The thermal voltage k_B T / e is 25.26 mV at 293.15 K."""
    assert np.isclose(make_gas().thermal_voltage, 0.02526, rtol=1e-3)
