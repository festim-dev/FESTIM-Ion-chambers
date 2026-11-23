"""The carrier gas: its charge carriers, their reactions and its ionisation data.

A :class:`Gas` only holds user input. The package ships no gas data: the coefficients
are yours to supply. The documentation shows a complete example for air.
"""

from dataclasses import dataclass, field

from ion_chambers.carrier_reaction import CarrierReaction
from ion_chambers.charge_carrier import ChargeCarrier
from ion_chambers.constants import K_B, K_B_EV


@dataclass
class Gas:
    """The gas filling the chamber, at a given pressure and temperature.

    All coefficients are the values *at this pressure and temperature*: the class does
    not rescale them. The air example in the documentation shows how to write a
    function that builds a gas whose coefficients follow the pressure.

    Args:
        name: the name of the gas
        w_value: the mean energy spent per ion pair (eV)
        density: the mass density (kg/m3), which sets the range of the betas
        carriers: the charge carriers
        ionisation_products: the two carriers created by each ionisation, one
            negative and one positive
        reactions: the reactions between carriers
        pressure: the pressure (Pa)
        temperature: the temperature (K)
    """

    name: str
    w_value: float
    density: float
    carriers: list[ChargeCarrier]
    ionisation_products: tuple[ChargeCarrier, ChargeCarrier]
    reactions: list[CarrierReaction] = field(default_factory=list)
    pressure: float = 101325.0
    temperature: float = 293.15

    def __post_init__(self) -> None:
        for attribute in ("w_value", "density", "pressure", "temperature"):
            value = getattr(self, attribute)
            if isinstance(value, bool) or not isinstance(value, int | float):
                raise TypeError(
                    f"{attribute} must be a float, not {type(value).__name__}"
                )
            if value <= 0:
                raise ValueError(f"{attribute} must be positive, not {value}")

        for carrier in self.carriers:
            if not isinstance(carrier, ChargeCarrier):
                raise TypeError(
                    f"carriers must be ion_chambers.ChargeCarrier objects, "
                    f"not {type(carrier).__name__}"
                )

        names = [carrier.name for carrier in self.carriers]
        if len(set(names)) != len(names):
            raise ValueError(f"carriers must have unique names, not {names}")

        for reaction in self.reactions:
            if not isinstance(reaction, CarrierReaction):
                raise TypeError(
                    f"reactions must be ion_chambers.CarrierReaction objects, "
                    f"not {type(reaction).__name__}"
                )
            for species in reaction.reactant + reaction.product:
                if not self._has_carrier(species):
                    raise ValueError(
                        f"a reaction involves '{species.name}', which is not one of "
                        f"the carriers of the gas {names}"
                    )

        total_charge = 0
        for carrier in self.ionisation_products:
            if not self._has_carrier(carrier):
                carrier_name = getattr(carrier, "name", carrier)
                raise ValueError(
                    f"ionisation_products contains '{carrier_name}', which is not "
                    f"one of the carriers of the gas {names}"
                )
            total_charge += carrier.charge
        if total_charge != 0:
            raise ValueError("ionisation_products must have opposite charges")

    def _has_carrier(self, carrier: object) -> bool:
        """Whether this exact carrier object belongs to the gas."""
        for own_carrier in self.carriers:
            if own_carrier is carrier:
                return True
        return False

    @property
    def number_density(self) -> float:
        """The number density of the gas molecules, ideal gas (1/m3)."""
        return self.pressure / (K_B * self.temperature)

    @property
    def thermal_voltage(self) -> float:
        """``k_B T / e`` (V), the ratio D / mobility of a thermal carrier."""
        return K_B_EV * self.temperature
