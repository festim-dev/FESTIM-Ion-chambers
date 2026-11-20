"""The carrier gas: its charge carriers, their reactions and its ionisation data.

Everything here is plain user input. A :class:`Gas` holds numbers (or callables of the
reduced electric field); no DOLFINx object is created until
:meth:`ion_chambers.IonChamberModel.initialise`.

Ready-made examples live in :mod:`ion_chambers.library`.
"""

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from ion_chambers.constants import K_B, K_B_EV

Coefficient = float | Callable[[np.ndarray], np.ndarray]
"""A transport or rate coefficient.

Either a constant, or a callable taking the reduced electric field E/N in townsend
(1 Td = 1e-21 V m2) as a numpy array and returning the coefficient in SI units.
"""


def _check_coefficient(name: str, value: object) -> None:
    if isinstance(value, bool) or not (
        isinstance(value, int | float) or callable(value)
    ):
        raise TypeError(
            f"{name} must be a float or a callable of the reduced field, "
            f"not {type(value).__name__}"
        )


@dataclass
class ChargeCarrier:
    """A charged species transported by drift and diffusion.

    Args:
        name: the name of the species, used to refer to it in reactions
        charge: the charge number (-1 for electrons and negative ions, +1 for positive
            ions)
        mobility: the mobility at the gas conditions (m2/V/s), positive
        diffusivity: the diffusion coefficient (m2/s). ``None`` uses the Einstein
            relation ``D = mobility * k_B T / e``, which is right for thermal ions but
            not for electrons heated by the field.
    """

    name: str
    charge: int
    mobility: Coefficient
    diffusivity: Coefficient | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.name, str):
            raise TypeError(f"name must be a str, not {type(self.name).__name__}")
        if isinstance(self.charge, bool) or not isinstance(self.charge, int):
            raise TypeError(f"charge must be an int, not {type(self.charge).__name__}")
        if self.charge == 0:
            raise ValueError("charge must be non-zero: neutral species do not drift")
        _check_coefficient("mobility", self.mobility)
        if self.diffusivity is not None:
            _check_coefficient("diffusivity", self.diffusivity)


@dataclass
class Reaction:
    """A volume reaction between charge carriers, with a mass-action rate.

    The rate is ``rate * product of the reactant densities``, so ``rate`` is in 1/s for
    one reactant (attachment) and m3/s for two (recombination).

    Args:
        reactants: the names of the reacting carriers
        rate: the rate coefficient at the gas conditions
        products: the names of the carriers produced. Empty when the charges
            neutralise each other (recombination).
    """

    reactants: list[str]
    rate: Coefficient
    products: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        for attribute in ("reactants", "products"):
            value = getattr(self, attribute)
            if not isinstance(value, list):
                raise TypeError(
                    f"{attribute} must be a list of str, not {type(value).__name__}"
                )
        if not self.reactants:
            raise ValueError("reactants must name at least one carrier")
        _check_coefficient("rate", self.rate)

    @property
    def is_recombination(self) -> bool:
        """Whether the reaction removes charge pairs (it has no product)."""
        return not self.products


@dataclass
class Gas:
    """The gas filling the chamber, at a given pressure and temperature.

    All coefficients are the values *at this pressure and temperature*: the class does
    not rescale them. The functions in :mod:`ion_chambers.library` show how to build a
    gas whose coefficients follow the pressure.

    Args:
        name: the name of the gas
        w_value: the mean energy spent per ion pair (eV)
        density: the mass density (kg/m3), which sets the range of the betas
        carriers: the charge carriers
        reactions: the reactions between carriers
        pressure: the pressure (Pa)
        temperature: the temperature (K)
        ionisation_products: the names of the two carriers created by each ionisation
    """

    name: str
    w_value: float
    density: float
    carriers: list[ChargeCarrier]
    reactions: list[Reaction] = field(default_factory=list)
    pressure: float = 101325.0
    temperature: float = 293.15
    ionisation_products: tuple[str, str] = ("e", "P")

    def __post_init__(self) -> None:
        for attribute in ("w_value", "density", "pressure", "temperature"):
            value = getattr(self, attribute)
            if isinstance(value, bool) or not isinstance(value, int | float):
                raise TypeError(
                    f"{attribute} must be a float, not {type(value).__name__}"
                )
            if value <= 0:
                raise ValueError(f"{attribute} must be positive, not {value}")
        names = [carrier.name for carrier in self.carriers]
        if len(set(names)) != len(names):
            raise ValueError(f"carriers must have unique names, not {names}")
        for reaction in self.reactions:
            for name in reaction.reactants + reaction.products:
                if name not in names:
                    raise ValueError(
                        f"reaction refers to '{name}', which is not one of the "
                        f"carriers {names}"
                    )
        for name in self.ionisation_products:
            if name not in names:
                raise ValueError(
                    f"ionisation_products refers to '{name}', which is not one of the "
                    f"carriers {names}"
                )
        if sum(self.carrier(name).charge for name in self.ionisation_products) != 0:
            raise ValueError("ionisation_products must have opposite charges")

    def carrier(self, name: str) -> ChargeCarrier:
        """Returns the carrier called ``name``."""
        for carrier in self.carriers:
            if carrier.name == name:
                return carrier
        raise ValueError(f"'{name}' is not one of the carriers of {self.name}")

    @property
    def number_density(self) -> float:
        """The number density of the gas molecules, ideal gas (1/m3)."""
        return self.pressure / (K_B * self.temperature)

    @property
    def thermal_voltage(self) -> float:
        """``k_B T / e`` (V), the ratio D / mobility of a thermal carrier."""
        return K_B_EV * self.temperature
