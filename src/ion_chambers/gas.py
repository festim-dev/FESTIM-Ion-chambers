"""The carrier gas: its charge carriers, their reactions and its ionisation data.

Everything here is user input. The package ships no gas data: the coefficients are
yours to supply. The documentation shows a complete example for air.
"""

from __future__ import annotations

from collections.abc import Callable

import festim as F
import numpy as np

from ion_chambers.constants import K_B, K_B_EV

Coefficient = float | Callable[[np.ndarray], np.ndarray]
"""A transport or rate coefficient.

Either a constant, or a callable taking the reduced electric field E/N in townsend
(1 Td = 1e-21 V m2) as a numpy array and returning the coefficient in SI units.
"""


def check_coefficient(name: str, value: object) -> None:
    """Raises a TypeError unless ``value`` is a number or a callable.

    Args:
        name: the name of the attribute, used in the error message
        value: the value to check
    """
    is_number = isinstance(value, int | float) and not isinstance(value, bool)
    if not is_number and not callable(value):
        raise TypeError(
            f"{name} must be a float or a callable of the reduced field, "
            f"not {type(value).__name__}"
        )


def check_positive_number(name: str, value: object) -> None:
    """Raises an error unless ``value`` is a positive number.

    Args:
        name: the name of the attribute, used in the error message
        value: the value to check
    """
    is_number = isinstance(value, int | float) and not isinstance(value, bool)
    if not is_number:
        raise TypeError(f"{name} must be a float, not {type(value).__name__}")
    if value <= 0:
        raise ValueError(f"{name} must be positive, not {value}")


class ChargeCarrier(F.Species):
    """A charged species transported by drift and diffusion.

    A :class:`festim.Species` with a charge and a mobility. It can be used wherever
    FESTIM expects a species, for instance in ``festim.SpeciesExport``.

    A carrier holds the solution of the problem it belongs to, like any FESTIM species,
    so do not share one carrier between two problems that are alive at the same time.

    Args:
        name: the name of the species
        charge: the charge number (-1 for electrons and negative ions, +1 for positive
            ions)
        mobility: the mobility at the gas conditions (m2/V/s), positive. A float, or a
            callable of the reduced field E/N in townsend.
        diffusivity: the diffusion coefficient (m2/s). ``None`` uses the Einstein
            relation ``D = mobility * k_B T / e``, which is right for thermal ions but
            not for electrons heated by the field.

    Attributes:
        name: the name of the species
        charge: the charge number
        mobility: the mobility at the gas conditions (m2/V/s)
        diffusivity: the diffusion coefficient (m2/s), or ``None``
    """

    def __init__(
        self,
        name: str,
        charge: int,
        mobility: Coefficient,
        diffusivity: Coefficient | None = None,
    ) -> None:
        if not isinstance(name, str):
            raise TypeError(f"name must be a str, not {type(name).__name__}")
        super().__init__(name=name)

        if isinstance(charge, bool) or not isinstance(charge, int):
            raise TypeError(f"charge must be an int, not {type(charge).__name__}")
        if charge == 0:
            raise ValueError("charge must be non-zero: neutral species do not drift")
        self.charge = charge

        check_coefficient("mobility", mobility)
        self.mobility = mobility

        if diffusivity is not None:
            check_coefficient("diffusivity", diffusivity)
        self.diffusivity = diffusivity


class CarrierReaction(F.GenericReaction):
    """A volume reaction between charge carriers, with a mass-action rate.

    A :class:`festim.GenericReaction` that differs from its parent in two ways:

    * It needs no volume subdomain when it is created. The transport problem owns the
      gas volume and assigns it during ``initialise``.
    * Its rate coefficient may be a callable of the reduced electric field E/N, which
      FESTIM's rate coefficients do not know about. The transport problem evaluates
      ``rate`` in the local field during ``initialise`` and stores the result as
      ``forward_rate``.

    The rate of the reaction is ``rate`` times the densities of the reactants, so
    ``rate`` is in 1/s for one reactant (attachment) and m3/s for two (recombination).

    Args:
        reactant: the reacting carrier, or a list of them
        rate: the rate coefficient at the gas conditions. A float, or a callable of
            the reduced field E/N in townsend.
        product: the carrier produced, or a list of them. ``None`` when the charges
            neutralise each other (recombination).

    Attributes:
        reactant: the reacting carriers, as a list
        product: the carriers produced, as a list (empty for recombination)
        rate: the rate coefficient as given by the user
        forward_rate: the rate coefficient FESTIM uses, set by the transport problem
        volume: the gas volume, set by the transport problem
    """

    def __init__(
        self,
        reactant: ChargeCarrier | list[ChargeCarrier],
        rate: Coefficient,
        product: ChargeCarrier | list[ChargeCarrier] | None = None,
    ) -> None:
        check_coefficient("rate", rate)
        self.rate = rate
        # volume and forward_rate are placeholders until the transport problem sets
        # them
        super().__init__(
            volume=None, reactant=reactant, product=product, forward_rate=0.0
        )

    @property
    def volume(self) -> F.VolumeSubdomain | None:
        return self._volume

    @volume.setter
    def volume(self, value: F.VolumeSubdomain | None) -> None:
        if value is not None and not isinstance(value, F.VolumeSubdomain):
            raise TypeError(
                f"volume must be a festim.VolumeSubdomain or None, "
                f"not {type(value).__name__}"
            )
        self._volume = value

    @property
    def is_recombination(self) -> bool:
        """Whether the reaction removes charge pairs (it has no product)."""
        return len(self.product) == 0


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

    Attributes:
        name: the name of the gas
        w_value: the mean energy spent per ion pair (eV)
        density: the mass density (kg/m3)
        carriers: the charge carriers
        ionisation_products: the two carriers created by each ionisation
        reactions: the reactions between carriers
        pressure: the pressure (Pa)
        temperature: the temperature (K)
    """

    def __init__(
        self,
        name: str,
        w_value: float,
        density: float,
        carriers: list[ChargeCarrier],
        ionisation_products: tuple[ChargeCarrier, ChargeCarrier],
        reactions: list[CarrierReaction] | None = None,
        pressure: float = 101325.0,
        temperature: float = 293.15,
    ) -> None:
        check_positive_number("w_value", w_value)
        check_positive_number("density", density)
        check_positive_number("pressure", pressure)
        check_positive_number("temperature", temperature)

        self.name = name
        self.w_value = w_value
        self.density = density
        self.pressure = pressure
        self.temperature = temperature
        self.carriers = carriers
        self.ionisation_products = ionisation_products
        self.reactions = reactions or []

        self._check_carriers()
        self._check_reactions()
        self._check_ionisation_products()

    def _check_carriers(self) -> None:
        names = []
        for carrier in self.carriers:
            if not isinstance(carrier, ChargeCarrier):
                raise TypeError(
                    f"carriers must be ion_chambers.ChargeCarrier objects, "
                    f"not {type(carrier).__name__}"
                )
            if carrier.name in names:
                raise ValueError(
                    f"carriers must have unique names, but '{carrier.name}' is used "
                    f"twice"
                )
            names.append(carrier.name)

    def _check_reactions(self) -> None:
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
                        f"the carriers of the gas"
                    )

    def _check_ionisation_products(self) -> None:
        total_charge = 0
        for carrier in self.ionisation_products:
            if not self._has_carrier(carrier):
                carrier_name = getattr(carrier, "name", carrier)
                raise ValueError(
                    f"ionisation_products contains '{carrier_name}', which is not "
                    f"one of the carriers of the gas"
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
