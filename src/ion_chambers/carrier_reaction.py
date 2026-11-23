"""A reaction between charge carriers in the gas."""

from __future__ import annotations

import festim as F

from ion_chambers.charge_carrier import ChargeCarrier
from ion_chambers.coefficient import Coefficient, check_coefficient


class CarrierReaction(F.GenericReaction):
    """A volume reaction between charge carriers, with a mass-action rate.

    A :class:`festim.GenericReaction` that differs from its parent in two ways:

    * It needs no volume subdomain when it is created. The model creates the gas
      volume during ``initialise`` and assigns it then.
    * Its rate coefficient may be a callable of the reduced electric field E/N, which
      FESTIM's rate coefficients do not know about. The model evaluates ``rate`` in
      the local field during ``initialise`` and stores the result as
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
        forward_rate: the rate coefficient FESTIM uses, set by the model
        volume: the gas volume, set by the model
    """

    def __init__(
        self,
        reactant: ChargeCarrier | list[ChargeCarrier],
        rate: Coefficient,
        product: ChargeCarrier | list[ChargeCarrier] | None = None,
    ) -> None:
        check_coefficient("rate", rate)
        self.rate = rate
        # volume and forward_rate are placeholders until the model sets them
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
