"""A charged species that drifts in the electric field."""

from __future__ import annotations

import festim as F

from ion_chambers.coefficient import Coefficient, check_coefficient


class ChargeCarrier(F.Species):
    """A charged species transported by drift and diffusion.

    A :class:`festim.Species` with a charge and a mobility. It can be used wherever
    FESTIM expects a species.

    A carrier holds the solution of the model it belongs to, like any FESTIM species,
    so do not share one carrier between two models that are alive at the same time.

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
