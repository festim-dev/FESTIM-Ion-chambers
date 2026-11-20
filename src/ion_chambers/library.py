"""Example gases.

These are starting points with documented, simple assumptions, not reference data.
Check the coefficients against the conditions of your experiment, or build your own
:class:`ion_chambers.Gas`.
"""

from ion_chambers.constants import ROOM_TEMPERATURE, STANDARD_PRESSURE
from ion_chambers.gas import ChargeCarrier, Gas, Reaction
from ion_chambers.geometry import CoaxialChamber


def air(
    pressure: float = STANDARD_PRESSURE, temperature: float = ROOM_TEMPERATURE
) -> Gas:
    """Room air: electrons, positive cluster ions and O2- negative ions.

    Free electrons attach to O2 within nanoseconds, so nearly all the negative charge
    is carried by slow ions and ion-ion recombination limits the collection.

    Values at 1 atm and 293.15 K:

    ====================================  ==========  =======
    quantity                              value       unit
    ====================================  ==========  =======
    W value                               33.97       eV
    density                               1.205       kg/m3
    positive ion mobility                 2.0e-4      m2/V/s
    negative ion mobility                 2.2e-4      m2/V/s
    electron mobility (low field)         0.04        m2/V/s
    electron characteristic energy D/mu   1.0         V
    attachment frequency to O2            1.0e8       1/s
    ion-ion recombination                 2.0e-12     m3/s
    electron-ion recombination            2.0e-13     m3/s
    ====================================  ==========  =======

    Away from these conditions the coefficients follow the relative number density
    ``n = (p / 1 atm) (293.15 K / T)``: the density scales as ``n``, the mobilities as
    ``1 / n`` and the three-body attachment frequency as ``n**2``. The recombination
    coefficients and the electron characteristic energy are kept constant, which is a
    simplification: neither holds over a wide pressure range.

    Args:
        pressure: the pressure (Pa)
        temperature: the temperature (K)

    Returns:
        the gas
    """
    n = (pressure / STANDARD_PRESSURE) * (ROOM_TEMPERATURE / temperature)
    electron_mobility = 0.04 / n
    characteristic_energy = 1.0
    return Gas(
        name="air",
        w_value=33.97,
        density=1.205 * n,
        pressure=pressure,
        temperature=temperature,
        carriers=[
            ChargeCarrier(
                name="e",
                charge=-1,
                mobility=electron_mobility,
                diffusivity=electron_mobility * characteristic_energy,
            ),
            ChargeCarrier(name="P", charge=1, mobility=2.0e-4 / n),
            ChargeCarrier(name="N", charge=-1, mobility=2.2e-4 / n),
        ],
        reactions=[
            Reaction(reactants=["e"], products=["N"], rate=1.0e8 * n**2),
            Reaction(reactants=["P", "N"], rate=2.0e-12),
            Reaction(reactants=["P", "e"], rate=2.0e-13),
        ],
    )


def reference_chamber() -> CoaxialChamber:
    """The generic 1 L coaxial chamber of the prototype study.

    A 1.5 mm radius rod in a 5 cm radius can, about 12.7 cm long. It is *not* the
    geometry of a real monitor: it is the reference case the prototype results in the
    documentation were computed for.
    """
    return CoaxialChamber.from_volume(volume=1e-3, rod_radius=1.5e-3, can_radius=0.05)
