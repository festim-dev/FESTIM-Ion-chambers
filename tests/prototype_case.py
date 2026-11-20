"""The case the original prototype was run for: a 1 L coaxial chamber in room air.

The regression tests in test_model.py pin currents computed with exactly these
numbers, so do not change them.
"""

import ion_chambers as ic


def air() -> ic.Gas:
    """Room air at 1 atm and 293.15 K, with the coefficients of the prototype."""
    electrons = ic.ChargeCarrier(name="e", charge=-1, mobility=0.04, diffusivity=0.04)
    positive_ions = ic.ChargeCarrier(name="P", charge=1, mobility=2.0e-4)
    negative_ions = ic.ChargeCarrier(name="N", charge=-1, mobility=2.2e-4)

    attachment = ic.Reaction(reactants=["e"], products=["N"], rate=1.0e8)
    ion_ion_recombination = ic.Reaction(reactants=["P", "N"], rate=2.0e-12)
    electron_ion_recombination = ic.Reaction(reactants=["P", "e"], rate=2.0e-13)

    return ic.Gas(
        name="air",
        w_value=33.97,
        density=1.205,
        carriers=[electrons, positive_ions, negative_ions],
        reactions=[attachment, ion_ion_recombination, electron_ion_recombination],
    )


def reference_chamber() -> ic.CoaxialChamber:
    """A 1.5 mm radius rod in a 5 cm radius can holding 1 L of gas."""
    return ic.CoaxialChamber.from_volume(
        volume=1e-3, rod_radius=1.5e-3, can_radius=0.05
    )
