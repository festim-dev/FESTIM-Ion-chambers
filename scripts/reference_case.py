"""The reference case of the prototype study: a 1 L coaxial chamber in room air.

This file is an example of how to define a gas and a chamber and how to set up the
two problems. The numbers are reasonable starting values, not validated reference
data. Copy it and adapt it to your own gas and conditions.

The other scripts in this folder import it:

    from reference_case import air, build_transport_problem, reference_chamber
"""

import ion_chambers as ic

STANDARD_PRESSURE = 101325.0  # Pa
ROOM_TEMPERATURE = 293.15  # K


def air(pressure=STANDARD_PRESSURE, temperature=ROOM_TEMPERATURE):
    """Room air: electrons, positive cluster ions and O2- negative ions.

    Free electrons attach to O2 within nanoseconds, so nearly all the negative charge
    is carried by slow ions, and ion-ion recombination limits the collection.

    Args:
        pressure: the pressure (Pa)
        temperature: the temperature (K)

    Returns:
        the gas, an ion_chambers.Gas
    """
    # The number density of the gas relative to 1 atm and 293.15 K. The coefficients
    # below are given at those conditions and scaled with this ratio.
    relative_density = (pressure / STANDARD_PRESSURE) * (ROOM_TEMPERATURE / temperature)

    # Mobilities are inversely proportional to the gas density.
    electron_mobility = 0.04 / relative_density  # m2/V/s, at low field
    positive_ion_mobility = 2.0e-4 / relative_density  # m2/V/s
    negative_ion_mobility = 2.2e-4 / relative_density  # m2/V/s

    # Electrons are hotter than the gas, so their diffusion coefficient does not follow
    # the Einstein relation. It is the mobility times the characteristic energy.
    electron_characteristic_energy = 1.0  # V
    electron_diffusivity = electron_mobility * electron_characteristic_energy  # m2/s

    # Attachment to O2 needs a third molecule to carry away the energy, so its
    # frequency goes with the square of the gas density.
    attachment_frequency = 1.0e8 * relative_density**2  # 1/s

    # Kept constant here. This is a simplification: neither coefficient is constant
    # over a wide pressure range.
    ion_ion_recombination = 2.0e-12  # m3/s
    electron_ion_recombination = 2.0e-13  # m3/s

    electrons = ic.ChargeCarrier(
        name="e",
        charge=-1,
        mobility=electron_mobility,
        diffusivity=electron_diffusivity,
    )
    positive_ions = ic.ChargeCarrier(name="P", charge=1, mobility=positive_ion_mobility)
    negative_ions = ic.ChargeCarrier(
        name="N", charge=-1, mobility=negative_ion_mobility
    )

    # electron + O2 -> negative ion
    attachment = ic.CarrierReaction(
        reactant=electrons, product=negative_ions, rate=attachment_frequency
    )
    # positive ion + negative ion -> neutral
    ion_ion = ic.CarrierReaction(
        reactant=[positive_ions, negative_ions], rate=ion_ion_recombination
    )
    # positive ion + electron -> neutral
    electron_ion = ic.CarrierReaction(
        reactant=[positive_ions, electrons], rate=electron_ion_recombination
    )

    return ic.Gas(
        name="air",
        w_value=33.97,  # eV per ion pair, dry air
        density=1.205 * relative_density,  # kg/m3
        pressure=pressure,
        temperature=temperature,
        carriers=[electrons, positive_ions, negative_ions],
        reactions=[attachment, ion_ion, electron_ion],
        ionisation_products=(electrons, positive_ions),
    )


def reference_chamber():
    """A generic 1 L coaxial chamber. It is not the geometry of a real monitor.

    A 1.5 mm radius rod in a 5 cm radius can. The length, about 12.7 cm, follows from
    the volume.

    Returns:
        the geometry, an ion_chambers.CoaxialChamber
    """
    return ic.CoaxialChamber.from_volume(
        volume=1e-3, rod_radius=1.5e-3, can_radius=0.05
    )


def build_transport_problem(chamber, bias, resolution):
    """Solves the potential of a chamber and sets up the transport problem in air.

    Args:
        chamber: the geometry, an ion_chambers.CoaxialChamber
        bias: the potential of the anode relative to the cathode (V)
        resolution: the number of mesh cells across the electrode gap

    Returns:
        the transport problem, an ion_chambers.ChargeTransportProblem that is not
        yet initialised. It has no exports and its tritium activity is zero: give
        it exports before initialising, and set the activity before running.
    """
    # 1. the electric potential
    electrostatics = ic.ElectrostaticProblem()
    electrostatics.mesh = chamber.create_mesh(resolution=resolution)
    electrostatics.subdomains = [chamber.gas, chamber.anode, chamber.cathode]
    electrostatics.boundary_conditions = [
        ic.FixedPotentialBC(subdomain=chamber.anode, value=bias),
        ic.FixedPotentialBC(subdomain=chamber.cathode, value=0.0),
    ]
    electrostatics.initialise()
    electrostatics.run()

    # 2. the transport of the charge carriers in that potential
    transport = ic.ChargeTransportProblem()
    transport.electrostatic_problem = electrostatics
    transport.gas = air()
    transport.ionisation_sources = [ic.TritiumInGas(activity=0.0, geometry=chamber)]
    return transport
