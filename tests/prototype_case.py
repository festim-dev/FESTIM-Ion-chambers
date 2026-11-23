"""The case the original prototype was run for: a 1 L coaxial chamber in room air.

The regression tests pin currents computed with exactly these numbers and this mesh,
so do not change them.
"""

from mpi4py import MPI

import festim as F
import numpy as np
from dolfinx import mesh as dmesh

import ion_chambers as ic


def air() -> ic.Gas:
    """Room air at 1 atm and 293.15 K, with the coefficients of the prototype."""
    electrons = ic.ChargeCarrier(name="e", charge=-1, mobility=0.04, diffusivity=0.04)
    positive_ions = ic.ChargeCarrier(name="P", charge=1, mobility=2.0e-4)
    negative_ions = ic.ChargeCarrier(name="N", charge=-1, mobility=2.2e-4)

    attachment = ic.CarrierReaction(
        reactant=electrons, product=negative_ions, rate=1.0e8
    )
    ion_ion_recombination = ic.CarrierReaction(
        reactant=[positive_ions, negative_ions], rate=2.0e-12
    )
    electron_ion_recombination = ic.CarrierReaction(
        reactant=[positive_ions, electrons], rate=2.0e-13
    )

    return ic.Gas(
        name="air",
        w_value=33.97,
        density=1.205,
        carriers=[electrons, positive_ions, negative_ions],
        ionisation_products=(electrons, positive_ions),
        reactions=[attachment, ion_ion_recombination, electron_ion_recombination],
    )


def reference_chamber() -> ic.CoaxialChamber:
    """A 1.5 mm radius rod in a 5 cm radius can holding 1 L of gas."""
    return ic.CoaxialChamber.from_volume(
        volume=1e-3, rod_radius=1.5e-3, can_radius=0.05
    )


def prototype_mesh(
    chamber: ic.CoaxialChamber, radial_cells: int, axial_cells: int
) -> F.Mesh:
    """The graded mesh of the original prototype, which the pinned currents were
    computed on: the radial spacing grows in proportion to r and the axial spacing
    is smallest at the end caps."""
    mesh = dmesh.create_rectangle(
        MPI.COMM_WORLD,
        [[0, -1], [1, 1]],
        [radial_cells, axial_cells],
        dmesh.CellType.triangle,
        diagonal=dmesh.DiagonalType.crossed,
    )
    x = mesh.geometry.x
    s = x[:, 0].copy()  # 0 on the rod, 1 on the can
    t = x[:, 1].copy()  # -1 on the bottom cap, 1 on the top cap
    x[:, 0] = chamber.rod_radius * (chamber.can_radius / chamber.rod_radius) ** s
    x[:, 1] = 0.5 * chamber.length * (1 + np.tanh(1.6 * t) / np.tanh(1.6))
    return F.Mesh(mesh, coordinate_system="cylindrical")


def solve_potential(
    chamber: ic.CoaxialChamber,
    radial_cells: int,
    axial_cells: int,
    bias: float = 500.0,
) -> ic.ElectrostaticProblem:
    """The potential of a chamber with its anode at ``bias`` and a grounded cathode."""
    electrostatics = ic.ElectrostaticProblem()
    electrostatics.mesh = prototype_mesh(
        chamber, radial_cells=radial_cells, axial_cells=axial_cells
    )
    electrostatics.subdomains = [chamber.gas, chamber.anode, chamber.cathode]
    electrostatics.boundary_conditions = [
        ic.FixedPotentialBC(subdomain=chamber.anode, value=bias),
        ic.FixedPotentialBC(subdomain=chamber.cathode, value=0.0),
    ]
    electrostatics.initialise()
    electrostatics.run()
    return electrostatics
