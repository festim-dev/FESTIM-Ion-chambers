"""Verification: the transport problem against known answers and under mesh refinement.

Two cases are used. A parallel-plate chamber with a uniform source and two kinds of
ions has an analytical collection efficiency (Boag). The reference coaxial chamber at
an activity too low for recombination must collect every pair that is created.
"""

from mpi4py import MPI

import dolfinx
import festim as F
import numpy as np
import pytest
from prototype_case import air, reference_chamber, solve_potential
from scipy.constants import elementary_charge

import ion_chambers as ic

pytestmark = pytest.mark.slow

# ---- the parallel-plate chamber -----------------------------------------------------

GAP = 0.01  # m, distance between the plates
HEIGHT = 0.002  # m, extent of the mesh along the plates
BIAS = 100.0  # V
POSITIVE_MOBILITY = 2.0e-4  # m2/V/s
NEGATIVE_MOBILITY = 2.2e-4  # m2/V/s
RECOMBINATION_RATE = 2.0e-12  # m3/s
W_VALUE = 33.97  # eV


def boag_collection_efficiency(pair_rate: float) -> float:
    """The collection efficiency of a parallel-plate chamber under continuous,
    uniform irradiation, for a small loss to recombination.

    J. W. Boag, in Radiation Dosimetry vol. 2 (Attix and Roesch, 1966):
    f = 1 / (1 + m^2 d^4 q / (6 V^2)) with m^2 = alpha / (e mu+ mu-) and q the
    charge of one sign created per volume and time, here e times the pair rate.

    Args:
        pair_rate: the production rate of ion pairs (1/m3/s)
    """
    loss = (
        RECOMBINATION_RATE
        * pair_rate
        * GAP**4
        / (6 * POSITIVE_MOBILITY * NEGATIVE_MOBILITY * BIAS**2)
    )
    return 1 / (1 + loss)


def solve_parallel_plates(
    cells: int, pair_rate: float, settings: F.Settings | None = None
) -> dict[str, ic.ChamberQuantity]:
    """Solves the parallel-plate chamber and returns its quantities.

    Args:
        cells: the number of cells between the plates
        pair_rate: the production rate of ion pairs (1/m3/s)
        settings: the settings of the transport problem, steady state by default
    """
    mesh = dolfinx.mesh.create_rectangle(
        MPI.COMM_WORLD,
        [[0.0, 0.0], [GAP, HEIGHT]],
        [cells, 2],
        dolfinx.mesh.CellType.triangle,
        diagonal=dolfinx.mesh.DiagonalType.crossed,
    )
    gas_volume = F.VolumeSubdomain(
        id=1, material=None, locator=lambda x: np.full(x.shape[1], True)
    )
    positive_plate = F.SurfaceSubdomain(id=1, locator=lambda x: np.isclose(x[0], 0.0))
    grounded_plate = F.SurfaceSubdomain(id=2, locator=lambda x: np.isclose(x[0], GAP))

    electrostatics = ic.ElectrostaticProblem()
    electrostatics.mesh = F.Mesh(mesh)
    electrostatics.subdomains = [gas_volume, positive_plate, grounded_plate]
    electrostatics.boundary_conditions = [
        ic.FixedPotentialBC(subdomain=positive_plate, value=BIAS),
        ic.FixedPotentialBC(subdomain=grounded_plate, value=0.0),
    ]
    electrostatics.initialise()
    electrostatics.run()

    positive_ions = ic.ChargeCarrier(name="P", charge=1, mobility=POSITIVE_MOBILITY)
    negative_ions = ic.ChargeCarrier(name="N", charge=-1, mobility=NEGATIVE_MOBILITY)
    recombination = ic.CarrierReaction(
        reactant=[positive_ions, negative_ions], rate=RECOMBINATION_RATE
    )
    gas = ic.Gas(
        name="two kinds of ions",
        w_value=W_VALUE,
        density=1.205,
        carriers=[positive_ions, negative_ions],
        ionisation_products=(negative_ions, positive_ions),
        reactions=[recombination],
    )

    quantities = {
        "current": ic.ChamberCurrent(surface=positive_plate),
        "collection_efficiency": ic.CollectionEfficiency(volume=gas_volume),
    }
    activity = pair_rate / ic.beta.ion_pairs_per_decay(W_VALUE)
    transport = ic.ChargeTransportProblem(
        electrostatic_problem=electrostatics,
        gas=gas,
        ionisation_sources=[ic.TritiumInGas(activity=activity)],
        exports=list(quantities.values()),
        settings=settings,
    )
    transport.show_progress_bar = False
    transport.initialise()
    transport.run()
    return quantities


def collected_share(cells: int, pair_rate: float) -> float:
    """The current of the parallel-plate chamber over the current if every pair
    were collected. The mesh is a unit deep, so the gas volume is GAP * HEIGHT."""
    quantities = solve_parallel_plates(cells=cells, pair_rate=pair_rate)
    ideal_current = elementary_charge * pair_rate * GAP * HEIGHT
    return quantities["current"].value / ideal_current


@pytest.mark.parametrize("pair_rate", [1e14, 1e15])
def test_recombination_loss_matches_boag(pair_rate):
    """The share of the pairs lost to recombination between parallel plates agrees
    with Boag's formula to 1%."""
    quantities = solve_parallel_plates(cells=200, pair_rate=pair_rate)
    loss = 1 - quantities["collection_efficiency"].value
    expected = 1 - boag_collection_efficiency(pair_rate)
    assert np.isclose(loss, expected, rtol=0.01)


def test_parallel_plate_current_converges_at_first_order():
    """Halving the cell size halves the change in the current: the error of the
    diffusion floor is first order in the cell size."""
    coarse = collected_share(cells=50, pair_rate=1e14)
    medium = collected_share(cells=100, pair_rate=1e14)
    fine = collected_share(cells=200, pair_rate=1e14)
    order = np.log2((medium - coarse) / (fine - medium))
    assert np.isclose(order, 1.0, atol=0.05)


def test_extrapolated_parallel_plate_current_matches_boag():
    """The current extrapolated to zero cell size from two meshes is the current
    of Boag's formula."""
    medium = collected_share(cells=100, pair_rate=1e14)
    fine = collected_share(cells=200, pair_rate=1e14)
    extrapolated = 2 * fine - medium
    assert np.isclose(extrapolated, boag_collection_efficiency(1e14), rtol=1e-4)


def test_transient_current_reaches_the_steady_current():
    """A transient run from an empty chamber, several ion transit times long (one is
    5 ms), ends at the current of the steady-state problem."""
    settings = F.Settings(
        atol=1e-30,
        rtol=1e-10,
        max_iterations=40,
        transient=True,
        final_time=0.03,
        stepsize=0.001,
    )
    transient = solve_parallel_plates(cells=50, pair_rate=1e14, settings=settings)
    steady = solve_parallel_plates(cells=50, pair_rate=1e14)
    assert np.isclose(transient["current"].data[-1], steady["current"].value, rtol=1e-6)


# ---- the reference coaxial chamber --------------------------------------------------

RADIAL_CELLS = [30, 60, 120]  # each mesh twice as fine as the one before


@pytest.fixture(scope="module")
def coaxial_collected_shares() -> list[float]:
    """The current over the ideal current of the reference chamber, without wall
    losses and at 1 uCi/m3 where nothing recombines, on each mesh of RADIAL_CELLS."""
    chamber = reference_chamber()
    shares = []
    for radial_cells in RADIAL_CELLS:
        current = ic.ChamberCurrent(surface=chamber.anode)
        ideal_current = ic.IdealCurrent(volume=chamber.gas)
        transport = ic.ChargeTransportProblem(
            electrostatic_problem=solve_potential(
                chamber, radial_cells=radial_cells, axial_cells=radial_cells * 3 // 5
            ),
            gas=air(),
            ionisation_sources=[ic.TritiumInGas(activity=3.7e4)],
            exports=[current, ideal_current],
        )
        transport.initialise()
        transport.run()
        shares.append(current.value / ideal_current.value)
    return shares


def test_coaxial_current_converges_at_first_order(coaxial_collected_shares):
    """Halving the cell size of the coaxial chamber roughly halves the change in
    the current."""
    coarse, medium, fine = coaxial_collected_shares
    order = np.log2((medium - coarse) / (fine - medium))
    assert np.isclose(order, 1.0, atol=0.1)


def test_extrapolated_coaxial_current_is_the_ideal_current(coaxial_collected_shares):
    """With no wall loss and no recombination, the current extrapolated to zero cell
    size is the charge of all the pairs created per second."""
    medium = coaxial_collected_shares[1]
    fine = coaxial_collected_shares[2]
    extrapolated = 2 * fine - medium
    assert np.isclose(extrapolated, 1.0, rtol=2e-3)
