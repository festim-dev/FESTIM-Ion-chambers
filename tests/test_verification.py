"""Verification: the transport problem against known answers and under mesh refinement.

Three cases are used. A parallel-plate chamber with a uniform source and two kinds of
ions has an analytical collection efficiency (Boag). A manufactured solution between
two electrodes checks the densities themselves, on a cartesian and on a cylindrical
mesh. The reference coaxial chamber at an activity too low for recombination must
collect every pair that is created.
"""

from mpi4py import MPI

import dolfinx
import festim as F
import numpy as np
import pytest
import ufl
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
        stepsize=0.003,
    )
    transient = solve_parallel_plates(cells=50, pair_rate=1e14, settings=settings)
    steady = solve_parallel_plates(cells=50, pair_rate=1e14)
    assert np.isclose(transient["current"].data[-1], steady["current"].value, rtol=1e-6)


# ---- a manufactured solution ---------------------------------------------------------
# Densities are chosen, the sources that make them a solution of the equations are
# worked out from them, and the solve must give the chosen densities back. The gas is
# held between an electrode at x = INNER and one at x = OUTER: two plates on a
# cartesian mesh, two coaxial cylinders on a cylindrical one.

INNER = 0.01  # m
OUTER = 0.05  # m
MMS_BIAS = 500.0  # V, on the inner electrode
MMS_ROWS = 4  # cells along the electrodes
REFERENCE_DENSITY = 1e14  # 1/m3


def chosen_positive_density(x):
    position = (x[0] - INNER) / (OUTER - INNER)  # 0 on the inner electrode, 1 outer
    return REFERENCE_DENSITY * (1 + 0.5 * ufl.sin(ufl.pi * position))


def chosen_negative_density(x):
    position = (x[0] - INNER) / (OUTER - INNER)
    return REFERENCE_DENSITY * (1 + 0.5 * ufl.cos(ufl.pi * position))


def manufactured_error(
    cells: int,
    coordinate_system: str,
    positive_diffusivity: float,
    negative_diffusivity: float,
    electrode_conditions: str = "fixed density",
) -> float:
    """Solves for the chosen densities and returns the larger of the two relative L2
    errors.

    The chosen densities vary between the electrodes only, so the mesh is a strip
    MMS_ROWS cells high. Its height shrinks with the cell size, which keeps the
    cells square: the errors are those of a mesh as high as it is wide.

    Args:
        cells: the number of cells between the electrodes
        coordinate_system: "cartesian" or "cylindrical"
        positive_diffusivity: the diffusion coefficient of the positive ions (m2/s)
        negative_diffusivity: the diffusion coefficient of the negative ions (m2/s)
        electrode_conditions: "fixed density" sets the chosen densities on the
            electrodes. "absorbing" puts an AbsorbingElectrodeBC there instead, with
            a flux that makes up the difference to the chosen densities.
    """
    is_cylindrical = coordinate_system == "cylindrical"

    cell_size = (OUTER - INNER) / cells
    mesh = dolfinx.mesh.create_rectangle(
        MPI.COMM_WORLD,
        [[INNER, 0.0], [OUTER, MMS_ROWS * cell_size]],
        [cells, MMS_ROWS],
        dolfinx.mesh.CellType.triangle,
        diagonal=dolfinx.mesh.DiagonalType.crossed,
    )
    gas_volume = F.VolumeSubdomain(
        id=1, material=None, locator=lambda x: np.full(x.shape[1], True)
    )
    inner = F.SurfaceSubdomain(id=1, locator=lambda x: np.isclose(x[0], INNER))
    outer = F.SurfaceSubdomain(id=2, locator=lambda x: np.isclose(x[0], OUTER))

    electrostatics = ic.ElectrostaticProblem()
    electrostatics.mesh = F.Mesh(mesh, coordinate_system=coordinate_system)
    electrostatics.subdomains = [gas_volume, inner, outer]
    electrostatics.boundary_conditions = [
        ic.FixedPotentialBC(subdomain=inner, value=MMS_BIAS),
        ic.FixedPotentialBC(subdomain=outer, value=0.0),
    ]
    electrostatics.initialise()
    electrostatics.run()

    positive_ions = ic.ChargeCarrier(
        name="P", charge=1, mobility=POSITIVE_MOBILITY, diffusivity=positive_diffusivity
    )
    negative_ions = ic.ChargeCarrier(
        name="N",
        charge=-1,
        mobility=NEGATIVE_MOBILITY,
        diffusivity=negative_diffusivity,
    )
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

    # The exact field, along x and pointing from the inner to the outer electrode:
    # uniform between plates, as 1/r between cylinders.
    def field(x):
        if is_cylindrical:
            return MMS_BIAS / (x[0] * np.log(OUTER / INNER))
        return MMS_BIAS / (OUTER - INNER)

    # the divergence of a flux along x that depends on x only
    def divergence(flux, x):
        if is_cylindrical:
            return flux.dx(0) + flux / x[0]
        return flux.dx(0)

    def recombined(x):
        density_product = chosen_positive_density(x) * chosen_negative_density(x)
        return RECOMBINATION_RATE * density_product

    # the drift velocity and the flux of each carrier, both along x
    def positive_velocity(x):
        return POSITIVE_MOBILITY * field(x)

    def negative_velocity(x):
        return -NEGATIVE_MOBILITY * field(x)

    def positive_flux(x):
        density = chosen_positive_density(x)
        return -positive_diffusivity * density.dx(0) + density * positive_velocity(x)

    def negative_flux(x):
        density = chosen_negative_density(x)
        return -negative_diffusivity * density.dx(0) + density * negative_velocity(x)

    # the source of each carrier: the divergence of its flux, plus what recombines
    def positive_source(x):
        return divergence(positive_flux(x), x) + recombined(x)

    def negative_source(x):
        return divergence(negative_flux(x), x) + recombined(x)

    # No ionisation source: the sources and the conditions on the electrodes are
    # given as FESTIM objects. The chosen densities vary along x only, so the two
    # other sides need no condition.
    transport = ic.ChargeTransportProblem(electrostatic_problem=electrostatics, gas=gas)
    transport.sources = [
        F.ParticleSource(
            value=positive_source, volume=gas_volume, species=positive_ions
        ),
        F.ParticleSource(
            value=negative_source, volume=gas_volume, species=negative_ions
        ),
    ]

    # for each carrier: its chosen density, its velocity and its flux
    carriers = [
        (positive_ions, chosen_positive_density, positive_velocity, positive_flux),
        (negative_ions, chosen_negative_density, negative_velocity, negative_flux),
    ]
    # for each electrode: the x component of the normal pointing out of the gas
    electrodes = [(inner, -1.0), (outer, 1.0)]

    for electrode, normal in electrodes:
        for carrier, chosen_density, velocity, flux in carriers:
            if electrode_conditions == "fixed density":
                fixed_density = F.FixedConcentrationBC(
                    subdomain=electrode, value=chosen_density, species=carrier
                )
                transport.boundary_conditions.append(fixed_density)
                continue

            # An absorbing electrode lets c max(v.n, 0) out of the gas. The chosen
            # densities have another flux there, so the difference is brought in.
            def difference_to_absorbing(
                x,
                normal=normal,
                chosen_density=chosen_density,
                velocity=velocity,
                flux=flux,
            ):
                absorbed = chosen_density(x) * ufl.max_value(velocity(x) * normal, 0.0)
                return absorbed - flux(x) * normal

            absorbing = ic.AbsorbingElectrodeBC(subdomain=electrode, species=carrier)
            make_up = F.ParticleFluxBC(
                subdomain=electrode, value=difference_to_absorbing, species=carrier
            )
            transport.boundary_conditions.append(absorbing)
            transport.boundary_conditions.append(make_up)

    transport.initialise()
    transport.run()

    x = ufl.SpatialCoordinate(mesh)
    weight = electrostatics.mesh.coordinate_system.integration_weight(mesh)

    def relative_error(carrier, chosen_density) -> float:
        difference = carrier.post_processing_solution - chosen_density(x)
        squared_error = dolfinx.fem.assemble_scalar(
            dolfinx.fem.form(difference**2 * weight * ufl.dx)
        )
        squared_norm = dolfinx.fem.assemble_scalar(
            dolfinx.fem.form(chosen_density(x) ** 2 * weight * ufl.dx)
        )
        return np.sqrt(squared_error / squared_norm)

    positive_error = relative_error(positive_ions, chosen_positive_density)
    negative_error = relative_error(negative_ions, chosen_negative_density)
    return max(positive_error, negative_error)


@pytest.mark.parametrize("coordinate_system", ["cartesian", "cylindrical"])
def test_manufactured_solution_converges_at_second_order_without_the_floor(
    coordinate_system,
):
    """With diffusion strong enough that the floor does nothing, halving the cell
    size divides the error of the densities by four."""
    coarse = manufactured_error(
        cells=16,
        coordinate_system=coordinate_system,
        positive_diffusivity=0.02,
        negative_diffusivity=0.03,
    )
    fine = manufactured_error(
        cells=32,
        coordinate_system=coordinate_system,
        positive_diffusivity=0.02,
        negative_diffusivity=0.03,
    )
    assert np.isclose(np.log2(coarse / fine), 2.0, atol=0.05)


@pytest.mark.parametrize("coordinate_system", ["cartesian", "cylindrical"])
def test_manufactured_solution_converges_at_first_order_with_the_floor(
    coordinate_system,
):
    """With the weak diffusion of real ions the floor sets the diffusion, and
    halving the cell size halves the error of the densities."""
    coarse = manufactured_error(
        cells=64,
        coordinate_system=coordinate_system,
        positive_diffusivity=5.0e-6,
        negative_diffusivity=5.5e-6,
    )
    fine = manufactured_error(
        cells=128,
        coordinate_system=coordinate_system,
        positive_diffusivity=5.0e-6,
        negative_diffusivity=5.5e-6,
    )
    assert np.isclose(np.log2(coarse / fine), 1.0, atol=0.1)


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


@pytest.mark.parametrize("coordinate_system", ["cartesian", "cylindrical"])
def test_absorbing_electrodes_converge_at_second_order_without_the_floor(
    coordinate_system,
):
    """With absorbing electrodes in place of fixed densities, and diffusion strong
    enough that the floor does nothing, halving the cell size divides the error of
    the densities by four. The meshes are finer than for fixed densities: on the
    cylindrical mesh the order is 1.92 between 16 and 32 cells and rises to 2,
    because the drift velocity on the electrode comes from the solved potential."""
    coarse = manufactured_error(
        cells=64,
        coordinate_system=coordinate_system,
        positive_diffusivity=0.02,
        negative_diffusivity=0.03,
        electrode_conditions="absorbing",
    )
    fine = manufactured_error(
        cells=128,
        coordinate_system=coordinate_system,
        positive_diffusivity=0.02,
        negative_diffusivity=0.03,
        electrode_conditions="absorbing",
    )
    assert np.isclose(np.log2(coarse / fine), 2.0, atol=0.05)


@pytest.mark.parametrize("coordinate_system", ["cartesian", "cylindrical"])
def test_absorbing_electrodes_converge_at_first_order_with_the_floor(
    coordinate_system,
):
    """With absorbing electrodes and the weak diffusion of real ions, halving the
    cell size halves the error of the densities."""
    coarse = manufactured_error(
        cells=64,
        coordinate_system=coordinate_system,
        positive_diffusivity=5.0e-6,
        negative_diffusivity=5.5e-6,
        electrode_conditions="absorbing",
    )
    fine = manufactured_error(
        cells=128,
        coordinate_system=coordinate_system,
        positive_diffusivity=5.0e-6,
        negative_diffusivity=5.5e-6,
        electrode_conditions="absorbing",
    )
    assert np.isclose(np.log2(coarse / fine), 1.0, atol=0.1)
