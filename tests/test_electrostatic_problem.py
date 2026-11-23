import festim as F
import numpy as np
import pytest
from dolfinx import fem
from prototype_case import solve_potential

import ion_chambers as ic
from ion_chambers.charge_transport_problem import cell_field_length
from ion_chambers.constants import EPS0


@pytest.fixture
def long_chamber() -> ic.CoaxialChamber:
    """Long enough for the field at mid-height to be that of an infinite coax."""
    return ic.CoaxialChamber(rod_radius=1e-3, can_radius=0.02, length=0.4)


def test_is_a_festim_problem():
    """An ElectrostaticProblem is a FESTIM problem, like HeatTransferProblem."""
    assert isinstance(ic.ElectrostaticProblem(), F.ProblemBase)


def test_fixed_potential_bc_is_a_festim_dirichlet_bc(chamber):
    """A FixedPotentialBC is a FESTIM Dirichlet boundary condition."""
    bc = ic.FixedPotentialBC(subdomain=chamber.anode, value=500.0)
    assert isinstance(bc, F.DirichletBCBase)


def test_rejects_boundary_conditions_of_another_type(chamber):
    """A boundary condition that is not a FixedPotentialBC raises a TypeError."""
    temperature_bc = F.FixedTemperatureBC(subdomain=chamber.anode, value=500.0)
    with pytest.raises(TypeError, match="boundary_conditions must be a list of"):
        ic.ElectrostaticProblem(boundary_conditions=[temperature_bc])


def test_potential_is_none_before_initialise():
    """The potential does not exist until the problem is initialised."""
    assert ic.ElectrostaticProblem().potential is None


def test_potential_is_the_solution(long_chamber):
    """After a run, potential is the solution function u."""
    electrostatics = solve_potential(long_chamber, resolution=10, bias=300.0)
    assert electrostatics.potential is electrostatics.u


def test_potential_equals_the_bias_on_the_anode(long_chamber):
    """The largest value of the potential is the bias, reached on the anode."""
    electrostatics = solve_potential(long_chamber, resolution=20, bias=300.0)
    assert np.isclose(electrostatics.u.x.array.max(), 300.0)


def test_potential_is_zero_on_the_cathode(long_chamber):
    """The smallest value of the potential is zero, reached on the cathode."""
    electrostatics = solve_potential(long_chamber, resolution=20, bias=300.0)
    assert np.isclose(electrostatics.u.x.array.min(), 0.0, atol=1e-9)


def test_field_at_mid_height_is_that_of_an_infinite_coax(long_chamber):
    """Far from the end caps the field matches V / (r ln(R/a)) within 0.2%."""
    electrostatics = solve_potential(long_chamber, resolution=40, bias=300.0)
    strength, _ = cell_field_length(electrostatics.u)

    mesh = electrostatics.mesh.mesh
    centroids = mesh.geometry.x[:, :2][mesh.geometry.dofmaps[0]].mean(axis=1)
    middle = np.abs(centroids[:, 1] - 0.2) < 0.02
    a, R = long_chamber.rod_radius, long_chamber.can_radius
    exact = 300.0 / (centroids[middle, 0] * np.log(R / a))
    assert np.allclose(strength[middle], exact, rtol=2e-3)


def test_cell_length_is_bounded_by_the_cell_size(long_chamber):
    """The cell length along the field never exceeds twice the longest cell edge."""
    electrostatics = solve_potential(long_chamber, resolution=10, bias=300.0)
    _, length = cell_field_length(electrostatics.u)
    mesh = electrostatics.mesh.mesh
    vertices = mesh.geometry.x[:, :2][mesh.geometry.dofmaps[0]]
    longest_edge = np.linalg.norm(vertices[:, 1] - vertices[:, 0], axis=1).max()
    assert length.max() <= 2 * longest_edge


def test_uniform_charge_between_grounded_electrodes(long_chamber):
    """A uniform charge between grounded electrodes of an infinite coax gives the
    analytic potential at mid-height, within 1%."""
    charge_density = 1e-6  # C/m3
    mesh = long_chamber.create_mesh(resolution=40)
    electrostatics = ic.ElectrostaticProblem(
        mesh=mesh,
        subdomains=[long_chamber.gas, long_chamber.anode, long_chamber.cathode],
        boundary_conditions=[
            ic.FixedPotentialBC(subdomain=long_chamber.anode, value=0.0),
            ic.FixedPotentialBC(subdomain=long_chamber.cathode, value=0.0),
        ],
        charge_density=fem.Constant(mesh.mesh, charge_density),
    )
    electrostatics.initialise()
    electrostatics.run()

    # phi(r) = k (R^2 - r^2) - k (R^2 - a^2) ln(R/r) / ln(R/a), k = rho / (4 eps0)
    a, R = long_chamber.rod_radius, long_chamber.can_radius
    r = 0.01
    k = charge_density / (4 * EPS0)
    exact = k * (R**2 - r**2) - k * (R**2 - a**2) * np.log(R / r) / np.log(R / a)

    point = np.array([[r, 0.2, 0.0]])
    cell = _cell_containing(mesh.mesh, point)
    assert np.isclose(electrostatics.u.eval(point, cell)[0], exact, rtol=1e-2)


def _cell_containing(mesh, point: np.ndarray) -> np.ndarray:
    """The index of the cell that contains a point."""
    from dolfinx import geometry

    tree = geometry.bb_tree(mesh, mesh.topology.dim)
    candidates = geometry.compute_collisions_points(tree, point)
    colliding = geometry.compute_colliding_cells(mesh, candidates, point)
    return colliding.links(0)[:1]


def test_potential_export_writes_a_file(long_chamber, tmp_path):
    """A PotentialExport writes the potential to the file it was given."""
    filename = tmp_path / "potential.bp"
    electrostatics = ic.ElectrostaticProblem(
        mesh=long_chamber.create_mesh(resolution=6),
        subdomains=[long_chamber.gas, long_chamber.anode, long_chamber.cathode],
        boundary_conditions=[
            ic.FixedPotentialBC(subdomain=long_chamber.anode, value=100.0),
            ic.FixedPotentialBC(subdomain=long_chamber.cathode, value=0.0),
        ],
        exports=[ic.PotentialExport(filename=filename)],
    )
    electrostatics.initialise()
    electrostatics.run()
    assert filename.exists()


def test_rejects_exports_of_another_type(long_chamber):
    """An export that is not a PotentialExport raises a TypeError at initialise."""
    electrostatics = ic.ElectrostaticProblem(
        mesh=long_chamber.create_mesh(resolution=6),
        subdomains=[long_chamber.gas, long_chamber.anode, long_chamber.cathode],
        boundary_conditions=[
            ic.FixedPotentialBC(subdomain=long_chamber.anode, value=100.0),
            ic.FixedPotentialBC(subdomain=long_chamber.cathode, value=0.0),
        ],
        exports=[ic.IdealCurrent()],
    )
    with pytest.raises(TypeError, match="exports must be a list of"):
        electrostatics.initialise()
