"""Stabilisation, projection, potential solver and integration helpers."""

import numpy as np
import pytest
import ufl
from dolfinx import fem

import ion_chambers as ic
from ion_chambers.integration import assemble, volume_weight
from ion_chambers.projection import lumped_projection
from ion_chambers.stabilisation import cell_field_length


@pytest.fixture
def long_chamber() -> ic.CoaxialChamber:
    """Long enough for the field at mid-height to be that of an infinite coax."""
    return ic.CoaxialChamber(rod_radius=1e-3, can_radius=0.02, length=0.4)


@pytest.mark.parametrize(
    "diffusivity, expected",
    [(1e-6, 5e-4), (1.0, 1.0)],  # floor active, floor inactive
)
def test_floor_is_the_larger_of_physical_and_upwind(diffusivity, expected):
    floor = ic.DiffusionFloor()
    D = floor.apply(
        diffusivity=np.array([diffusivity]),
        mobility=np.array([1e-4]),
        field_strength=np.array([1e4]),
        cell_length=np.array([1e-3]),
    )
    assert np.isclose(D[0], expected)


def test_floor_factor_zero_switches_it_off():
    D = ic.DiffusionFloor(factor=0.0).apply(
        np.array([1e-6]), np.array([1e-4]), np.array([1e4]), np.array([1e-3])
    )
    assert D[0] == 1e-6


def test_floor_rejects_negative_factor():
    with pytest.raises(ValueError, match="factor must not be negative"):
        ic.DiffusionFloor(factor=-1.0)


def test_volume_weight_integrates_to_the_chamber_volume(chamber):
    mesh = chamber.build_mesh(10)
    volume = assemble(volume_weight(mesh, "cylindrical") * ufl.dx(domain=mesh))
    assert np.isclose(volume, chamber.volume)


def test_volume_weight_rejects_unknown_coordinate_system(chamber):
    with pytest.raises(ValueError, match="coordinate_system must be"):
        volume_weight(chamber.build_mesh(4), "spherical")


def test_lumped_projection_reproduces_a_constant(chamber):
    mesh = chamber.build_mesh(8)
    V1 = fem.functionspace(mesh, ("Lagrange", 1))
    values = lumped_projection(mesh, V1, lambda r, z: np.full_like(r, 3.0))
    assert np.allclose(values, 3.0)


def test_lumped_projection_conserves_the_integral_of_a_thin_layer(chamber):
    """A layer much thinner than the cells keeps its volume integral."""
    mesh = chamber.build_mesh(8)
    V1 = fem.functionspace(mesh, ("Lagrange", 1))
    thickness = 5e-3

    f = fem.Function(V1)
    f.x.array[:] = lumped_projection(
        mesh, V1, lambda r, z: (r > chamber.can_radius - thickness).astype(float), k=40
    )
    integral = assemble(volume_weight(mesh, "cylindrical") * f * ufl.dx)
    R = chamber.can_radius
    exact = np.pi * (R**2 - (R - thickness) ** 2) * chamber.length
    assert np.isclose(integral, exact, rtol=0.05)


def test_potential_equals_the_bias_on_the_anode(long_chamber):
    mesh = long_chamber.build_mesh(20)
    potential = ic.LaplaceField().solve(mesh, long_chamber, bias=300.0)
    assert np.isclose(potential.x.array.max(), 300.0)


def test_potential_is_zero_on_the_cathode(long_chamber):
    mesh = long_chamber.build_mesh(20)
    potential = ic.LaplaceField().solve(mesh, long_chamber, bias=300.0)
    assert np.isclose(potential.x.array.min(), 0.0, atol=1e-9)


def test_field_at_mid_height_is_that_of_an_infinite_coax(long_chamber):
    """E(r) = V / (r ln(R / a)) far from the end caps."""
    mesh = long_chamber.build_mesh(40)
    potential = ic.LaplaceField().solve(mesh, long_chamber, bias=300.0)
    strength, _ = cell_field_length(mesh, potential)

    centroids = mesh.geometry.x[:, :2][mesh.geometry.dofmaps[0]].mean(axis=1)
    middle = np.abs(centroids[:, 1] - 0.2) < 0.02
    a, R = long_chamber.rod_radius, long_chamber.can_radius
    exact = 300.0 / (centroids[middle, 0] * np.log(R / a))
    assert np.allclose(strength[middle], exact, rtol=2e-3)


def test_cell_length_is_bounded_by_the_cell_size(long_chamber):
    mesh = long_chamber.build_mesh(10)
    potential = ic.LaplaceField().solve(mesh, long_chamber, bias=300.0)
    _, length = cell_field_length(mesh, potential)
    vertices = mesh.geometry.x[:, :2][mesh.geometry.dofmaps[0]]
    longest_edge = np.linalg.norm(vertices[:, 1] - vertices[:, 0], axis=1).max()
    assert length.max() <= 2 * longest_edge
