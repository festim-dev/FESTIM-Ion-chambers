import numpy as np
import pytest

import ion_chambers as ic


def test_volume(chamber):
    assert np.isclose(chamber.volume, 1e-3)


def test_from_volume_sets_the_length(chamber):
    expected = 1e-3 / (np.pi * (0.05**2 - 1.5e-3**2))
    assert np.isclose(chamber.length, expected)


@pytest.mark.parametrize("aspect_ratio", [0.4, 1.27, 7.0])
def test_from_aspect_ratio_keeps_the_volume(aspect_ratio):
    geometry = ic.CoaxialChamber.from_aspect_ratio(1e-3, 1.5e-3, aspect_ratio)
    assert np.isclose(geometry.volume, 1e-3)


@pytest.mark.parametrize("aspect_ratio", [0.4, 1.27, 7.0])
def test_from_aspect_ratio_sets_length_over_diameter(aspect_ratio):
    geometry = ic.CoaxialChamber.from_aspect_ratio(1e-3, 1.5e-3, aspect_ratio)
    assert np.isclose(geometry.length / (2 * geometry.can_radius), aspect_ratio)


@pytest.mark.parametrize(
    "attribute", ["rod_radius", "can_radius", "length", "feedthrough_gap"]
)
def test_rejects_non_positive_dimensions(attribute):
    dimensions = {"rod_radius": 1e-3, "can_radius": 0.05, "length": 0.1}
    dimensions[attribute] = -1.0
    with pytest.raises(ValueError, match=f"{attribute} must be positive"):
        ic.CoaxialChamber(**dimensions)


def test_rejects_wrong_type():
    with pytest.raises(TypeError, match="length must be a float, not str"):
        ic.CoaxialChamber(rod_radius=1e-3, can_radius=0.05, length="long")


def test_rejects_rod_wider_than_can():
    with pytest.raises(ValueError, match="must be smaller than can_radius"):
        ic.CoaxialChamber(rod_radius=0.06, can_radius=0.05, length=0.1)


def test_rejects_insulator_covering_the_caps():
    with pytest.raises(ValueError, match="leaves no cathode"):
        ic.CoaxialChamber(
            rod_radius=1e-3, can_radius=0.05, length=0.1, feedthrough_gap=0.05
        )


def test_geometries_with_the_same_dimensions_are_equal():
    """Equality and hashing ignore the mesh options: they key the wall-loss cache."""
    a = ic.CoaxialChamber(1e-3, 0.05, 0.1, axial_cells=10)
    assert a == ic.CoaxialChamber(1e-3, 0.05, 0.1, axial_cells=20)


@pytest.mark.parametrize(
    "point, expected",
    [
        ((1.5e-3, 0.05), True),  # on the rod
        ((0.05, 0.05), False),  # on the can
        ((0.02, 0.05), False),  # in the gas
    ],
)
def test_is_anode(chamber, point, expected):
    assert chamber.is_anode(np.array(point)) == expected


@pytest.mark.parametrize(
    "point, expected",
    [
        ((0.05, 0.05), True),  # on the can
        ((0.02, 0.0), True),  # on the bottom cap
        ((0.02, "top"), True),  # on the top cap
        ((0.003, 0.0), False),  # on the insulating ring
        ((1.5e-3, 0.05), False),  # on the rod
        ((0.02, 0.05), False),  # in the gas
    ],
)
def test_is_cathode(chamber, point, expected):
    r, z = point
    z = chamber.length if z == "top" else z
    assert chamber.is_cathode(np.array([r, z])) == expected


@pytest.mark.parametrize(
    "direction, expected",
    [
        ((1, 0, 0), 0.05 - 0.02),  # outwards to the can
        ((-1, 0, 0), 0.02 - 1.5e-3),  # inwards to the rod
        ((0, 0, 1), "to the top cap"),
        ((0, 0, -1), 0.03),  # down to the bottom cap
        ((0, 1, 0), np.sqrt(0.05**2 - 0.02**2)),  # tangentially, missing the rod
    ],
)
def test_distance_to_wall(chamber, direction, expected):
    if expected == "to the top cap":
        expected = chamber.length - 0.03
    distance = chamber.distance_to_wall(
        np.array([[0.02, 0.03]]), np.array([direction], dtype=float)
    )
    assert np.isclose(distance[0, 0], expected)


def test_mesh_spans_the_gas_radially(chamber):
    r = chamber.build_mesh(10).geometry.x[:, 0]
    assert np.allclose([r.min(), r.max()], [1.5e-3, 0.05])


def test_mesh_spans_the_gas_axially(chamber):
    z = chamber.build_mesh(10).geometry.x[:, 1]
    assert np.allclose([z.min(), z.max()], [0.0, chamber.length])


def test_axial_cells_override():
    geometry = ic.CoaxialChamber(1e-3, 0.05, 0.1, axial_cells=7)
    mesh = geometry.build_mesh(5)
    # a crossed rectangle mesh has four triangles per quadrilateral
    assert mesh.topology.index_map(2).size_local == 4 * 5 * 7


def test_wall_clustered_grid_spans_the_gas(chamber):
    r, z = chamber.wall_clustered_grid(11, 13)
    assert np.allclose([r[0], r[-1], z[0], z[-1]], [1.5e-3, 0.05, 0, chamber.length])
