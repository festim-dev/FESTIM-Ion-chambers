import numpy as np
import pytest

import ion_chambers as ic
from ion_chambers.wall_loss import efficiency_grid, sample_decays, wall_efficiency

AIR_DENSITY = 1.205  # kg/m3


@pytest.fixture
def huge_chamber() -> ic.CoaxialChamber:
    """So large that its end cap is a flat wall on the scale of the beta range."""
    return ic.CoaxialChamber(rod_radius=1e-3, can_radius=1.0, length=2.0)


def test_directions_are_unit_vectors():
    _, directions = sample_decays(500)
    assert np.allclose(np.linalg.norm(directions, axis=1), 1.0)


def test_directions_are_isotropic():
    _, directions = sample_decays(20000)
    assert np.allclose(directions.mean(axis=0), 0.0, atol=0.02)


def test_sampled_energies_follow_the_spectrum():
    energies, _ = sample_decays(20000)
    assert np.isclose(energies.mean(), ic.beta.mean_energy(), rtol=0.02)


def test_half_the_energy_is_lost_at_a_flat_wall(huge_chamber):
    """A decay on a flat wall sends half of its betas straight into it."""
    f = wall_efficiency(np.array([[0.5, 0.0]]), huge_chamber, AIR_DENSITY, 20000)
    assert np.isclose(f[0], 0.5, atol=0.01)


def test_nothing_is_lost_beyond_the_maximum_range(huge_chamber):
    f = wall_efficiency(np.array([[0.5, 0.02]]), huge_chamber, AIR_DENSITY)
    assert np.isclose(f[0], 1.0)


def test_efficiency_increases_away_from_the_wall(huge_chamber):
    distances = np.linspace(0, 5e-3, 12)
    points = np.column_stack([np.full_like(distances, 0.5), distances])
    f = wall_efficiency(points, huge_chamber, AIR_DENSITY, 5000)
    assert np.all(np.diff(f) >= 0)


def test_thinner_gas_loses_more_to_the_wall(huge_chamber):
    point = np.array([[0.5, 1e-3]])
    dense = wall_efficiency(point, huge_chamber, AIR_DENSITY)
    assert wall_efficiency(point, huge_chamber, AIR_DENSITY / 2) < dense


def test_efficiency_grid_is_cached(chamber):
    first = efficiency_grid(chamber, AIR_DENSITY, nr=8, nz=8, n_rays=50)
    assert efficiency_grid(chamber, AIR_DENSITY, nr=8, nz=8, n_rays=50) is first


def test_efficiency_grid_shape(chamber):
    _, _, f = efficiency_grid(chamber, AIR_DENSITY, nr=8, nz=9, n_rays=50)
    assert f.shape == (8, 9)
