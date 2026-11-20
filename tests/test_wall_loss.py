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
    """Every sampled emission direction has unit length."""
    _, directions = sample_decays(500)
    assert np.allclose(np.linalg.norm(directions, axis=1), 1.0)


def test_directions_are_isotropic():
    """The sampled directions average to zero, as isotropic emission must."""
    _, directions = sample_decays(20000)
    assert np.allclose(directions.mean(axis=0), 0.0, atol=0.02)


def test_sampled_energies_follow_the_spectrum():
    """The mean of the sampled energies is within 2% of the mean of the spectrum."""
    energies, _ = sample_decays(20000)
    assert np.isclose(energies.mean(), ic.beta.mean_energy(), rtol=0.02)


def test_half_the_energy_is_lost_at_a_flat_wall(huge_chamber):
    """A decay on a flat wall deposits half of its energy in the gas, within 1%."""
    f = wall_efficiency(np.array([[0.5, 0.0]]), huge_chamber, AIR_DENSITY, 20000)
    assert np.isclose(f[0], 0.5, atol=0.01)


def test_nothing_is_lost_beyond_the_maximum_range(huge_chamber):
    """A decay 2 cm from any wall, beyond the reach of every beta, deposits all of
    its energy in the gas."""
    f = wall_efficiency(np.array([[0.5, 0.02]]), huge_chamber, AIR_DENSITY)
    assert np.isclose(f[0], 1.0)


def test_efficiency_increases_away_from_the_wall(huge_chamber):
    """The efficiency never decreases when moving away from a wall."""
    distances = np.linspace(0, 5e-3, 12)
    points = np.column_stack([np.full_like(distances, 0.5), distances])
    f = wall_efficiency(points, huge_chamber, AIR_DENSITY, 5000)
    assert np.all(np.diff(f) >= 0)


def test_thinner_gas_loses_more_to_the_wall(huge_chamber):
    """Halving the gas density lowers the efficiency 1 mm from a wall."""
    point = np.array([[0.5, 1e-3]])
    dense = wall_efficiency(point, huge_chamber, AIR_DENSITY)
    assert wall_efficiency(point, huge_chamber, AIR_DENSITY / 2) < dense


def test_efficiency_grid_is_cached(chamber):
    """Asking for the same efficiency grid twice returns the stored result without
    recomputing."""
    first = efficiency_grid(chamber, AIR_DENSITY, nr=8, nz=8, n_rays=50)
    assert efficiency_grid(chamber, AIR_DENSITY, nr=8, nz=8, n_rays=50) is first


def test_efficiency_grid_shape(chamber):
    """The efficiency grid has one value per pair of r and z grid points."""
    _, _, f = efficiency_grid(chamber, AIR_DENSITY, nr=8, nz=9, n_rays=50)
    assert f.shape == (8, 9)
