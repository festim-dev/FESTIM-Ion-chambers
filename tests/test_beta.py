import numpy as np
import pytest

import ion_chambers as ic
from ion_chambers import beta
from ion_chambers.beta import efficiency_grid, sample_decays, wall_efficiency

AIR_DENSITY = 1.205  # kg/m3


@pytest.fixture
def huge_chamber() -> ic.CoaxialChamber:
    """So large that its end cap is a flat wall on the scale of the beta range."""
    return ic.CoaxialChamber(rod_radius=1e-3, can_radius=1.0, length=2.0)


def test_spectrum_is_normalised():
    """The beta spectrum integrates to one over the energy."""
    E, N = beta.beta_spectrum()
    assert np.isclose(np.trapezoid(N, E), 1.0)


def test_spectrum_is_positive():
    """The beta spectrum is never negative."""
    _, N = beta.beta_spectrum()
    assert np.all(N >= 0)


def test_mean_energy_matches_literature():
    """The mean beta energy is within 0.5% of the literature value of 5.69 keV."""
    assert np.isclose(beta.mean_energy(), 5.69e3, rtol=5e-3)


def test_ion_pairs_per_decay_in_air():
    """A decay that stops in air (W = 33.97 eV) makes about 168 ion pairs."""
    assert np.isclose(beta.ion_pairs_per_decay(33.97), 168, rtol=5e-3)


def test_ion_pairs_scale_inversely_with_w_value():
    """Halving the W value doubles the number of ion pairs per decay."""
    assert np.isclose(
        beta.ion_pairs_per_decay(20.0), 2 * beta.ion_pairs_per_decay(40.0)
    )


def test_range_increases_with_energy():
    """A more energetic beta has a longer range."""
    ranges = beta.mass_range(np.logspace(2, np.log10(beta.Q_BETA), 50))
    assert np.all(np.diff(ranges) > 0)


def test_maximum_range_in_air_is_a_few_millimetres():
    """A beta at the endpoint energy travels between 4 and 6 mm in air at 1 atm."""
    range_m = beta.mass_range(beta.Q_BETA) / 1.205
    assert 4e-3 < range_m < 6e-3


@pytest.mark.parametrize("energy", [50.0, 1e3, 5.7e3, 18e3])
def test_energy_from_mass_range_inverts_mass_range(energy):
    """Converting an energy to a range and back returns the energy within 0.1%."""
    assert np.isclose(
        beta.energy_from_mass_range(beta.mass_range(energy)), energy, rtol=1e-3
    )


def test_no_energy_left_at_zero_range():
    """A beta with no range left has no energy left."""
    assert beta.energy_from_mass_range(0.0) == 0.0


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
    assert np.isclose(energies.mean(), beta.mean_energy(), rtol=0.02)


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


def test_efficiency_grid_is_reused(chamber):
    """Asking for the same efficiency grid twice returns the stored result without
    recomputing."""
    first = efficiency_grid(chamber, AIR_DENSITY, nr=8, nz=8, n_rays=50)
    assert efficiency_grid(chamber, AIR_DENSITY, nr=8, nz=8, n_rays=50) is first


def test_efficiency_grid_is_shared_by_chambers_of_the_same_shape():
    """Two chamber objects with the same dimensions share one stored efficiency
    grid."""
    first = efficiency_grid(
        ic.CoaxialChamber(1e-3, 0.03, 0.1), AIR_DENSITY, nr=8, nz=8, n_rays=50
    )
    second = efficiency_grid(
        ic.CoaxialChamber(1e-3, 0.03, 0.1), AIR_DENSITY, nr=8, nz=8, n_rays=50
    )
    assert second is first


def test_efficiency_grid_shape(chamber):
    """The efficiency grid has one value per pair of r and z grid points."""
    _, _, f = efficiency_grid(chamber, AIR_DENSITY, nr=8, nz=9, n_rays=50)
    assert f.shape == (8, 9)
