import numpy as np
import pytest

from ion_chambers import beta


def test_spectrum_is_normalised():
    E, N = beta.beta_spectrum()
    assert np.isclose(np.trapezoid(N, E), 1.0)


def test_spectrum_is_positive():
    _, N = beta.beta_spectrum()
    assert np.all(N >= 0)


def test_mean_energy_matches_literature():
    """The mean tritium beta energy is 5.69 keV."""
    assert np.isclose(beta.mean_energy(), 5.69e3, rtol=5e-3)


def test_ion_pairs_per_decay_in_air():
    assert np.isclose(beta.ion_pairs_per_decay(33.97), 168, rtol=5e-3)


def test_ion_pairs_scale_inversely_with_w_value():
    assert np.isclose(
        beta.ion_pairs_per_decay(20.0), 2 * beta.ion_pairs_per_decay(40.0)
    )


def test_range_increases_with_energy():
    ranges = beta.mass_range(np.logspace(2, np.log10(beta.Q_BETA), 50))
    assert np.all(np.diff(ranges) > 0)


def test_maximum_range_in_air_is_a_few_millimetres():
    range_m = beta.mass_range(beta.Q_BETA) / 1.205
    assert 4e-3 < range_m < 6e-3


@pytest.mark.parametrize("energy", [50.0, 1e3, 5.7e3, 18e3])
def test_energy_from_mass_range_inverts_mass_range(energy):
    assert np.isclose(
        beta.energy_from_mass_range(beta.mass_range(energy)), energy, rtol=1e-3
    )


def test_no_energy_left_at_zero_range():
    assert beta.energy_from_mass_range(0.0) == 0.0
