import numpy as np
import pytest

from ion_chambers import beta


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
