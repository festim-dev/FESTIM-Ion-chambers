"""The tritium beta spectrum and the range of the betas in matter."""

import numpy as np

Q_BETA = 18.591e3
"""Endpoint energy of the tritium beta spectrum (eV)."""

_ELECTRON_REST_ENERGY = 510.999e3  # eV
_FINE_STRUCTURE = 1 / 137.036


def beta_spectrum(n: int = 4000) -> tuple[np.ndarray, np.ndarray]:
    """The tritium beta spectrum: allowed shape with the Fermi function.

    Args:
        n: the number of energy points

    Returns:
        the energies (eV) and the spectrum, normalised to unit integral
    """
    E = np.linspace(1.0, Q_BETA - 1.0, n)
    W = E + _ELECTRON_REST_ENERGY
    p = np.sqrt(W**2 - _ELECTRON_REST_ENERGY**2)
    eta = 2 * _FINE_STRUCTURE * W / p  # Z = 2 daughter (He-3)
    fermi = 2 * np.pi * eta / (1 - np.exp(-2 * np.pi * eta))
    N = fermi * p * W * (Q_BETA - E) ** 2
    return E, N / np.trapezoid(N, E)


def mean_energy() -> float:
    """The mean energy of a tritium beta (eV), about 5.7 keV."""
    E, N = beta_spectrum()
    return float(np.trapezoid(E * N, E))


def ion_pairs_per_decay(w_value: float) -> float:
    """The mean number of ion pairs made by one decay that stops in the gas.

    Args:
        w_value: the mean energy spent per ion pair in the gas (eV)
    """
    return mean_energy() / w_value


def mass_range(energy: np.ndarray | float) -> np.ndarray:
    """The Katz-Penfold practical range of an electron, as a mass thickness.

    Dividing by the mass density of the medium gives the range as a length.

    Args:
        energy: the electron energy (eV)

    Returns:
        the range (kg/m2)
    """
    E = np.asarray(energy) / 1e6
    range_mg_cm2 = 412.0 * E ** (1.265 - 0.0954 * np.log(E))
    return range_mg_cm2 * 1e-2


_ENERGY_GRID = np.logspace(0, np.log10(Q_BETA), 600)
_RANGE_GRID = mass_range(_ENERGY_GRID)


def energy_from_mass_range(remaining_range: np.ndarray) -> np.ndarray:
    """The energy of an electron that still has ``remaining_range`` to travel.

    The inverse of :func:`mass_range`.

    Args:
        remaining_range: the remaining range (kg/m2)

    Returns:
        the energy (eV), zero once the range is used up
    """
    return np.interp(remaining_range, _RANGE_GRID, _ENERGY_GRID, left=0.0)
