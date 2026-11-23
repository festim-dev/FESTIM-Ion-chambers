"""Tritium betas: their spectrum, their range, and the energy they lose to walls.

Tritium betas travel up to about 5 mm in air at 1 atm, so a decay close to a wall
leaves only part of its energy in the gas. The share that stays in the gas, the
wall-loss efficiency f in [0, 1], is computed by a straight-track Monte Carlo:

* tracks are straight, with the Katz-Penfold practical range (scattering is folded
  into the practical range), and lose energy along the way following the range-energy
  relation
* emission is isotropic, and the energy is counted at the point of decay (the track
  is short compared with the drift gap)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ion_chambers.geometry import Geometry

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


def sample_decays(n_rays: int, seed: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Draws beta energies from the tritium spectrum and isotropic directions.

    Args:
        n_rays: the number of betas
        seed: the seed of the random generator

    Returns:
        the energies (eV), shape (n_rays,), and unit direction vectors, shape
        (n_rays, 3)
    """
    rng = np.random.default_rng(seed)

    # energies: invert the cumulative distribution of the spectrum
    E, N = beta_spectrum()
    cumulative = np.cumsum(N)
    cumulative /= cumulative[-1]
    energies = np.interp(rng.random(n_rays), cumulative, E)

    # directions: uniform on the sphere
    cos_polar = 2 * rng.random(n_rays) - 1
    azimuth = 2 * np.pi * rng.random(n_rays)
    sin_polar = np.sqrt(1 - cos_polar**2)
    directions = np.column_stack(
        [sin_polar * np.cos(azimuth), sin_polar * np.sin(azimuth), cos_polar]
    )
    return energies, directions


def wall_efficiency(
    points: np.ndarray,
    geometry: Geometry,
    density: float,
    n_rays: int = 3000,
    seed: int = 1,
) -> np.ndarray:
    """The fraction of the beta energy deposited in the gas, for decays at ``points``.

    Args:
        points: the decay positions in mesh coordinates, shape (n, 2)
        geometry: the chamber geometry
        density: the mass density of the gas (kg/m3)
        n_rays: the number of betas followed from each point
        seed: the seed of the random generator

    Returns:
        the efficiency at each point, between 0 and 1
    """
    energies, directions = sample_decays(n_rays, seed)
    ranges = mass_range(energies)

    efficiency = np.empty(len(points))
    chunk_size = 400  # bounds the (points x rays) arrays held in memory
    for start in range(0, len(points), chunk_size):
        chunk = slice(start, start + chunk_size)
        path_length = geometry.distance_to_wall(points[chunk], directions)
        # the range a beta still has when it reaches the wall, zero if it stopped
        remaining_range = np.clip(ranges[None, :] - path_length * density, 0, None)
        energy_at_wall = energy_from_mass_range(remaining_range)
        deposited_in_gas = energies[None, :] - energy_at_wall
        efficiency[chunk] = deposited_in_gas.mean(axis=1) / energies.mean()
    return efficiency


_EFFICIENCY_GRID_CACHE: dict = {}


def efficiency_grid(
    geometry: Geometry, density: float, nr: int = 110, nz: int = 110, n_rays: int = 3000
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The wall-loss efficiency tabulated on a grid refined towards the walls.

    The Monte Carlo is the slow part of setting up a problem and depends only on the
    shape of the chamber and the gas density, so results are kept and reused: a bias
    or activity sweep pays for it once.

    Args:
        geometry: the chamber geometry
        density: the mass density of the gas (kg/m3)
        nr: the number of grid points along r
        nz: the number of grid points along z
        n_rays: the number of betas followed from each point

    Returns:
        the r axis, the z axis and the efficiency, shape (nr, nz)
    """
    key = (type(geometry).__name__, geometry.dimensions, density, nr, nz, n_rays)
    if key in _EFFICIENCY_GRID_CACHE:
        return _EFFICIENCY_GRID_CACHE[key]

    r, z = geometry.wall_clustered_grid(nr, nz)
    rr, zz = np.meshgrid(r, z, indexing="ij")
    points = np.column_stack([rr.ravel(), zz.ravel()])
    efficiency = wall_efficiency(points, geometry, density, n_rays=n_rays)
    result = (r, z, efficiency.reshape(rr.shape))

    _EFFICIENCY_GRID_CACHE[key] = result
    return result
