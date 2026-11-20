"""Wall losses: the share of the beta energy that ends up in the gas.

Tritium betas travel up to about 5 mm in air at 1 atm, so a decay close to a wall
leaves only part of its energy in the gas. The local efficiency f in [0, 1] is
computed by a straight-track Monte Carlo:

* tracks are straight, with the Katz-Penfold practical range (scattering is folded
  into the practical range), and lose energy along the way following the range-energy
  relation
* emission is isotropic, and the energy is counted at the point of decay (the track
  is short compared with the drift gap)
"""

from functools import lru_cache

import numpy as np

from ion_chambers.beta import beta_spectrum, energy_from_mass_range, mass_range
from ion_chambers.geometry import Geometry


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
    E, N = beta_spectrum()
    cdf = np.cumsum(N)
    cdf /= cdf[-1]
    energies = np.interp(rng.random(n_rays), cdf, E)
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
    chunk = 400  # bounds the (points x rays) arrays held in memory
    for start in range(0, len(points), chunk):
        block = slice(start, start + chunk)
        path = geometry.distance_to_wall(points[block], directions)
        remaining = np.clip(ranges[None, :] - path * density, 0, None)
        deposited = energies[None, :] - energy_from_mass_range(remaining)
        efficiency[block] = deposited.mean(axis=1) / energies.mean()
    return efficiency


@lru_cache(maxsize=16)
def efficiency_grid(
    geometry: Geometry, density: float, nr: int = 110, nz: int = 110, n_rays: int = 3000
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The wall-loss efficiency tabulated on a grid refined towards the walls.

    The Monte Carlo is the slow part of setting up a model and depends only on the
    geometry and the gas density, so results are cached: a bias or activity sweep pays
    for it once.

    Args:
        geometry: the chamber geometry
        density: the mass density of the gas (kg/m3)
        nr: the number of grid points along r
        nz: the number of grid points along z
        n_rays: the number of betas followed from each point

    Returns:
        the r axis, the z axis and the efficiency, shape (nr, nz)
    """
    r, z = geometry.wall_clustered_grid(nr, nz)
    rr, zz = np.meshgrid(r, z, indexing="ij")
    points = np.column_stack([rr.ravel(), zz.ravel()])
    f = wall_efficiency(points, geometry, density, n_rays=n_rays)
    return r, z, f.reshape(rr.shape)
