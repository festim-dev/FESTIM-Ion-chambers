"""Ionisation source for a tritium-in-air coaxial ionisation chamber.

Tritium betas (Q = 18.59 keV, mean 5.69 keV) have ranges of up to ~5 mm in air,
so decays close to a wall deposit only part of their energy in the gas. This module
computes, by straight-line Monte Carlo, the local energy-deposition efficiency
f(r, z) in [0, 1] for a coaxial chamber (anode rod radius a, cathode radius R,
length L). The ion-pair source is then S(r, z) = A * (E_mean / W_air) * f(r, z).

Assumptions (prototype level):
  * straight tracks with the Katz-Penfold practical range (scattering is folded into
    the practical range); energy loss along the track from the range-energy relation
  * isotropic emission; the deposited energy is counted where the decay happens
    (the ~mm track length is small compared with the drift gap)
"""

import numpy as np

Q_BETA = 18.591e3  # eV, tritium endpoint
W_AIR = 33.97  # eV per ion pair in dry air
ME_C2 = 510.999e3  # eV
ALPHA_FS = 1 / 137.036
RHO_AIR = 1.205e-3  # g/cm3 at 20 C, 1 atm


def beta_spectrum(n=4000):
    """Tritium beta spectrum N(E) on a grid (allowed shape with Fermi function)."""
    E = np.linspace(1.0, Q_BETA - 1.0, n)
    W = E + ME_C2
    p = np.sqrt(W**2 - ME_C2**2)
    eta = 2 * ALPHA_FS * W / p  # Z = 2 daughter (He-3)
    fermi = 2 * np.pi * eta / (1 - np.exp(-2 * np.pi * eta))
    N = fermi * p * W * (Q_BETA - E) ** 2
    return E, N / np.trapezoid(N, E)


def range_air_m(E_eV):
    """Katz-Penfold practical range in air, metres (E in eV)."""
    E = np.asarray(E_eV) / 1e6
    R_mg_cm2 = 412.0 * E ** (1.265 - 0.0954 * np.log(E))
    return R_mg_cm2 * 1e-3 / RHO_AIR * 1e-2


E_GRID = np.logspace(0, np.log10(Q_BETA), 600)
R_GRID = range_air_m(E_GRID)


def energy_from_range(rng_m):
    return np.interp(rng_m, R_GRID, E_GRID, left=0.0)


def mean_energy():
    E, N = beta_spectrum()
    return np.trapezoid(E * N, E)


def ion_pairs_per_decay():
    return mean_energy() / W_AIR


def _sample(n_rays, seed=1):
    rng = np.random.default_rng(seed)
    E, N = beta_spectrum()
    cdf = np.cumsum(N)
    cdf /= cdf[-1]
    E0 = np.interp(rng.random(n_rays), cdf, E)
    mu = 2 * rng.random(n_rays) - 1
    ph = 2 * np.pi * rng.random(n_rays)
    s = np.sqrt(1 - mu**2)
    return E0, s * np.cos(ph), s * np.sin(ph), mu


def _dist_to_cyl(x0, ux, uy, rho, inside):
    """Distance along ray from (x0, 0) to cylinder of radius rho (2D projection)."""
    A = ux**2 + uy**2
    B = 2 * x0[:, None] * ux[None, :]
    C = (x0**2 - rho**2)[:, None]
    disc = B**2 - 4 * A[None, :] * C
    with np.errstate(invalid="ignore", divide="ignore"):
        sq = np.sqrt(np.where(disc > 0, disc, np.nan))
        if inside:  # point inside the cylinder: exit at the + root
            s = (-B + sq) / (2 * A[None, :])
        else:  # point outside a solid rod: hit at the - root if positive
            s = (-B - sq) / (2 * A[None, :])
            s = np.where(s > 0, s, np.inf)
    return np.where(np.isnan(s), np.inf, s)


def wall_efficiency(r, z, a, R, L, n_rays=3000, seed=1):
    """Energy-deposition efficiency f at points (r, z) (1D arrays of equal length)."""
    E0, ux, uy, uz = _sample(n_rays, seed)
    R0 = range_air_m(E0)
    out = np.empty(len(r))
    chunk = 400
    for i in range(0, len(r), chunk):
        ri, zi = r[i : i + chunk], z[i : i + chunk]
        s_out = _dist_to_cyl(ri, ux, uy, R, inside=True)
        s_rod = _dist_to_cyl(ri, ux, uy, a, inside=False)
        with np.errstate(divide="ignore"):
            s_z = np.where(
                uz[None, :] > 0,
                (L - zi[:, None]) / uz[None, :],
                np.where(uz[None, :] < 0, -zi[:, None] / uz[None, :], np.inf),
            )
        s = np.minimum(np.minimum(s_out, s_rod), s_z)
        dep = E0[None, :] - energy_from_range(np.clip(R0[None, :] - s, 0, None))
        out[i : i + chunk] = dep.mean(axis=1) / E0.mean()
    return out


def efficiency_grid(a, R, L, nr=110, nz=110, n_rays=3000):
    """f on a tensor grid clustered at every wall (rod, can, both caps)."""
    s = np.linspace(-1, 1, nr)
    rg = a + 0.5 * (R - a) * (1 + np.sin(0.5 * np.pi * s))
    t = np.linspace(-1, 1, nz)
    zg = 0.5 * L * (1 + np.sin(0.5 * np.pi * t))  # clustered at both caps
    RR, ZZ = np.meshgrid(rg, zg, indexing="ij")
    f = wall_efficiency(RR.ravel(), ZZ.ravel(), a, R, L, n_rays=n_rays)
    return rg, zg, f.reshape(RR.shape)


def volume_efficiency(a, R, L, n_pts=20000, n_rays=1500, seed=3):
    """Direct check: volume-averaged f from uniformly sampled decay points."""
    rng = np.random.default_rng(seed)
    r = np.sqrt(a**2 + (R**2 - a**2) * rng.random(n_pts))
    z = L * rng.random(n_pts)
    return wall_efficiency(r, z, a, R, L, n_rays=n_rays).mean()


if __name__ == "__main__":
    print(f"mean beta energy {mean_energy():.0f} eV")
    print(f"ion pairs / decay {ion_pairs_per_decay():.1f}")
    print(
        f"max range {range_air_m(Q_BETA) * 1e3:.2f} mm, "
        f"range at mean {range_air_m(mean_energy()) * 1e3:.2f} mm"
    )
    # planar-wall check: f vs distance from a flat wall (huge chamber)
    for d in [0, 0.2e-3, 0.5e-3, 1e-3, 2e-3, 5e-3]:
        f = wall_efficiency(np.array([0.5]), np.array([d]), 1e-3, 1.0, 2.0, 20000)
        print(f"d = {d * 1e3:4.1f} mm  f = {f[0]:.3f}")
