"""Coaxial tritium ionisation chamber in room air, FESTIM 2.2 (drift terms).

Geometry (2D axisymmetric, r horizontal, z vertical):
    anode rod   r = a           at +V
    cathode     r = R, and the end caps z = 0, L for r > a + g     at 0 V
    insulators  end caps for a < r < a + g (rod feedthroughs), no flux

Species (number densities, m^-3):
    e   free electrons      z = -1, mu_e
    P   positive ions       z = +1, mu_p
    N   negative ions (O2-) z = -1, mu_n

    S(r, z) = A * (E_mean / W_air) * f_wall(r, z)        ion pairs / m3 / s
    e  -> N            attachment to O2, nu_att
    P + N -> 0         ion-ion recombination, alpha
    P + e -> 0         electron-ion recombination, beta

Transport: J = -D grad(c) + c v with v = -z mu grad(phi), through
festim.ElectromigrationTerm. The field is the fixed Laplace field (no space charge).

Numerics: FESTIM does not stabilise drift-dominated transport (ion cell Peclet numbers
reach ~100), so each species gets a cell-local diffusion floor
D = max(D_phys, mu |E| h_E / 2), stored as a DG0 field (CellwiseDiffusionMaterial),
with h_E the cell length along the field. The floor vanishes where the field is weak,
so it cannot fake losses in low-field corners. The drift velocity is kept at exactly
-z mu grad(phi) by MobilityElectromigration, which uses the mobility instead of the
Einstein relation with the (floored) D. The remaining error is first order in h; the
sweep scripts Richardson-extrapolate the current from N_r = 150 and 300.

Quick check:  python -m ion_chambers.chamber_model
"""

import time

import festim as F
import numpy as np
import ufl
from dolfinx import fem
from dolfinx import mesh as dmesh
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
from scipy.interpolate import RegularGridInterpolator

from . import beta_source as bs

E_CHARGE = 1.602176634e-19
EPS0 = 8.8541878128e-12
T_GAS = 293.15  # K
KT_E = F.k_B * T_GAS  # thermal voltage, V
BQ_PER_CI = 3.7e10

PHYS = {
    "mu_p": 2.0e-4,  # m2/V/s, positive cluster ions in air
    "mu_n": 2.2e-4,  # m2/V/s, O2- (hydrated) in air
    "mu_e": 0.04,  # m2/V/s, electrons in air, low E/N
    "eps_k_e": 1.0,  # V, electron characteristic energy (D_e / mu_e)
    "nu_att": 1.0e8,  # 1/s, three-body attachment to O2 at 1 atm
    "alpha": 2.0e-12,  # m3/s, ion-ion recombination in air
    "beta": 2.0e-13,  # m3/s, electron-ion recombination
}


class MobilityElectromigration(F.ElectromigrationTerm):
    """Electromigration with the drift velocity set by a mobility, v = -z mu grad(phi).

    FESTIM's ElectromigrationTerm ties the mobility to D through Einstein's relation.
    Here D carries a cell-local numerical floor, so the drift is decoupled from it.
    """

    def __init__(self, species, charge, mobility, potential, subdomain):
        super().__init__(
            species=species, charge=charge, potential=potential, subdomain=subdomain
        )
        self.mobility = mobility

    def drift_velocity(self, D, temperature):
        return -self.charge * self.mobility * ufl.grad(self.potential.fenics_object)


class CellwiseDiffusionMaterial(F.Material):
    """Material whose per-species diffusivity is a DG0 field."""

    def __init__(self, D_fields, **kw):
        super().__init__(
            D_0={k: 1.0 for k in D_fields}, E_D={k: 0.0 for k in D_fields}, **kw
        )
        self.D_fields = D_fields

    def get_diffusion_coefficient(self, mesh=None, temperature=None, species=None):
        return self.D_fields[species.name]


def cell_velocity_length(m, phi):
    """|E| at cell centroids and the cell length along E (SUPG-style h)."""
    P = m.geometry.x[:, :2][m.geometry.dofmaps[0]]  # (ncells, 3, 2)
    cen = P.mean(axis=1)
    pts = np.c_[cen, np.zeros(len(cen))]
    cells = np.arange(len(cen), dtype=np.int32)
    Vv = fem.functionspace(m, ("DG", 1, (2,)))
    E = fem.Function(Vv)
    E.interpolate(fem.Expression(-ufl.grad(phi), Vv.element.interpolation_points))
    Ec = E.eval(pts, cells)
    Emag = np.linalg.norm(Ec, axis=1)
    vhat = Ec / np.maximum(Emag[:, None], 1e-300)
    edges = np.stack([P[:, 1] - P[:, 0], P[:, 2] - P[:, 1], P[:, 0] - P[:, 2]], axis=1)
    h_v = np.max(np.abs(np.einsum("cej,cj->ce", edges, vhat)), axis=1)
    return Emag, h_v


def lumped_projection(m, V1, func, k=6):
    """Nodal values f_i = int f phi_i r dx / int phi_i r dx, sampled on k^2 points/cell.

    Conserves int f r dx exactly, so a wall layer thinner than the cells is not smeared
    into extra loss the way nodal interpolation would.
    """
    P = m.geometry.x[:, :2][m.geometry.dofmaps[0]]  # (nc, 3, 2)
    dofs = V1.dofmap.list  # same local vertex order as the P1 geometry
    lam = []
    for i in range(k):
        for j in range(k - i):
            lam.append(((i + 1 / 3) / k, (j + 1 / 3) / k))  # upward sub-triangles
            if i + j < k - 1:
                lam.append(((i + 2 / 3) / k, (j + 2 / 3) / k))  # downward ones
    lam = np.array(lam)  # (nq, 2): barycentric (l1, l2), l0 = 1 - l1 - l2
    B = np.c_[1 - lam.sum(1), lam]  # (nq, 3) = phi_i at sample points
    X = np.einsum("qv,cvd->cqd", B, P)  # (nc, nq, 2)
    e1, e2 = P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]
    area = 0.5 * np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0])
    f = func(X[..., 0].ravel(), X[..., 1].ravel()).reshape(X.shape[:2])
    w = (area[:, None] / len(lam)) * X[..., 0]  # r-weighted
    num = np.zeros(V1.dofmap.index_map.size_local + V1.dofmap.index_map.num_ghosts)
    den = np.zeros_like(num)
    np.add.at(num, dofs, np.einsum("cq,qv->cv", w * f, B))
    np.add.at(den, dofs, np.einsum("cq,qv->cv", w, B))
    return num / den


def build_mesh(a, R, L, N_r, N_z, gamma=1.6):
    m = dmesh.create_rectangle(
        MPI.COMM_WORLD,
        [[0, -1], [1, 1]],
        [N_r, N_z],
        dmesh.CellType.triangle,
        diagonal=dmesh.DiagonalType.crossed,
    )
    x = m.geometry.x
    s, t = x[:, 0].copy(), x[:, 1].copy()
    x[:, 0] = a * (R / a) ** s  # h_r proportional to r: |v| h_r uniform
    x[:, 1] = 0.5 * L * (1 + np.tanh(gamma * t) / np.tanh(gamma))  # refine caps
    return m


def solve_potential(m, a, R, L, g, V):
    """Laplace field in (r, z), P2, axisymmetric."""
    Vh = fem.functionspace(m, ("Lagrange", 2))
    tol = 1e-9 * R
    fdim = m.topology.dim - 1
    anode = dmesh.locate_entities_boundary(
        m, fdim, lambda x: np.isclose(x[0], a, atol=tol)
    )
    cath = dmesh.locate_entities_boundary(
        m, fdim, lambda x: _is_cathode(x, a, R, L, g, tol)
    )
    bcs = [
        fem.dirichletbc(
            fem.Constant(m, V), fem.locate_dofs_topological(Vh, fdim, anode), Vh
        ),
        fem.dirichletbc(
            fem.Constant(m, 0.0), fem.locate_dofs_topological(Vh, fdim, cath), Vh
        ),
    ]
    u, w = ufl.TrialFunction(Vh), ufl.TestFunction(Vh)
    r = ufl.SpatialCoordinate(m)[0]
    a_form = r * ufl.dot(ufl.grad(u), ufl.grad(w)) * ufl.dx
    L_form = fem.Constant(m, 0.0) * w * ufl.dx
    prob = LinearProblem(
        a_form,
        L_form,
        bcs=bcs,
        petsc_options_prefix="phi_",
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    )
    phi = prob.solve()
    return phi


def space_charge_solver(m, a, R, L, g):
    """Perturbation potential phi1 from a charge density: -div(eps0 r grad phi1) = r rho."""
    Vh = fem.functionspace(m, ("Lagrange", 2))
    V1 = fem.functionspace(m, ("Lagrange", 1))
    tol = 1e-9 * R
    fdim = m.topology.dim - 1
    facets = dmesh.locate_entities_boundary(
        m,
        fdim,
        lambda x: np.isclose(x[0], a, atol=tol) | _is_cathode(x, a, R, L, g, tol),
    )
    bc = fem.dirichletbc(
        fem.Constant(m, 0.0), fem.locate_dofs_topological(Vh, fdim, facets), Vh
    )
    rho = fem.Function(V1)
    u, w = ufl.TrialFunction(Vh), ufl.TestFunction(Vh)
    r = ufl.SpatialCoordinate(m)[0]
    prob = LinearProblem(
        r * ufl.dot(ufl.grad(u), ufl.grad(w)) * ufl.dx,
        r * rho / EPS0 * w * ufl.dx,
        bcs=[bc],
        petsc_options_prefix="sc_",
        petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
    )
    return prob, rho


def _is_cathode(x, a, R, L, g, tol):
    cap = np.isclose(x[1], 0, atol=tol) | np.isclose(x[1], L, atol=tol)
    return np.isclose(x[0], R, atol=tol) | (cap & (x[0] >= a + g - tol))


def chamber_volume(a, R, L):
    return np.pi * (R**2 - a**2) * L


class Chamber:
    def __init__(
        self,
        a,
        R,
        L,
        V=500.0,
        g=3e-3,
        N_r=120,
        N_z=None,
        phys=None,
        n_rays=3000,
        floor_factor=1.0,
        nz_factor=0.35,
    ):
        self.floor_factor = floor_factor
        self.a, self.R, self.L, self.V, self.g = a, R, L, V, g
        self.phys = dict(PHYS, **(phys or {}))
        p = self.phys
        if N_z is None:
            # centre spacing ~ outer radial spacing, capped
            h_out = R * np.log(R / a) / N_r
            N_z = int(np.clip(L / h_out * nz_factor, 0.6 * N_r, 1.0 * N_r))
        self.N_r, self.N_z = N_r, N_z
        self.vol = chamber_volume(a, R, L)
        self.pairs_per_decay = bs.ion_pairs_per_decay()

        t0 = time.perf_counter()
        m = build_mesh(a, R, L, N_r, N_z)
        self.dmesh = m
        self.phi = solve_potential(m, a, R, L, g, V)
        self._sc_prob, self._sc_rho = space_charge_solver(m, a, R, L, g)

        # ---- cell-local diffusion floors: D = max(D_phys, mu |E| h_v / 2) -----
        D_phys = {
            "P": p["mu_p"] * KT_E,
            "N": p["mu_n"] * KT_E,
            "e": p["mu_e"] * p["eps_k_e"],
        }
        mu = {"P": p["mu_p"], "N": p["mu_n"], "e": p["mu_e"]}
        z = {"P": 1, "N": -1, "e": -1}
        self.mu, self.z, self.D_phys = mu, z, D_phys
        Emag, h_v = cell_velocity_length(m, self.phi)
        V0 = fem.functionspace(m, ("DG", 0))
        self.D, self.D_stats = {}, {}
        for k, value in D_phys.items():
            Dk = fem.Function(V0, name=f"D_{k}")
            vals = np.maximum(value, self.floor_factor * 0.5 * mu[k] * Emag * h_v)
            Dk.x.array[:] = vals
            self.D[k] = Dk
            self.D_stats[k] = (vals.min(), np.median(vals), vals.max())

        # ---- wall-loss efficiency field f(r, z) -----------------------------
        rg, zg, fg = bs.efficiency_grid(a, R, L, n_rays=n_rays)
        interp = RegularGridInterpolator(
            (rg, zg), fg, bounds_error=False, fill_value=None
        )
        V1 = fem.functionspace(m, ("Lagrange", 1))
        self.f_wall = fem.Function(V1)
        self.f_wall.x.array[:] = lumped_projection(
            m,
            V1,
            lambda r, z: np.clip(
                interp(np.c_[np.clip(r, a, R), np.clip(z, 0, L)]), 0, 1
            ),
        )
        self.f_grid = (rg, zg, fg)
        self.S0 = fem.Constant(m, 0.0)  # pairs/m3/s for f = 1

        self._build_problem()
        self.t_setup = time.perf_counter() - t0

    # -------------------------------------------------------------------------
    def _build_problem(self):
        a, R, L, g, p = self.a, self.R, self.L, self.g, self.phys
        model = F.HydrogenTransportProblem()
        model.mesh = F.Mesh(self.dmesh, coordinate_system="cylindrical")

        self.e = F.Species("e")
        self.P = F.Species("P")
        self.N = F.Species("N")
        sp = {"e": self.e, "P": self.P, "N": self.N}
        model.species = list(sp.values())

        air = CellwiseDiffusionMaterial({k: self.D[k] for k in sp}, name="air")
        tol = 1e-9 * R
        gas = F.VolumeSubdomain(id=1, material=air, locator=lambda x: x[0] >= a - tol)
        anode = F.SurfaceSubdomain(
            id=1, locator=lambda x: np.isclose(x[0], a, atol=tol)
        )
        cathode = F.SurfaceSubdomain(
            id=2, locator=lambda x: _is_cathode(x, a, R, L, g, tol)
        )
        model.subdomains = [gas, anode, cathode]
        self.surfaces = {"anode": anode, "cathode": cathode}

        model.temperature = T_GAS

        source = self.S0 * self.f_wall
        model.sources = [
            F.ParticleSource(value=source, volume=gas, species=self.e),
            F.ParticleSource(value=source, volume=gas, species=self.P),
        ]
        model.reactions = [
            F.ArrheniusReaction(
                reactant=self.e,
                product=self.N,
                k_0=p["nu_att"],
                E_k=0.0,
                p_0=0.0,
                E_p=0.0,
                volume=gas,
            ),
            F.ArrheniusReaction(
                reactant=[self.P, self.N],
                product=None,
                k_0=p["alpha"],
                E_k=0.0,
                volume=gas,
            ),
            F.ArrheniusReaction(
                reactant=[self.P, self.e],
                product=None,
                k_0=p["beta"],
                E_k=0.0,
                volume=gas,
            ),
        ]
        model.drift_terms = [
            MobilityElectromigration(
                species=sp[k],
                charge=self.z[k],
                mobility=self.mu[k],
                potential=self.phi,
                subdomain=gas,
            )
            for k in sp
        ]
        # electrodes absorb what the field drives into them (outflow) and inject
        # nothing (c = 0 for the species the field pushes away)
        model.boundary_conditions = [
            F.OutflowBC(subdomain=anode, species=self.e),
            F.OutflowBC(subdomain=anode, species=self.N),
            F.FixedConcentrationBC(subdomain=anode, value=0.0, species=self.P),
            F.OutflowBC(subdomain=cathode, species=self.P),
            F.FixedConcentrationBC(subdomain=cathode, value=0.0, species=self.e),
            F.FixedConcentrationBC(subdomain=cathode, value=0.0, species=self.N),
        ]
        model.settings = F.Settings(
            atol=1e-30, rtol=1e-10, max_iterations=40, transient=False
        )
        model.initialise()
        self.model = model

    # -------------------------------------------------------------------------
    def solve(self, activity_Bq_m3):
        """Steady state for a given tritium activity concentration (Bq/m3)."""
        self.S0.value = activity_Bq_m3 * self.pairs_per_decay
        self.model.run()
        return self.readout(activity_Bq_m3)

    def fields(self):
        return {k: getattr(self, k).post_processing_solution for k in ("e", "P", "N")}

    def readout(self, activity):
        m = self.dmesh
        r = ufl.SpatialCoordinate(m)[0]
        n = ufl.FacetNormal(m)
        c = self.fields()
        gphi = ufl.grad(self.phi)
        ds = self.model.ds
        dx = ufl.dx(domain=m)
        q = {"P": 1, "N": -1, "e": -1}

        def flux(k, surf_id):
            v = -self.z[k] * self.mu[k] * gphi
            J = -self.D[k] * ufl.grad(c[k]) + c[k] * v
            return 2 * np.pi * _asm(r * ufl.dot(J, n) * ds(surf_id))

        G = 2 * np.pi * _asm(r * self.S0 * self.f_wall * dx)
        Rec = (
            2
            * np.pi
            * _asm(
                r
                * (
                    self.phys["alpha"] * c["P"] * c["N"]
                    + self.phys["beta"] * c["P"] * c["e"]
                )
                * dx
            )
        )
        I_an = -E_CHARGE * sum(q[k] * flux(k, 1) for k in q)  # negative charge in
        I_ca = E_CHARGE * sum(q[k] * flux(k, 2) for k in q)
        I_vol = E_CHARGE * (G - Rec)
        # Shockley-Ramo: I = e * sum_s q_s int J_s . (-grad phi_w), phi_w = phi / V.
        # Exact in steady state (div of the net charge flux vanishes) and free of
        # boundary-gradient errors: this is the reading.
        Jq = sum(
            q[k] * (-self.D[k] * ufl.grad(c[k]) - self.z[k] * self.mu[k] * c[k] * gphi)
            for k in q
        )
        I_ramo = E_CHARGE * 2 * np.pi * _asm(r * ufl.dot(Jq, -gphi / self.V) * dx)
        I_ideal = E_CHARGE * activity * self.pairs_per_decay * self.vol
        rho = 2 * np.pi * _asm(r * abs(c["P"] - c["N"] - c["e"]) * dx) / self.vol
        # field distortion from the space charge (Poisson solve on the converged
        # densities): ||grad phi1|| / ||grad phi0|| in the r-weighted L2 norm
        self._sc_rho.interpolate(
            fem.Expression(
                E_CHARGE * (c["P"] - c["N"] - c["e"]),
                self._sc_rho.function_space.element.interpolation_points,
            )
        )
        phi1 = self._sc_prob.solve()
        e1 = _asm(r * ufl.dot(ufl.grad(phi1), ufl.grad(phi1)) * dx)
        e0 = _asm(r * ufl.dot(gphi, gphi) * dx)
        sc_field = np.sqrt(e1 / e0)
        sc_pot = np.max(np.abs(phi1.x.array)) / self.V
        rho_max = np.max(np.abs(c["P"].x.array - c["N"].x.array - c["e"].x.array))
        return {
            "activity": activity,
            "I": I_ramo,
            "I_anode": I_an,
            "I_cathode": I_ca,
            "I_balance": I_vol,
            "I_ideal": I_ideal,
            "generation": G,
            "recombination": Rec,
            "wall_eff": G * E_CHARGE / I_ideal if I_ideal > 0 else np.nan,
            "collection_eff": (G - Rec) / G if G > 0 else np.nan,
            # space-charge field relative to the applied one
            "space_charge": sc_field,
            "space_charge_pot": sc_pot,
            "space_charge_crude": E_CHARGE
            * rho_max
            * (self.R - self.a) ** 2
            / (EPS0 * self.V),
            "n_mean": rho,
        }


def _asm(form):
    return MPI.COMM_WORLD.allreduce(fem.assemble_scalar(fem.form(form)), op=MPI.SUM)


if __name__ == "__main__":
    # reference 1 L chamber on a coarse mesh: a quick end-to-end check
    a, Vol = 1.5e-3, 1e-3
    R = 0.05
    L = Vol / (np.pi * (R**2 - a**2))
    ch = Chamber(a, R, L, V=500.0, N_r=100)
    print(f"setup {ch.t_setup:.1f}s  N_r={ch.N_r}  N_z={ch.N_z}")
    for k, (lo, med, hi) in ch.D_stats.items():
        print(f"  D_{k}: min {lo:.3g}  median {med:.3g}  max {hi:.3g} m2/s")
    for A in [3.7e4, 3.7e10, 3.7e13]:
        t0 = time.perf_counter()
        out = ch.solve(A)
        print(
            f"{time.perf_counter() - t0:.1f}s", {k: f"{v:.4g}" for k, v in out.items()}
        )
