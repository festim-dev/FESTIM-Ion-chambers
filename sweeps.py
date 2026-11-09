"""Response curve (reference 1 L chamber) and L/D sweep at fixed 1 L volume."""
import sys, gc, json, time
import numpy as np
from chamber_model import Chamber, chamber_volume

A_ROD, VOL, VBIAS = 1.5e-3, 1e-3, 500.0
NR = (150, 300)


def R_for_aspect(ar):
    # pi (R^2 - a^2) L = VOL with L = 2 R ar
    from scipy.optimize import brentq
    return brentq(lambda R: chamber_volume(A_ROD, R, 2 * R * ar) - VOL, 5e-3, 1.0)


def run(R, L, activities, tag):
    rows = {n: [] for n in NR}
    for n in NR:
        ch = Chamber(A_ROD, R, L, V=VBIAS, N_r=n, N_z=None)
        for A in activities:
            o = ch.solve(A); rows[n].append(o)
            print(f"[{tag}] N_r={n} A={A:.2e} I/Iid={o['I']/o['I_ideal']:.4f} "
                  f"coll={o['collection_eff']:.4f} sc={o['space_charge']:.3g}", flush=True)
        wall = ch.f_grid
        del ch; gc.collect()
    out = []
    for c, f in zip(rows[NR[0]], rows[NR[1]]):
        r = dict(f)
        r["I_extrap"] = 2 * f["I"] - c["I"]  # first-order error in h
        r["I_coarse"], r["I_fine"] = c["I"], f["I"]
        out.append(r)
    return out


if __name__ == "__main__":
    which = sys.argv[1]
    if which == "response":
        R = 0.05
        L = VOL / (np.pi * (R**2 - A_ROD**2))
        acts = np.r_[np.logspace(3, 8, 3), np.logspace(9, 14, 11)]
        res = dict(R=R, L=L, a=A_ROD, V=VBIAS, rows=run(R, L, acts, "ref"))
        json.dump(res, open("response.json", "w"), indent=1)
    else:
        ars = [float(x) for x in sys.argv[2].split(",")]
        acts = np.r_[3.7e7, np.logspace(8, 14, 13)]
        allres = []
        for ar in ars:
            R = R_for_aspect(ar); L = 2 * R * ar
            t0 = time.perf_counter()
            allres.append(dict(aspect=ar, R=R, L=L, a=A_ROD, V=VBIAS,
                               rows=run(R, L, acts, f"L/D={ar}")))
            print(f"== L/D={ar} R={R*100:.2f} cm L={L*100:.2f} cm done in "
                  f"{time.perf_counter()-t0:.0f}s", flush=True)
            json.dump(allres, open(f"aspect_{sys.argv[3]}.json", "w"), indent=1)
