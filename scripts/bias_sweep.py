"""Saturation curve: chamber current vs bias for the reference 1 L chamber."""

import functools
import json
import sys
import time

import numpy as np

from ion_chambers import beta_source as bs

bs.efficiency_grid = functools.lru_cache(maxsize=4)(
    bs.efficiency_grid
)  # same geometry every time
from sweeps import A_ROD, RESULTS, VOL

from ion_chambers.chamber_model import Chamber

R = 0.05
L = VOL / (np.pi * (R**2 - A_ROD**2))
CI = 3.7e10
ACTS = [1 * CI, 10 * CI, 100 * CI]
volts = [float(v) for v in sys.argv[1].split(",")]
RESULTS.mkdir(exist_ok=True)
out = []
for V in volts:
    t0 = time.perf_counter()
    rec = {"V": V, "rows": {}}
    for n in (150, 300):
        ch = Chamber(A_ROD, R, L, V=V, N_r=n)
        for A in ACTS:
            o = ch.solve(A)
            rec["rows"].setdefault(f"{A:.3e}", {})[n] = o
        del ch
    for d in rec["rows"].values():
        d["I_extrap"] = 2 * d[300]["I"] - d[150]["I"]
    print(
        f"V={V:g}: "
        + " ".join(
            f"{k}:{d['I_extrap'] / d[300]['I_ideal']:.4f}(sc {d[300]['space_charge']:.3f})"
            for k, d in rec["rows"].items()
        )
        + f"  {time.perf_counter() - t0:.0f}s",
        flush=True,
    )
    out.append(rec)
    json.dump(out, open(RESULTS / f"bias_{sys.argv[2]}.json", "w"), default=float)
