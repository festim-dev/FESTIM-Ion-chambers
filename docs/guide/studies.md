# Parameter studies

## Mesh extrapolation

The current computed on one mesh carries an error proportional to the cell size (see
[Stabilisation](transport.md#stabilisation)). Solving on two meshes and extrapolating
to zero cell size removes most of it. {func}`~ion_chambers.richardson` does the
arithmetic:

```python
import ion_chambers as ic

extrapolated = ic.richardson(coarse=current_150, fine=current_300)
```

## Response curves

{func}`~ion_chambers.response_curve` computes the current against activity on two
meshes and extrapolates each point. You give it two transport problems that describe
the same chamber on a coarse and on a fine mesh. Each needs a `ChamberCurrent` among
its exports.

`scripts/reference_case.py` has a function, `build_transport_problem`, that sets up
both problems for a given chamber, bias and resolution. With it:

```python
import numpy as np
from reference_case import build_transport_problem, reference_chamber

import ion_chambers as ic

chamber = reference_chamber()
activities = np.logspace(8, 13, 6)  # Bq/m3

rows = ic.response_curve(
    coarse_problem=build_transport_problem(chamber, bias=500.0, resolution=150),
    fine_problem=build_transport_problem(chamber, bias=500.0, resolution=300),
    activities=activities,
)
ic.save_json(rows, "results/my_study.json")
```

Each row holds `activity`, `current_coarse`, `current_fine`, `current_extrapolated`,
and the value of every other quantity exported by the fine problem, under its title.

To study another parameter, loop over it and build the two problems again in each
pass. `scripts/bias_sweep.py` does this for the bias and
`scripts/aspect_ratio_sweep.py` for the chamber shape. The same pattern works for the
pressure or the gas.

## Ready-made scripts

The `scripts/` folder holds the studies behind the [prototype results](../results.md).
Each writes a JSON file to `results/`.

```bash
python scripts/quick_check.py          # 5 s end-to-end check, every step spelled out
python scripts/response_curve.py       # current against activity
python scripts/bias_sweep.py           # saturation curve
python scripts/aspect_ratio_sweep.py   # chamber shape at fixed volume
```

Each script is a plain sequence of steps with its parameters at the top. Edit those to
change the case. With the default resolutions of 150 and 300 cells, expect seconds to
minutes per case.
