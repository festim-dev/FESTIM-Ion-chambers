# Parameter studies

## Mesh extrapolation

The current computed on one mesh carries an error proportional to the cell size (see
[Stabilisation](field.md#stabilisation)). Solving on two meshes and extrapolating to
zero cell size removes most of it. {func}`~ion_chambers.richardson` does the
arithmetic:

```python
import ion_chambers as ic

extrapolated = ic.richardson(coarse=current_150, fine=current_300)
```

## Response curves

{func}`~ion_chambers.response_curve` computes the current against activity on two
meshes and extrapolates each point. You give it two models of the same chamber that
differ only in their resolution:

```python
import numpy as np

import ion_chambers as ic


def build_model(resolution):
    return ic.IonChamberModel(
        geometry=chamber,  # your geometry
        gas=gas,  # your gas
        sources=ic.TritiumInGas(activity=0.0),
        bias=500.0,
        resolution=resolution,
    )


activities = np.logspace(8, 13, 6)  # Bq/m3
rows = ic.response_curve(
    coarse_model=build_model(resolution=150),
    fine_model=build_model(resolution=300),
    activities=activities,
)
ic.save_json(rows, "results/my_study.json")
```

Each row holds the fine-mesh reading plus `activity`, `current_coarse`, `current_fine`
and `current_extrapolated`.

To study another parameter, loop over it and build the two models again in each pass.
`scripts/bias_sweep.py` does this for the bias and `scripts/aspect_ratio_sweep.py` for
the chamber shape. The same pattern works for the pressure or the gas.

## Ready-made scripts

The `scripts/` folder holds the studies behind the [prototype results](../results.md).
Each writes a JSON file to `results/`.

```bash
python scripts/quick_check.py                  # 5 s end-to-end check
python scripts/response_curve.py               # current against activity
python scripts/bias_sweep.py 250 500 1000      # saturation curve
python scripts/aspect_ratio_sweep.py 0.4 1.27 7  # chamber shape at fixed volume
```

Run any of them with `--help` for the options. With the default resolutions of 150 and
300 cells, expect seconds to minutes per case.
