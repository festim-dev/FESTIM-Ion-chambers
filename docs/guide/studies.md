# Parameter studies

The package has no function for parameter studies. A study is a loop in your own
script, which keeps every step visible and easy to change. This page shows the two
patterns you need, and the `examples/` folder holds complete scripts to copy.

## Mesh extrapolation

The current computed on one mesh carries an error proportional to the cell size (see
[Stabilisation](transport.md#stabilisation)). Solving on two meshes and extrapolating
to zero cell size removes most of it. This is Richardson extrapolation of a
first-order error:

```python
difference = fine_current.value - coarse_current.value
extrapolated_current = fine_current.value + difference / (refinement - 1)
```

Here `refinement` is how many times more cells the fine mesh has in both directions.
With twice as many, this is `2 * fine - coarse`. The difference between the fine and the extrapolated value is a
measure of the error that remains.

## Changing the activity between runs

A problem is initialised once. After that, the activity of a
{class}`~ion_chambers.TritiumInGas` source can be changed and the problem run again,
which is much cheaper than building a new problem:

```python
transport.initialise()

for activity in activities:
    tritium.activity = activity
    transport.run()
    print(activity, current.value)
```

The activity is the only thing that can change this way. Another bias, pressure, gas
or shape needs new problems, so loop over that parameter and build the problems inside
the loop.

## A complete response curve

`examples/response_curve.py` puts the two patterns together: the same chamber on a
coarse and on a fine mesh, solved at every activity and extrapolated. It uses
`build_transport_problem` from the
[example case](../example_case.md#the-script) to set up each problem.

```{literalinclude} ../../examples/response_curve.py
:language: python
```

It prints the extrapolated current at each activity:

```text
activity 1.00e+03 Bq/m3   current 2.6592e-17 A
activity 3.16e+05 Bq/m3   current 8.4092e-15 A
activity 1.00e+08 Bq/m3   current 2.6592e-12 A
activity 1.00e+09 Bq/m3   current 2.6587e-11 A
activity 3.16e+09 Bq/m3   current 8.4040e-11 A
activity 1.00e+10 Bq/m3   current 2.6540e-10 A
activity 3.16e+10 Bq/m3   current 8.3576e-10 A
activity 1.00e+11 Bq/m3   current 2.6099e-09 A
activity 3.16e+11 Bq/m3   current 7.9693e-09 A
activity 1.00e+12 Bq/m3   current 2.3125e-08 A
activity 3.16e+12 Bq/m3   current 6.1242e-08 A
activity 1.00e+13 Bq/m3   current 1.4420e-07 A
activity 3.16e+13 Bq/m3   current 3.0617e-07 A
activity 1.00e+14 Bq/m3   current 6.0818e-07 A
saved results/response.json
```

## Ready-made examples

The `examples/` folder holds the studies of the
[example case](../example_case.md#studies-of-this-case). They come in pairs: a sweep that solves the model and writes a JSON file to `results/`,
and a post-processing script that reads that file, prints the numbers of interest and
saves a figure.

| Sweep | Post-processing | What it gives |
|---|---|---|
| `response_curve.py`: current against activity | `calibration.py` | the sensitivity, the end of the linear range, and the activity that corresponds to a measured current |
| `bias_sweep.py`: current against bias | `saturation_curve.py` | the bias for 99% collection, and the saturation current from a fit of 1/I against 1/V² |
| `aspect_ratio_sweep.py`: chamber shape at fixed volume | `shape_comparison.py` | the sensitivity and the linear range of each shape |

```bash
python examples/response_curve.py
python examples/calibration.py
```

Each script is a plain sequence of steps with its parameters at the top. Edit those to
change the case.

### Bias

`saturation_curve.py` prints, for each activity of `bias_sweep.py`, the bias that
collects 99% of the charge and the saturation current recovered by the fit:

```text
activity 1 Ci/m3
  99% collection at 432 V
  saturation current, fitted  9.8426e-10 A
  saturation current, true    9.8353e-10 A
  error of the fit            +0.07%
activity 10 Ci/m3
  99% collection at 1346 V
  saturation current, fitted  9.8290e-09 A
  saturation current, true    9.8353e-09 A
  error of the fit            -0.06%
activity 100 Ci/m3
  99% collection needs more than 2000 V
  saturation current, fitted  9.5463e-08 A
  saturation current, true    9.8353e-08 A
  error of the fit            -2.94%
saved results/saturation_curve.png
```

```{figure} ../images/guide_saturation_curve.png
:alt: Two plots. On the left the collection efficiency rises with the bias, more slowly at higher activity. On the right the inverse of the current against the inverse square of the bias, with a straight line fitted to the points at high bias.

The figure it saves. Left: the collection efficiency against the bias, with the 99%
target as a dashed line. Right: the saturation current over the current against
1/bias², with the line fitted from 500 V up.
```

### Chamber shape

`shape_comparison.py` prints, for each length-to-diameter ratio of
`aspect_ratio_sweep.py`:

```text
L/D     sensitivity (A per Bq/m3)   wall efficiency   linear up to (Ci/m3)
0.4     2.6425e-20                  0.9886            0.201
1.27    2.6592e-20                  0.9894            7.82
7       2.6470e-20                  0.9851            169
saved results/shape_comparison.png
```

```{figure} ../images/guide_shape_comparison.png
:alt: The current relative to the linear response against the activity, for three chamber shapes. The flat chamber leaves the linear response first and the long one last.

The figure it saves: the current relative to the linear response against the
activity, for each shape. The dashed line is 5% below linear.
```

## From a current to an activity

A monitor is used the other way round from the model: the current is measured and the
activity is wanted. `calibration.py` shows the two ways of doing that from a response
curve. The linear calibration divides the current by the sensitivity. Above the linear
range it reads too low, because part of the charge recombined, and reading the
activity off the full curve corrects for it.

```{literalinclude} ../../examples/calibration.py
:language: python
```

It prints the sensitivity, the two limits, and each measured current converted both
ways:

```text
sensitivity  2.659e-20 A per Bq/m3
             0.984 fA per uCi/m3
response 5% below linear at  2.92e+11 Bq/m3 (7.9 Ci/m3)
space-charge field reaches 0.1 at  1.51e+12 Bq/m3 (40.8 Ci/m3)

current (A)   activity, linear (Ci/m3)   activity, full curve (Ci/m3)
1.0e-12       0.001016                 0.001016
1.0e-09       1.016                    1.025
1.0e-08       10.16                    10.92
5.0e-08       50.82                    67.25

saved results/calibration.png
```

```{figure} ../images/guide_calibration.png
:alt: Two plots against the activity. On the left the chamber current follows the linear calibration and then falls below it. On the right the ratio of the two stays at 1 and then drops.

The figure it saves. Left: the model current and the linear calibration. Right: their
ratio, with the 5% limit as a dashed line and the activity where it is reached as a
dotted line.
```

The figures on this page were drawn with the colours of this documentation; the
scripts as written use matplotlib's defaults.
