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
change the case. With the default meshes of 150 and 300 radial cells, a sweep takes a
few minutes. The post-processing scripts take a second.

## From a current to an activity

A monitor is used the other way round from the model: the current is measured and the
activity is wanted. `calibration.py` shows the two ways of doing that from a response
curve. The linear calibration divides the current by the sensitivity. Above the linear
range it reads too low, because part of the charge recombined, and reading the
activity off the full curve corrects for it.

```{literalinclude} ../../examples/calibration.py
:language: python
```
