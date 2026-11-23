# Parameter studies

The package has no function for parameter studies. A study is a loop in your own
script, which keeps every step visible and easy to change. This page shows the two
patterns you need, and the scripts in `scripts/` are complete examples to copy.

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

`scripts/response_curve.py` puts the two patterns together: the same chamber on a
coarse and on a fine mesh, solved at every activity and extrapolated. It uses
`build_transport_problem` from the [worked example in air](../example_air.md) to set
up each problem.

```{literalinclude} ../../scripts/response_curve.py
:language: python
```

## Ready-made scripts

The `scripts/` folder holds the studies behind the [prototype results](../results.md).
Each writes a JSON file to `results/`.

```bash
python scripts/response_curve.py       # current against activity
python scripts/bias_sweep.py           # saturation curve
python scripts/aspect_ratio_sweep.py   # chamber shape at fixed volume
```

Each script is a plain sequence of steps with its parameters at the top. Edit those to
change the case. With the default meshes of 150 and 300 radial cells, expect seconds to
minutes per case.
