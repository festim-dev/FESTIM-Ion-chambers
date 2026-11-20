# Quickstart

A model is put together from four objects: a geometry, a gas, one or more ionisation
sources and a bias.

```python
import ion_chambers as ic

model = ic.IonChamberModel(
    geometry=ic.CoaxialChamber(rod_radius=1.5e-3, can_radius=0.05, length=0.127),
    gas=ic.library.air(),
    sources=ic.TritiumInGas(activity=3.7e10),  # Bq/m3, i.e. 1 Ci/m3
    bias=500.0,  # V
)
model.initialise()
reading = model.run()

print(f"current               {reading.current:.3e} A")
print(f"collection efficiency {reading.collection_efficiency:.4f}")
print(f"wall efficiency       {reading.wall_efficiency:.4f}")
```

Everything is in SI units: metres, volts, pascals, kelvin, Bq/m³.

## The two steps

Nothing is computed when you create the objects. They only hold your input.

`model.initialise()`
: builds the mesh, solves the electric potential, runs the wall-loss Monte Carlo and
  assembles the FESTIM problem. This takes a few seconds.

`model.run()`
: solves the steady state and returns a {class}`~ion_chambers.Reading`. This takes
  about a second.

## Changing the activity

The activity can be changed without initialising again, which makes a response curve
cheap:

```python
for activity in [3.7e8, 3.7e10, 3.7e12]:
    model.sources[0].activity = activity
    print(activity, model.run().current)
```

Changing anything else (geometry, gas, bias, resolution) needs a new call to
`model.initialise()`.

## Looking at the fields

After a run, the densities of the charge carriers and the potential are DOLFINx
functions:

```python
densities = model.concentrations()  # {"e": ..., "P": ..., "N": ...} in 1/m3
potential = model.potential  # V
```

## Mesh resolution

`resolution` is the number of cells across the electrode gap (120 by default). The
current converges at first order in the cell size, so a single mesh is not a converged
answer. For numbers you intend to quote, extrapolate from two meshes: see
[Parameter studies](guide/studies.md).
