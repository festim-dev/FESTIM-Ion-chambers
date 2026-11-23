# Quickstart

A model is put together from four things: a geometry, a gas, one or more ionisation
sources and a bias. Everything is in SI units: metres, volts, pascals, kelvin, Bq/m³.

## 1. The geometry

```python
import ion_chambers as ic

chamber = ic.CoaxialChamber(rod_radius=1.5e-3, can_radius=0.05, length=0.127)
```

## 2. The gas

The package ships no gas data, so you say what carries charge in your gas and how the
carriers react. This is room air, with illustrative values:

```python
electrons = ic.ChargeCarrier(name="e", charge=-1, mobility=0.04, diffusivity=0.04)
positive_ions = ic.ChargeCarrier(name="P", charge=1, mobility=2.0e-4)
negative_ions = ic.ChargeCarrier(name="N", charge=-1, mobility=2.2e-4)

attachment = ic.CarrierReaction(reactant=electrons, product=negative_ions, rate=1.0e8)
ion_ion_recombination = ic.CarrierReaction(
    reactant=[positive_ions, negative_ions], rate=2.0e-12
)
electron_ion_recombination = ic.CarrierReaction(
    reactant=[positive_ions, electrons], rate=2.0e-13
)

air = ic.Gas(
    name="air",
    w_value=33.97,  # eV per ion pair
    density=1.205,  # kg/m3
    pressure=101325.0,  # Pa
    temperature=293.15,  # K
    carriers=[electrons, positive_ions, negative_ions],
    ionisation_products=(electrons, positive_ions),
    reactions=[attachment, ion_ion_recombination, electron_ion_recombination],
)
```

[Gas](guide/gas.md) explains each entry, and
[Example: air](example_air.md) has a fuller version that follows the pressure.

## 3. The source

```python
tritium = ic.TritiumInGas(activity=3.7e10)  # Bq/m3, which is 1 Ci/m3
```

## 4. The model

```python
model = ic.IonChamberModel(
    geometry=chamber,
    gas=air,
    sources=tritium,
    bias=500.0,  # V
)
model.initialise()
reading = model.run()

print(f"current               {reading.current:.3e} A")
print(f"collection efficiency {reading.collection_efficiency:.4f}")
print(f"wall efficiency       {reading.wall_efficiency:.4f}")
```

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
    tritium.activity = activity
    reading = model.run()
    print(activity, reading.current)
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
