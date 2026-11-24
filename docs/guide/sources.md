# Ionisation sources

A source says where ion pairs are created and how fast. A transport problem takes a
list of them in `ionisation_sources`, and adds their contributions.

## Tritium in the gas

{class}`~ion_chambers.TritiumInGas` is tritium mixed uniformly into the gas:

```python
import ion_chambers as ic

tritium = ic.TritiumInGas(activity=3.7e10, geometry=chamber)  # Bq/m3
```

To convert from curies, `ic.constants.BQ_PER_CI` is 3.7 × 10¹⁰.

Each decay releases a beta with 5.7 keV on average. If all of that energy stays in the
gas, it makes `mean energy / W` ion pairs: about 168 in air.

## Wall losses

A tritium beta travels up to about 5 mm in air at 1 atm. A decay close to a wall
therefore leaves part of its energy in the metal, where it makes no signal. The source
accounts for this with a local efficiency f(r, z) between 0 and 1, computed by a Monte
Carlo: betas are drawn from the tritium spectrum, sent in random directions, and
followed in straight lines until they stop or hit a wall.

For the reference 1 L chamber in air, 98.9% of the beta energy stays in the gas. The
loss grows when the gas is thinner (lower pressure, lighter gas) or the chamber
smaller, because the range of the betas scales inversely with the gas density.

The source needs the chamber geometry for this, which is why it is given to the
source. Without a geometry, wall losses are ignored:

```python
ic.TritiumInGas(activity=3.7e10)  # every beta stops in the gas
ic.TritiumInGas(activity=3.7e10, geometry=chamber, n_rays=10000)  # more betas per point
```

The Monte Carlo is the slowest part of `transport.initialise()`. Its result depends
only on the shape of the chamber and the gas density, and is kept in memory, so
sweeping the bias or the activity pays for it once.

The method's assumptions are listed in [Theory](../theory.md#wall-losses).

## Changing the activity

```python
tritium.activity = 3.7e11
transport.run()
```

No call to `initialise()` is needed in between.

## Adding a source

A new source is a subclass of {class}`~ion_chambers.IonisationSource`. The planned one
is tritium held on the chamber walls, which is what gives a monitor its memory after an
exposure: see the [roadmap](../roadmap.md).
