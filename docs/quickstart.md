# Quickstart

A calculation has two parts, used like FESTIM problems:

1. an {class}`~ion_chambers.ElectrostaticProblem` gives the electric potential between
   the electrodes
2. a {class}`~ion_chambers.ChargeTransportProblem` moves the charge carriers of the gas
   in that potential and gives the chamber current

Everything is in SI units: metres, volts, pascals, kelvin, Bq/m³.

## 1. The geometry

```python
import ion_chambers as ic

chamber = ic.CoaxialChamber(rod_radius=1.5e-3, can_radius=0.05, length=0.127)
```

The chamber provides a mesh and three FESTIM subdomains: `chamber.gas`,
`chamber.anode` (the rod) and `chamber.cathode` (the can).

## 2. The electric potential

```python
electrostatics = ic.ElectrostaticProblem()
electrostatics.mesh = chamber.create_mesh(radial_cells=120, axial_cells=72)
electrostatics.subdomains = [chamber.gas, chamber.anode, chamber.cathode]
electrostatics.boundary_conditions = [
    ic.FixedPotentialBC(subdomain=chamber.anode, value=500.0),  # V
    ic.FixedPotentialBC(subdomain=chamber.cathode, value=0.0),
]
electrostatics.initialise()
electrostatics.run()
```

## 3. The gas

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
the [example case](example_case.md#the-script) has a fuller version that follows
the pressure.

## 4. The source

```python
tritium = ic.TritiumInGas(activity=3.7e10, geometry=chamber)  # Bq/m3, i.e. 1 Ci/m3
```

Giving the geometry makes the source account for the beta energy lost to the walls.

## 5. The transport of the charge carriers

Results are obtained through exports, as in FESTIM: choose the quantities you want,
add them to the problem, and read their `value` after the run.

```python
current = ic.ChamberCurrent(surface=chamber.anode)
collection_efficiency = ic.CollectionEfficiency(volume=chamber.gas)
wall_efficiency = ic.WallEfficiency(volume=chamber.gas)

transport = ic.ChargeTransportProblem()
transport.electrostatic_problem = electrostatics
transport.gas = air
transport.ionisation_sources = [tritium]
transport.exports = [current, collection_efficiency, wall_efficiency]
transport.initialise()
transport.run()

print(f"current               {current.value:.3e} A")
print(f"collection efficiency {collection_efficiency.value:.4f}")
print(f"wall efficiency       {wall_efficiency.value:.4f}")
```

`transport.initialise()` takes a few seconds, most of it for the wall-loss
calculation. `transport.run()` takes about a second.

## Changing the activity

The activity can be changed without initialising again, which makes a response curve
cheap:

```python
for activity in [3.7e8, 3.7e10, 3.7e12]:
    tritium.activity = activity
    transport.run()
    print(activity, current.value)
```

Every run adds a value to `current.data`.

Changing anything else (geometry, gas, potential, mesh) needs new problems. A
problem is initialised once.

## Looking at the fields

A `ChargeCarrier` is a FESTIM species, so after a run its density is available the
FESTIM way, and the potential is the solution of the electrostatic problem:

```python
electron_density = electrons.post_processing_solution  # 1/m3
potential = electrostatics.potential  # V
```

To write them to files for ParaView, see [Results and exports](guide/exports.md).

## Mesh resolution

`radial_cells` is the number of mesh cells from the rod to the can and `axial_cells`
the number from one end cap to the other. The current converges at first order in the cell size, so a single mesh is not a converged answer.
For numbers you intend to quote, extrapolate from two meshes: see
[Parameter studies](guide/studies.md).
