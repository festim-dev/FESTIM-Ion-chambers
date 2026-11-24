# Results and exports

Results are obtained the FESTIM way. You create the quantities you want, each with the
surface or the volume it is computed on, add them to the `exports` of the transport
problem, and read them after the run.

```python
current = ic.ChamberCurrent(surface=chamber.anode)
collection_efficiency = ic.CollectionEfficiency(volume=chamber.gas)

transport.exports = [current, collection_efficiency]
transport.initialise()
transport.run()

current.value  # the value of the last run
current.data  # the values of all the runs so far
```

To write a quantity to a file, give it a filename ending in `.csv` or `.txt`:

```python
current = ic.ChamberCurrent(surface=chamber.anode, filename="results/current.csv")
```

The quantities are children of FESTIM's `DerivedQuantity`. A FESTIM quantity is
computed from the field of one species, as in `festim.SurfaceFlux(field=..., surface=...)`.
The quantities here involve every charge carrier and the potential, so they take a
`surface` or a `volume` but no field.

## What the chamber reads

{class}`~ion_chambers.ChamberCurrent`
: The chamber current in amperes, from the Shockley–Ramo theorem. This is the number
  to compare with a measurement. `surface` is the electrode the current is read from.
  The current is positive when that electrode is the more positive one.

{class}`~ion_chambers.IdealCurrent`
: The current if every beta stopped in the gas and every charge reached an electrode:
  `e × activity × volume × pairs per decay`. `volume` is the gas, as for the
  quantities below.

The ratio of the two is the product of the two efficiencies below, up to
discretisation error.

{class}`~ion_chambers.WallEfficiency`
: The share of the beta energy deposited in the gas rather than in the walls.

{class}`~ion_chambers.CollectionEfficiency`
: The share of the ion pairs that escape recombination.

{class}`~ion_chambers.PairGeneration`, {class}`~ion_chambers.Recombination`
: The underlying rates, in ion pairs per second over the whole gas volume.

## How far to trust it

{class}`~ion_chambers.ElectrodeCurrent`
: The charge leaving the gas through an electrode per second, from the carrier fluxes
  at that surface. In exact arithmetic its size equals the chamber current. The
  difference shows the discretisation error of the mesh you used; the two typically
  agree to better than 0.01%. It is negative on the electrode that collects the
  negative carriers. For the flux of one carrier alone, use FESTIM's
  `F.SurfaceFlux(field=electrons, surface=chamber.anode)`, in carriers per second.

{class}`~ion_chambers.SpaceChargeField`
: The field the charges in the gas would create, relative to the applied field. The
  transport problem assumes this is small. Above about 0.1 the result is indicative
  only. See [Space charge](electrostatics.md#space-charge).

## Why the Ramo current

The current through an electrode can be computed directly from the flux of carriers
into it, but that needs the gradient of the densities at the boundary, where the
numerical error is largest. The Shockley–Ramo theorem gives the same current as an
integral over the whole gas volume, which is much less sensitive to the mesh.

## Fields

A `ChargeCarrier` is a FESTIM species, so FESTIM's field export writes its density to
a file that opens in ParaView:

```python
import festim as F

transport.exports = [
    current,
    F.SpeciesExport(field=electrons, filename="results/electrons.bp"),
]
```

The potential is exported from the electrostatic problem with
{class}`~ion_chambers.PotentialExport`.

:::{note}
Of FESTIM's own derived quantities, only `festim.SurfaceFlux` supports the cylindrical
meshes used here. FESTIM refuses the others, such as `festim.TotalVolume`, on a
cylindrical mesh.
:::

## Adding a quantity

A new quantity is a child of {class}`~ion_chambers.ChamberSurfaceQuantity` or
{class}`~ion_chambers.ChamberVolumeQuantity`, with a `title` and a
`compute(problem)` method that stores the number in `value` and appends it to
`data`, as FESTIM's quantities do. For instance, the number of
carriers of every kind in a volume:

```python
from scifem import assemble_scalar


class TotalCarriers(ic.ChamberVolumeQuantity):
    """The number of charge carriers in a volume."""

    @property
    def title(self):
        return f"Total carriers volume {self.volume.id}"

    def compute(self, problem):
        # the factor 2 pi r of a cylindrical mesh
        weight = problem.mesh.coordinate_system.integration_weight(problem.mesh.mesh)
        dx = problem.dx(self.volume.id)

        total_density = 0
        for carrier in problem.gas.carriers:
            total_density += carrier.post_processing_solution
        self.value = assemble_scalar(weight * total_density * dx)
        self.data.append(self.value)


total_carriers = TotalCarriers(volume=chamber.gas)
```
