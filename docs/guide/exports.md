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

The equations below use the symbols of the
[transport problem](transport.md#what-it-solves): $n_s$, $J_s$ and $z_s$ are the
density, the flux and the charge number of carrier $s$, $\varphi$ is the potential,
$S$ the ionisation source and $e$ the elementary charge. $\Omega$ is the `volume` of
the quantity and $\Gamma$ its `surface`.

### Chamber current

{class}`~ion_chambers.ChamberCurrent` is the chamber current in amperes, from the
Shockley–Ramo theorem:

$$
I = e \int_\Omega \sum_s z_s \, J_s \cdot \left( -\nabla \varphi_w \right) dV,
\qquad \varphi_w = \frac{\varphi}{V_\mathrm{readout} - V_\mathrm{other}}
$$

$\varphi_w$ is the weighting potential, which changes by 1 from the other electrode
to the readout electrode. `surface` is the electrode the current is read from. The
current is positive when that electrode is the more positive one. This is the number
a monitor shows.

### Ideal current

{class}`~ion_chambers.IdealCurrent` is the current if every beta stopped in the gas
and every charge reached an electrode:

$$
I_\mathrm{ideal} = e \, A \, \frac{\bar{E}}{W} \, V_\mathrm{gas}
$$

with $A$ the activity concentration, $\bar{E} / W$ the ion pairs per decay and
$V_\mathrm{gas}$ the size of the volume.

### Pair generation and recombination

{class}`~ion_chambers.PairGeneration` and {class}`~ion_chambers.Recombination` are the
ion pairs created and lost per second over the volume:

$$
G = \int_\Omega S \, dV, \qquad
R = \int_\Omega \sum_j k_j \prod_{s \in j} n_s \, dV
$$

The sum runs over the reactions $j$ without a product, the recombinations, each with
its rate coefficient $k_j$ and the densities of its reactants.

### Efficiencies

{class}`~ion_chambers.WallEfficiency` is the share of the beta energy deposited in
the gas rather than in the walls, and {class}`~ion_chambers.CollectionEfficiency` the
share of the ion pairs that escape recombination:

$$
\eta_\mathrm{wall} = \frac{e \, G}{I_\mathrm{ideal}}, \qquad
\eta_\mathrm{collection} = \frac{G - R}{G}
$$

The chamber current is the ideal current times the two efficiencies, up to
discretisation error:

$$
I = \eta_\mathrm{wall} \, \eta_\mathrm{collection} \, I_\mathrm{ideal}
$$

## How far to trust it

### Electrode current

{class}`~ion_chambers.ElectrodeCurrent` is the charge leaving the gas through an
electrode per second, from the carrier fluxes at that surface:

$$
I_\Gamma = e \int_\Gamma \sum_s z_s \, J_s \cdot n \, dS
$$

with $n$ the normal pointing out of the gas. In exact arithmetic its size equals the
chamber current. The difference shows the discretisation error of the mesh you used;
for the reference chamber the two agree to about 0.1% at 120 radial cells and 0.02%
at 300. It is negative on the electrode that collects the negative carriers. For the
flux of one carrier alone, use FESTIM's
`F.SurfaceFlux(field=electrons, surface=chamber.anode)`, in carriers per second.

### Space-charge field

{class}`~ion_chambers.SpaceChargeField` is the field the charges in the gas would
create, relative to the applied field:

$$
\frac{\lVert \nabla \varphi_1 \rVert}{\lVert \nabla \varphi \rVert},
\qquad
-\nabla \cdot (\varepsilon_0 \nabla \varphi_1) = e \sum_s z_s \, n_s,
\quad \varphi_1 = 0 \text{ on the electrodes}
$$

with $\lVert \cdot \rVert$ the L2 norm over the volume. The transport problem
assumes this ratio is small. Above about 0.1 the result is indicative only. See
[Space charge](electrostatics.md#space-charge).

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
