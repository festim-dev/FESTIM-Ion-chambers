# Geometry

The problems need a `festim.Mesh` of the gas, a volume subdomain for the gas and a
surface subdomain for each electrode. There are two ways to get them: let
{class}`~ion_chambers.CoaxialChamber` create them, or bring
[your own mesh](#using-your-own-mesh).

## Coaxial chamber

{class}`~ion_chambers.CoaxialChamber` is an anode rod on the axis of a closed
cylindrical can. The can and both end caps form the cathode. Where the rod passes
through a cap, an insulating ring separates the two.

```python
import ion_chambers as ic

chamber = ic.CoaxialChamber(
    rod_radius=1.5e-3,  # m
    can_radius=0.05,  # m
    length=0.127,  # m
    feedthrough_gap=3e-3,  # m, radial width of the insulating rings
)
print(chamber.volume)  # m3
```

The model is 2D axisymmetric: the mesh covers the (r, z) half section of the gas.

The chamber gives you what a FESTIM-style problem needs:

```python
mesh = chamber.create_mesh(radial_cells=120, axial_cells=72)  # a festim.Mesh
chamber.gas  # festim.VolumeSubdomain, the gas
chamber.anode  # festim.SurfaceSubdomain, the rod
chamber.cathode  # festim.SurfaceSubdomain, the can and end caps
```

Two shortcuts fix the gas volume and let a dimension follow:

```python
# given radii, the length follows
chamber = ic.CoaxialChamber.from_volume(volume=1e-3, rod_radius=1.5e-3, can_radius=0.05)

# given the length-to-diameter ratio, the radius and length follow
chamber = ic.CoaxialChamber.from_aspect_ratio(
    volume=1e-3, rod_radius=1.5e-3, aspect_ratio=7.0
)
```

The 1 L chamber used for the [prototype results](../results.md) is defined in
[Example: air and the reference chamber](../example_air.md).

### Mesh

`radial_cells` is the number of cells from the rod to the can and `axial_cells` the
number from one end cap to the other. The cells are of equal size in each direction.

## Using your own mesh

A mesh made with another tool is used exactly as in FESTIM. Read the mesh and its
cell and facet tags, for instance with
[io4dolfinx](https://github.com/scientificcomputing/io4dolfinx), and create one
subdomain for each tag you need. A subdomain only needs the id that the tag has in
the mesh, with no locator.

```python
from mpi4py import MPI

import festim as F
import io4dolfinx

import ion_chambers as ic

# the ids given to the gas and to the electrodes when the mesh was made
GAS_ID = 7
ANODE_ID = 11
CATHODE_ID = 12

mesh = io4dolfinx.read_mesh("chamber.bp", MPI.COMM_WORLD)
cell_tags = io4dolfinx.read_meshtags("chamber.bp", mesh, meshtag_name="cells")
facet_tags = io4dolfinx.read_meshtags("chamber.bp", mesh, meshtag_name="facets")

gas_volume = F.VolumeSubdomain(id=GAS_ID, material=None)
anode = F.SurfaceSubdomain(id=ANODE_ID)
cathode = F.SurfaceSubdomain(id=CATHODE_ID)

electrostatics = ic.ElectrostaticProblem()
electrostatics.mesh = F.Mesh(mesh, coordinate_system="cylindrical")
electrostatics.volume_meshtags = cell_tags
electrostatics.facet_meshtags = facet_tags
electrostatics.subdomains = [gas_volume, anode, cathode]
electrostatics.boundary_conditions = [
    ic.FixedPotentialBC(subdomain=anode, value=500.0),
    ic.FixedPotentialBC(subdomain=cathode, value=0.0),
]
electrostatics.initialise()
electrostatics.run()
```

The transport problem takes the mesh, the tags and the electrodes from the
electrostatic problem, so nothing more is needed there. Use your own surface
subdomains in the exports, for instance `ic.ChamberCurrent(surface=anode)`.

What the mesh must satisfy:

- It is 2D and made of triangles. For an axisymmetric chamber the coordinates are
  (r, z) and the coordinate system is `"cylindrical"`.
- Every cell belongs to the gas: the electrostatic problem has one volume subdomain.
- Every electrode has a facet tag and a `FixedPotentialBC`. Facets without a boundary
  condition, tagged or not, are insulators.

### Wall losses with your own mesh

The wall-loss calculation does not use the mesh. It needs a
{class}`~ion_chambers.Geometry` that knows where the walls are, given to the source as
`ic.TritiumInGas(activity=..., geometry=...)`. With your own mesh you can:

- leave `geometry` out, which ignores wall losses (about 1 % for the 1 L reference
  chamber in air at 1 atm, more for small chambers or low pressure);
- pass a `CoaxialChamber` of the same dimensions, if your chamber is close to one;
- write a subclass of `Geometry` for your shape.

## Adding a geometry

A new shape is a subclass of {class}`~ion_chambers.Geometry`. What the wall-loss
calculation needs from it is `dimensions`, `distance_to_wall` and
`wall_clustered_grid`. A geometry can also create its own mesh and subdomains, as
`CoaxialChamber` does, but it does not have to.

:::{admonition} Limitations
:class: note
- Only the coaxial chamber exists so far. The real 1 L monitor (a short, thick anode)
  and the 20 cc wire cage are not modelled yet.
- The numerical helpers assume a 2D triangular mesh, so a 3D geometry needs more than
  a new subclass.
- Wall losses are not computed from a mesh: a shape without a `Geometry` class runs
  without them. See the [roadmap](../roadmap.md).
:::
