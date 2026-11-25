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

```{raw} html
<div class="ic-figure ic-wide">
<svg viewBox="0 0 660 340" role="img" aria-label="On the left, a closed can with a rod on its axis, with one half section shaded. On the right, that half section in the r, z plane with the anode, the cathode, the insulating rings and the four dimensions of the chamber.">
        <rect x="110" y="70" width="70" height="180" fill="var(--fd-accent)" fill-opacity="0.3"/>
        <ellipse cx="110" cy="70" rx="70" ry="18" fill="none" stroke="var(--ic-ink)" stroke-width="2"/>
        <path d="M40 70 L40 250 A70 18 0 0 0 180 250 L180 70" fill="none" stroke="var(--ic-ink)" stroke-width="2"/>
        <path d="M40 250 A70 18 0 0 1 180 250" fill="none" stroke="var(--ic-muted)" stroke-width="1" stroke-dasharray="3 4"/>
        <path d="M110 70 L110 250" stroke="var(--ic-ink)" stroke-width="5" fill="none"/>
        <text x="110" y="296" font-size="13" fill="var(--ic-ink)" text-anchor="middle">a can with a rod on its axis</text>
        <text x="110" y="313" font-size="12" fill="var(--ic-muted)" text-anchor="middle">the mesh covers the shaded half section</text>
        <path d="M200 160 L252 160 M244 154 L252 160 L244 166" stroke="var(--ic-ink)" stroke-width="2" fill="none"/>

        <path d="M300 30 L300 290" stroke="var(--ic-muted)" stroke-width="1" stroke-dasharray="10 4 2 4" fill="none"/>
        <text x="294" y="44" font-size="12" fill="var(--ic-muted)" text-anchor="end">axis</text>
        <rect x="300" y="50" width="30" height="220" fill="var(--ic-ink)" fill-opacity="0.15"/>
        <rect x="330" y="50" width="230" height="220" fill="var(--fd-accent)" fill-opacity="0.3"/>
        <path d="M360 50 L560 50 L560 270 L360 270" fill="none" stroke="var(--ic-ink)" stroke-width="2"/>
        <path d="M330 50 L360 50 M330 270 L360 270" fill="none" stroke="var(--ic-muted)" stroke-width="2" stroke-dasharray="4 3"/>
        <path d="M330 50 L330 270" stroke="var(--ic-ink)" stroke-width="5" fill="none"/>

        <text x="342" y="150" font-size="13" fill="var(--ic-ink)" style="font-family: var(--font-stack--monospace)">chamber.anode</text>
        <text x="342" y="167" font-size="12" fill="var(--ic-muted)">the surface of the rod</text>
        <text x="548" y="76" font-size="13" fill="var(--ic-ink)" style="font-family: var(--font-stack--monospace)" text-anchor="end">chamber.cathode</text>
        <text x="548" y="93" font-size="12" fill="var(--ic-muted)" text-anchor="end">the can and the end caps</text>
        <text x="548" y="206" font-size="13" fill="var(--ic-ink)" style="font-family: var(--font-stack--monospace)" text-anchor="end">chamber.gas</text>
        <text x="548" y="223" font-size="12" fill="var(--ic-muted)" text-anchor="end">the volume that is meshed</text>
        <text x="342" y="258" font-size="12" fill="var(--ic-muted)">insulating ring</text>

        <path d="M330 36 L360 36 M330 31 L330 41 M360 31 L360 41" stroke="var(--ic-muted)" stroke-width="1" fill="none"/>
        <text x="368" y="40" font-size="12" fill="var(--ic-ink)" style="font-family: var(--font-stack--monospace)">feedthrough_gap</text>
        <path d="M300 288 L330 288 M300 283 L300 293 M330 283 L330 293" stroke="var(--ic-muted)" stroke-width="1" fill="none"/>
        <text x="338" y="292" font-size="12" fill="var(--ic-ink)" style="font-family: var(--font-stack--monospace)">rod_radius</text>
        <path d="M300 312 L560 312 M300 307 L300 317 M560 307 L560 317" stroke="var(--ic-muted)" stroke-width="1" fill="none"/>
        <text x="430" y="330" font-size="12" fill="var(--ic-ink)" style="font-family: var(--font-stack--monospace)" text-anchor="middle">can_radius</text>
        <path d="M580 50 L580 270 M575 50 L585 50 M575 270 L585 270" stroke="var(--ic-muted)" stroke-width="1" fill="none"/>
        <text x="590" y="164" font-size="12" fill="var(--ic-ink)" style="font-family: var(--font-stack--monospace)">length</text>

        <path d="M588 300 L628 300 M621 295 L628 300 L621 305" stroke="var(--ic-ink)" stroke-width="1.5" fill="none"/>
        <path d="M588 300 L588 262 M583 269 L588 262 L593 269" stroke="var(--ic-ink)" stroke-width="1.5" fill="none"/>
        <text x="634" y="304" font-size="12" fill="var(--ic-ink)">r</text>
        <text x="596" y="266" font-size="12" fill="var(--ic-ink)">z</text>
      </svg>
</div>
```

The chamber is a body of revolution, so the model is 2D axisymmetric: the mesh covers
one (r, z) half section of the gas, from the surface of the rod to the can, and every
integral carries the factor 2πr that turns the section back into the full volume. The
rod itself is not meshed. The insulating rings are drawn wider than they are.

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

The 1 L chamber of the [example case](../example_case.md) is defined in
[its script](../example_case.md#the-script).

### Mesh

`radial_cells` is the number of cells from the rod to the can and `axial_cells` the
number from one end cap to the other. The cells are of equal size in each direction.

To look at the mesh, plot it with matplotlib. The coordinates are (r, z) in metres:

```python
import matplotlib.pyplot as plt

mesh = chamber.create_mesh(radial_cells=30, axial_cells=18)
points = mesh.mesh.geometry.x
triangles = mesh.mesh.geometry.dofmaps[0]

plt.triplot(points[:, 0], points[:, 1], triangles, linewidth=0.5)
plt.gca().set_aspect("equal")
plt.xlabel("r (m)")
plt.ylabel("z (m)")
plt.show()
```

```{figure} ../images/guide_mesh.png
:alt: A triangular mesh of the half section of the chamber, from the rod on the left to the can on the right.
:width: 60%

The mesh of the chamber above with 30 radial and 18 axial cells.
```

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
- Only the coaxial chamber exists so far.
- The numerical helpers assume a 2D triangular mesh, so a 3D geometry needs more than
  a new subclass.
- Wall losses are not computed from a mesh: a shape without a `Geometry` class runs
  without them. See the [roadmap](../roadmap.md).
:::
