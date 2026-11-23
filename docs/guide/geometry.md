# Geometry

A geometry provides the mesh of the gas volume, the subdomains of the gas and of the
two electrodes, and the distance a beta can travel before it hits a wall.

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
mesh = chamber.create_mesh(resolution=120)  # a festim.Mesh, cylindrical
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

The mesh is graded so that cells are small where the field is strong: the radial
spacing grows in proportion to r, and the axial spacing is refined towards the end
caps. `resolution` is the number of radial cells. The number of axial cells follows
from it unless you set `axial_cells`.

## Adding a geometry

A new shape is a subclass of {class}`~ion_chambers.Geometry` that provides the `gas`,
`anode` and `cathode` subdomains and implements `create_mesh`, `distance_to_wall` and
`wall_clustered_grid`.

The problems themselves only need a `festim.Mesh` and FESTIM subdomains, so a mesh
from another source can be used directly. The geometry class is what the wall-loss
calculation needs.

:::{admonition} Limitations
:class: note
- Only the coaxial chamber exists so far. The real 1 L monitor (a short, thick anode)
  and the 20 cc wire cage are not modelled yet.
- The numerical helpers assume a 2D triangular mesh, so a 3D geometry needs more than
  a new subclass.
- Reading a mesh from a file is planned: see the [roadmap](../roadmap.md).
:::
