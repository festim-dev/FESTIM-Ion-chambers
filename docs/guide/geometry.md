# Geometry

A geometry tells the model four things: how to mesh the gas volume, where the two
electrodes are, how large the gas volume is, and how far a beta can travel before it
hits a wall.

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

Two shortcuts fix the gas volume and let a dimension follow:

```python
# given radii, the length follows
chamber = ic.CoaxialChamber.from_volume(volume=1e-3, rod_radius=1.5e-3, can_radius=0.05)

# given the length-to-diameter ratio, the radius and length follow
chamber = ic.CoaxialChamber.from_aspect_ratio(
    volume=1e-3, rod_radius=1.5e-3, aspect_ratio=7.0
)
```

`ic.library.reference_chamber()` returns the 1 L chamber used for the
[prototype results](../results.md).

### Mesh

The mesh is graded so that cells are small where the field is strong: the radial
spacing grows in proportion to r, and the axial spacing is refined towards the end
caps. The model's `resolution` is the number of radial cells. The number of axial cells
follows from it unless you set `axial_cells`.

### Bias polarity

The rod is the electrode held at the bias and the can is grounded. A negative `bias`
makes the rod collect the positive ions instead.

## Adding a geometry

A new shape is a subclass of {class}`~ion_chambers.Geometry` that implements
`build_mesh`, `is_anode`, `is_cathode`, `volume`, `distance_to_wall` and
`wall_clustered_grid`. The rest of the package only uses those methods.

:::{admonition} Limitations
:class: note
- Only the coaxial chamber exists so far. The real 1 L monitor (a short, thick anode)
  and the 20 cc wire cage are not modelled yet.
- The numerical helpers assume a 2D triangular mesh, so a 3D geometry needs more than
  a new subclass.
- Reading a mesh from a file is planned: see the [roadmap](../roadmap.md).
:::
