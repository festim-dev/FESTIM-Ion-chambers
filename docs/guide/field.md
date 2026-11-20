# Electric field and numerics

## Potential solver

The charge carriers drift in an electric potential that is computed once, before the
transport problem. Which equation the potential solves is chosen with the `field`
argument of the model.

{class}`~ion_chambers.LaplaceField` is the default and, for now, the only choice. It
solves Laplace's equation with the bias on the anode and 0 V on the cathode: the field
of the empty chamber.

```python
import ion_chambers as ic

model = ic.IonChamberModel(..., field=ic.LaplaceField())
```

This ignores the field created by the charges themselves. That is a good approximation
at low activity and breaks down at high activity, where slow ions accumulate. Every
reading tells you where you stand: `reading.space_charge_field` is the size of the
field the space charge would create, relative to the applied field. Treat results as
indicative once it exceeds about 0.1. For the reference chamber at 500 V that happens
near 37 Ci/m³.

A solver that feeds the space charge back into the potential is on the
[roadmap](../roadmap.md). It will be another subclass of
{class}`~ion_chambers.PotentialSolver`.

## Stabilisation

Ions in a chamber are carried by the field far more than they diffuse. Standard finite
elements produce oscillations in that regime, and FESTIM does not stabilise it. The
model therefore raises the diffusion coefficient in each cell to at least

$$
D = \max\left(D_\mathrm{phys},\ \tfrac{1}{2}\,\mu\,|E|\,h_E\right)
$$

where $h_E$ is the size of the cell along the field. This is the classic upwind
scheme, written as extra diffusion. Two properties matter:

- The extra diffusion vanishes where the field is weak, so it does not create
  artificial losses in the low-field corners of a chamber.
- The error it leaves shrinks in proportion to the cell size. This is why currents are
  [extrapolated from two meshes](studies.md#mesh-extrapolation).

The strength can be changed, mainly to test how sensitive a result is:

```python
model = ic.IonChamberModel(..., stabilisation=ic.DiffusionFloor(factor=0.5))
```

## Boundary conditions

An electrode absorbs the carriers that the field drives into it and emits none. The
insulating rings let nothing through. These are set automatically from the sign of
each carrier's charge and the sign of the bias.
