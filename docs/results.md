# Results so far

These results are for a generic chamber in air and have not been validated against
measurements. They show what the model produces and that it behaves sensibly. They are
not predictions for a real monitor.

The numbers on this page come from the [examples](guide/studies.md#ready-made-examples):
each section names the two scripts that produce it.

## Reference case

A 1 L coaxial chamber in room air: rod radius 1.5 mm, can radius 5.0 cm, length
12.7 cm (length/diameter 1.27), 3 mm insulating rings, anode at +500 V.
Both are defined in [Example: air and the reference chamber](example_air.md).

## Fields

```{figure} images/prototype_fields.png
:alt: Five maps in the half section of the chamber. The electric field, the density of free electrons, of positive ions and of negative ions, and the recombination rate.

The reference chamber at 100 Ci/m³, in the half section (r, z). The axis is on the
left and the anode rod is at r = 0.15 cm. The lines are electric field lines. Black
lines mark the can and the end caps, gold segments the insulating rings.
```

The field is strongest at the rod and weakest in the corners between the can and the
end caps. Electrons attach to oxygen within about 10 ns, a few micrometres from where
they are freed, so their density mirrors the source and is nearly uniform. Positive
ions pile up towards the can and negative ions towards the rod. The two overlap most
in the middle of the gap, and that is where they recombine.

:::{note}
The two figures on this page were made with the first version of the code, the
prototype that showed this model could be built on FESTIM. The package reproduces the
currents of that prototype on the same mesh, and the test suite checks this. At
100 Ci/m³ the space charge of the ions would distort the field, which the model does
not apply, so the maps show the pattern rather than exact values.
:::

## Response to tritium concentration

From `response_curve.py` and `calibration.py`.

| Quantity | Value |
|---|---|
| Sensitivity in the linear range | 2.66 × 10⁻²⁰ A per Bq/m³ |
| Beta energy deposited in the gas | 98.9% |
| Response 5% below linear at | 7.9 Ci/m³ |
| Space-charge field reaches 10% of the applied field at | 41 Ci/m³ |

The current is linear in the activity over eight decades, about 1% below the ideal
value because of wall losses. Ion–ion recombination then bends the curve. Beyond
41 Ci/m³ the fixed-field assumption no longer holds and the results are indicative
only.

## Chamber shape

From `aspect_ratio_sweep.py` and `shape_comparison.py`. At a fixed gas volume of 1 L
and a 500 V bias, the length-to-diameter ratio L/D sets a trade-off:

- Flat chambers have a long drift gap and corners that the field barely reaches. Ions
  linger there and recombine.
- Long, thin chambers collect quickly, which delays recombination, at the cost of
  slightly more beta energy lost to the can.

```{figure} images/prototype_shapes.png
:alt: Maps of the positive-ion density in three chambers of the same volume, one flat, one of the reference shape and one long and thin.

Positive-ion density at 100 Ci/m³ in three 1 L chambers drawn to one scale, with the
axis along the bottom of each. The lines are field lines from the rod. In the flat
chamber they bow out to the end caps and leave the outer corners nearly field-free.
```

| L/D | Response 5% below linear at |
|---|---|
| 0.4 | 0.2 Ci/m³ |
| 1.27 | 7.8 Ci/m³ |
| 7.0 | 169 Ci/m³ |

The linear range is 22 times wider at L/D = 7 than at 1.27. Trace-level sensitivity
peaks at L/D = 1.27, and the longest chamber gives up only 0.5% of it. (The shape
sweep uses fewer activities than the response curve, hence 7.8 here and 7.9 above.)

## Bias

From `bias_sweep.py` and `saturation_curve.py`. Saturation curves of the reference
chamber:

| Activity | Bias for 99% collection |
|---|---|
| 1 Ci/m³ | 432 V |
| 10 Ci/m³ | 1346 V |
| 100 Ci/m³ | more than 2000 V |

At 10 Ci/m³ the 500 V design point collects 94.0%. Fitting 1/I against 1/V² over the
model points from 500 to 2000 V, the way a measured saturation curve is analysed,
recovers the saturation current to 0.07% at 1 and 10 Ci/m³. At 100 Ci/m³, where the
chamber is far from saturation at these biases, the fit is 2.9% low.

## Checks

These are tests of the package (`tests/test_verification.py`), not comparisons with
measurements.

Parallel plates
: For a parallel-plate chamber with a uniform source and two kinds of ions, the share
  of the pairs lost to recombination agrees with Boag's formula to 1%, and the current
  extrapolated from two meshes agrees to 0.01%. A transient run ends at the
  steady-state current.

Manufactured solution
: With densities chosen in advance and the sources worked out from them, the solved
  densities converge to the chosen ones at second order when the diffusion floor is
  inactive, and at first order when it sets the diffusion, on a cartesian and on a
  cylindrical mesh.

Mesh convergence
: The error is first order in the cell size. Each quoted current is extrapolated from
  150 and 300 radial cells. The tests check the order on the parallel-plate and the
  coaxial chamber, and that the extrapolated coaxial current at trace level is the
  charge of all the pairs created.
