# Reference case: a chamber in air

One case runs through this documentation and the scripts in `examples/`: a generic
1 L coaxial chamber filled with room air. This page shows what the model gives for
it, then the script that defines it.

The case has not been validated against measurements. It shows what the model
produces and that it behaves sensibly. It is not a prediction for a real monitor.

## The fields

The chamber is a rod of radius 1.5 mm in a can of radius 5.0 cm and length 12.7 cm
(length/diameter 1.27), with 3 mm insulating rings and the anode at +500 V.

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
The figures on this page were drawn with the first version of the code, the
prototype that showed this model could be built on FESTIM. The package reproduces the
currents of that prototype on the same mesh, and the test suite checks this. At
100 Ci/m³ the space charge of the ions would distort the field, which the model does
not apply, so the maps show the pattern rather than exact values.
:::

## The script

The case is the file `examples/reference_case.py`. It defines the air, the chamber,
and a function that solves the potential and sets up the transport problem. The other
scripts in `examples/` import it.

The package does not ship any gas data. You describe the gas yourself, so that every
number in a result is one you chose and can check. Copy this file as a starting point
for your own gas and chamber.

:::{warning}
The values below are reasonable starting points, not validated reference data. The
way they are scaled with pressure and temperature is deliberately simple. Check both
against your conditions before relying on a result.
:::

```{literalinclude} ../examples/reference_case.py
:language: python
```

### Values at 1 atm and 293.15 K

| Quantity | Value | Unit |
|---|---|---|
| W value | 33.97 | eV |
| Density | 1.205 | kg/m³ |
| Positive ion mobility | 2.0 × 10⁻⁴ | m²/V/s |
| Negative ion mobility | 2.2 × 10⁻⁴ | m²/V/s |
| Electron mobility (low field) | 0.04 | m²/V/s |
| Electron characteristic energy D/μ | 1.0 | V |
| Attachment frequency to O₂ | 1.0 × 10⁸ | 1/s |
| Ion–ion recombination | 2.0 × 10⁻¹² | m³/s |
| Electron–ion recombination | 2.0 × 10⁻¹³ | m³/s |

### How the air scales with pressure and temperature

With `n` the gas density relative to 1 atm and 293.15 K:

- the mass density scales as `n`
- the mobilities scale as `1 / n`
- the attachment frequency scales as `n²`, because attachment to O₂ is a three-body
  process
- the recombination coefficients and the electron characteristic energy are kept
  constant, which does not hold over a wide pressure range

### Using it

From a script in the `examples/` folder:

```python
from reference_case import air, reference_chamber

import ion_chambers as ic

chamber = reference_chamber()
gas = air(pressure=50_000.0)  # Pa
```

## Studies of this case

The numbers below come from the [examples](guide/studies.md#ready-made-examples):
each section names the two scripts that produce it.

### Response to tritium concentration

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

### Chamber shape

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

### Bias

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
