# Results so far

These results are for a generic chamber in air and have not been validated against
measurements. They show what the model produces and that it behaves sensibly. They are
not predictions for a real monitor.

The full write-up with figures is the
<a href="prototype_report.html">prototype report</a> (October 2026).

:::{note}
The report was produced with the first version of the code, on FESTIM 2.2-rc.2. The
restructured package reproduces its currents on the same mesh, and the test suite
checks this.
:::

The numbers on this page were recomputed with the current package, on its uniform
mesh, with the [examples](guide/studies.md#ready-made-examples): each table names the
two scripts that produce it. They differ slightly from the report, which used a graded
mesh.

## Reference case

A 1 L coaxial chamber in room air: rod radius 1.5 mm, can radius 5.0 cm, length
12.7 cm (length/diameter 1.27), 3 mm insulating rings, anode at +500 V.
Both are defined in [Example: air and the reference chamber](example_air.md).

## Response to tritium concentration

From `response_curve.py` and `calibration.py`.

| Quantity | Value |
|---|---|
| Sensitivity in the linear range | 2.66 × 10⁻²⁰ A per Bq/m³ |
| Beta energy deposited in the gas | 98.9% |
| Response 5% below linear at | 7.9 Ci/m³ |
| Space-charge field reaches 10% of the applied field at | 41 Ci/m³ |

The current is linear in the activity over eight decades, 1.0% below the ideal value
because of wall losses. Ion–ion recombination then bends the curve. Beyond 41 Ci/m³
the fixed-field assumption no longer holds and the results are indicative only.

## Chamber shape

From `aspect_ratio_sweep.py` and `shape_comparison.py`. At a fixed gas volume of 1 L and a 500 V bias, the length-to-diameter ratio L/D sets a
trade-off:

- Flat chambers have a long drift gap and corners that the field barely reaches. Ions
  linger there and recombine.
- Long, thin chambers collect quickly, which delays recombination, at the cost of
  slightly more beta energy lost to the can.

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

Charge conservation
: The Ramo current and the integrated anode flux agree to better than 0.01%. At trace
  levels the extrapolated current matches the balance of generation minus
  recombination to 0.09%.

Recombination
: Boag's formula for a cylindrical chamber under continuous irradiation gives a
  collection efficiency of 0.779 at 100 Ci/m³ for the reference geometry. The model
  gives 0.705. The model's extra loss comes from the low-field corners near the end
  caps and the non-uniform source, which the one-dimensional formula ignores.

Parallel plates
: For a parallel-plate chamber with a uniform source and two kinds of ions, the share
  of the pairs lost to recombination agrees with Boag's formula to 1%, and the current
  extrapolated from two meshes agrees to 0.01%. A transient run ends at the
  steady-state current. These are tests (`tests/test_verification.py`).

Manufactured solution
: With densities chosen in advance and the sources worked out from them, the solved
  densities converge to the chosen ones at second order when the diffusion floor is
  inactive, and at first order when it sets the diffusion. Also tests.

Space charge
: The perturbation potential was verified against the analytic solution for a
  uniformly charged coaxial gap to 10⁻⁷.

Mesh convergence
: The error is first order in the cell size. Each quoted current is extrapolated from
  150 and 300 radial cells. The tests check the order on the parallel-plate and the
  coaxial chamber, and that the extrapolated coaxial current at trace level is the
  charge of all the pairs created.
