# ion_chambers

A FESTIM model of tritium ionisation chambers (in-line process monitors). The goal is
to reproduce existing measurements first, then use the validated model to explore
chamber design: the linear range, memory effects and the detection threshold.

**Status:** prototype. The model runs for a generic coaxial chamber in room air. It
has not been compared with measurements yet, because we don't have the real vessel
geometry or gas conditions in the model.

## Setup

```bash
conda env create -f environment.yml
conda activate ion-chambers
pip install -e .
```

This gives FESTIM 2.2rc2 (drift terms) on DOLFINx 0.11, from conda-forge.

Quick check, about 5 s:

```bash
python -m ion_chambers.chamber_model
```

## Layout

```
ion_chambers/        model package
  chamber_model.py   Chamber: geometry, mesh, field, species, FESTIM problem, readout
  beta_source.py     tritium beta spectrum and Monte Carlo wall-loss correction
scripts/             parameter sweeps (write JSON to results/, not tracked)
  sweeps.py          response curve and L/D sweep at fixed 1 L volume
  bias_sweep.py      saturation curve (current vs bias)
data/experimental/   measurements digitised from the TFG 2026 slides (see its README)
docs/
  prototype_report.html   results of the prototype model
```

## The model

The model is 2D axisymmetric (r, z). An anode rod sits on the axis, the cylindrical
can and end caps are the cathode, and insulating rings separate them where the rod
passes through the caps.

- **Species:** free electrons, positive ions and negative ions (O₂⁻), as number densities.
- **Transport:** Nernst–Planck, `J = -D∇c + c v` with `v = -zμ∇φ`, using FESTIM's `ElectromigrationTerm`. The field is the fixed Laplace field. Space charge isn't fed back into the field; instead, the distortion it would cause is computed afterwards and reported.
- **Source:** `S(r,z) = A · (Ē/W) · f_wall(r,z)`. Here `f_wall` is the fraction of beta energy deposited in the gas rather than the walls, from a straight-track Monte Carlo over the tritium spectrum. It comes to about 168 ion pairs per decay in air.
- **Reactions:** electron attachment to O₂, ion–ion recombination and electron–ion recombination.
- **Boundaries:** electrodes absorb what the field drives into them and inject nothing (`OutflowBC`, plus `c = 0` for the species pushed away).
- **Reading:** the anode current from the Shockley–Ramo theorem. It agrees with the integrated electrode flux to better than 0.01%.
- **Numerics:** a cell-local diffusion floor stabilises drift-dominated transport. Currents are Richardson-extrapolated from two meshes.

The module docstring in `chamber_model.py` and the report in `docs/` give the details.

### Running the sweeps

```bash
python scripts/sweeps.py response              # current vs activity, reference chamber
python scripts/sweeps.py aspect 0.4,1.27,7 run1  # L/D sweep at 1 L
python scripts/bias_sweep.py 250,500,1000 run1   # saturation curve
```

Each solve takes seconds to minutes. The extrapolation uses N_r = 150 and 300 radial cells.

## Known gaps before validation

- **Geometry:** the model uses a thin 1.5 mm rod in a 5 cm radius can. The real 1 L monitor appears to have a much larger anode, about 2.8 cm diameter and 9 cm long (slide 5 of the TFG deck). We need the real drawings.
- **Carrier gas:** the parameters are hard-wired to air. Most of the measurements are in N₂ or Ar, where electrons barely attach and stay free, which changes the physics. W, mobilities and attachment need to be per-gas.
- **Pressure:** everything is at 1 atm, but some data sets vary the pressure.
- **20 cc wire cage:** this is a different geometry (an open cage), not yet modelled.

## Possible extensions

- Couple space charge back into the field (needed above about 50 Ci/m³).
- Wall uptake and release of tritium (memory effect and background drift).
- Gas flow (advection) and transient response to a concentration step.
- Columnar recombination along beta tracks, and scattering in the wall-loss Monte Carlo.
