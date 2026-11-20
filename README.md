# ion_chambers

A [FESTIM](https://github.com/festim-dev/FESTIM) model of tritium ionisation chambers
(in-line process monitors). The goal is to reproduce existing measurements first, then
use the validated model to explore chamber design: the linear range, memory effects
and the detection threshold.

**Status:** prototype. The model runs for a generic coaxial chamber in room air. It
has not been compared with measurements yet, because the real vessel geometry and the
carrier gases of the experiments (N₂, Ar) are not in the model.

## Setup

```bash
conda env create -f environment.yml
conda activate ion-chambers
pip install -e .
```

This gives FESTIM 2.2 on DOLFINx 0.11, from conda-forge.

Quick check, about 5 s:

```bash
python scripts/quick_check.py
```

## Example

```python
import ion_chambers as ic

electrons = ic.ChargeCarrier(name="e", charge=-1, mobility=0.04, diffusivity=0.04)
positive_ions = ic.ChargeCarrier(name="P", charge=1, mobility=2.0e-4)
negative_ions = ic.ChargeCarrier(name="N", charge=-1, mobility=2.2e-4)

air = ic.Gas(
    name="air",
    w_value=33.97,  # eV per ion pair
    density=1.205,  # kg/m3
    carriers=[electrons, positive_ions, negative_ions],
    reactions=[
        ic.Reaction(reactants=["e"], products=["N"], rate=1.0e8),  # attachment
        ic.Reaction(reactants=["P", "N"], rate=2.0e-12),  # recombination
        ic.Reaction(reactants=["P", "e"], rate=2.0e-13),  # recombination
    ],
)

model = ic.IonChamberModel(
    geometry=ic.CoaxialChamber(rod_radius=1.5e-3, can_radius=0.05, length=0.127),
    gas=air,
    sources=ic.TritiumInGas(activity=3.7e10),  # Bq/m3
    bias=500.0,  # V
)
model.initialise()
reading = model.run()
print(reading.current, reading.collection_efficiency)
```

A model is composed of a geometry, a gas, one or more ionisation sources, a potential
solver and a bias. The package ships no gas data: you define the charge carriers and
their reactions, as above. `scripts/reference_case.py` is a fuller air example that
follows the pressure. The documentation explains each part.

## Documentation

The documentation is in `docs/` and covers the user guide, the theory, the results so
far and the roadmap. Build it locally with:

```bash
sphinx-build -b html docs docs/_build/html
```

then open `docs/_build/html/index.html`.

## Layout

```
src/ion_chambers/    the package
tests/               unit and integration tests (python -m pytest)
scripts/             parameter studies (write JSON to results/, not tracked)
data/experimental/   measurements digitised from the TFG 2026 slides (see its README)
docs/                documentation (Sphinx)
```

## Running the studies

```bash
python scripts/response_curve.py                 # current vs activity, reference chamber
python scripts/aspect_ratio_sweep.py 0.4 1.27 7  # L/D sweep at 1 L
python scripts/bias_sweep.py 250 500 1000        # saturation curve
```

Each solve takes seconds to minutes. Currents are extrapolated from 150 and 300 radial
cells.

## Known gaps before validation

- **Geometry:** the model uses a thin 1.5 mm rod in a 5 cm radius can. The real 1 L
  monitor appears to have a much larger anode, about 2.8 cm diameter and 9 cm long
  (slide 5 of the TFG deck). We need the real drawings.
- **Carrier gas:** the only worked example is air. Most of the measurements are in N₂
  or Ar, where electrons barely attach and stay free. The model accepts such gases,
  but we have not yet chosen transport data for them.
- **20 cc wire cage:** this is a different geometry (an open cage), not yet modelled.

See `docs/roadmap.md` for the plan.
