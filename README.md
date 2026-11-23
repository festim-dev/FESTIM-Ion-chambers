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

The package is used like FESTIM: an `ElectrostaticProblem` gives the electric
potential, then a `ChargeTransportProblem` moves the charge carriers in it.

```python
import ion_chambers as ic

chamber = ic.CoaxialChamber(rod_radius=1.5e-3, can_radius=0.05, length=0.127)

# 1. the electric potential
electrostatics = ic.ElectrostaticProblem()
electrostatics.mesh = chamber.create_mesh(radial_cells=120, axial_cells=72)
electrostatics.subdomains = [chamber.gas, chamber.anode, chamber.cathode]
electrostatics.boundary_conditions = [
    ic.FixedPotentialBC(subdomain=chamber.anode, value=500.0),  # V
    ic.FixedPotentialBC(subdomain=chamber.cathode, value=0.0),
]
electrostatics.initialise()
electrostatics.run()

# 2. the gas: you define the charge carriers and their reactions
electrons = ic.ChargeCarrier(name="e", charge=-1, mobility=0.04, diffusivity=0.04)
positive_ions = ic.ChargeCarrier(name="P", charge=1, mobility=2.0e-4)
negative_ions = ic.ChargeCarrier(name="N", charge=-1, mobility=2.2e-4)

air = ic.Gas(
    name="air",
    w_value=33.97,  # eV per ion pair
    density=1.205,  # kg/m3
    carriers=[electrons, positive_ions, negative_ions],
    ionisation_products=(electrons, positive_ions),
    reactions=[
        # attachment of electrons to O2
        ic.CarrierReaction(reactant=electrons, product=negative_ions, rate=1.0e8),
        # recombination
        ic.CarrierReaction(reactant=[positive_ions, negative_ions], rate=2.0e-12),
        ic.CarrierReaction(reactant=[positive_ions, electrons], rate=2.0e-13),
    ],
)

# 3. the transport of the charge carriers
current = ic.ChamberCurrent(surface=chamber.anode)
collection_efficiency = ic.CollectionEfficiency()

transport = ic.ChargeTransportProblem()
transport.electrostatic_problem = electrostatics
transport.gas = air
transport.ionisation_sources = [
    ic.TritiumInGas(activity=3.7e10, geometry=chamber)  # Bq/m3
]
transport.exports = [current, collection_efficiency]
transport.initialise()
transport.run()

print(current.value, collection_efficiency.value)
```

The classes are children of FESTIM's: the problems, `ChargeCarrier` (a species),
`CarrierReaction` (a reaction), `FixedPotentialBC` and the exports. The package ships
no gas data: `scripts/reference_case.py` is a fuller air example that follows the
pressure.

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
python scripts/response_curve.py       # current vs activity, reference chamber
python scripts/aspect_ratio_sweep.py   # L/D sweep at 1 L
python scripts/bias_sweep.py           # saturation curve
```

Each script is a plain sequence of steps with its parameters at the top. Each solve
takes seconds to minutes. Currents are extrapolated from 150 and 300 radial cells.

## Known gaps before validation

- **Geometry:** the model uses a thin 1.5 mm rod in a 5 cm radius can. The real 1 L
  monitor appears to have a much larger anode, about 2.8 cm diameter and 9 cm long
  (slide 5 of the TFG deck). We need the real drawings.
- **Carrier gas:** the only worked example is air. Most of the measurements are in N₂
  or Ar, where electrons barely attach and stay free. The model accepts such gases,
  but we have not yet chosen transport data for them.
- **20 cc wire cage:** this is a different geometry (an open cage), not yet modelled.

See `docs/roadmap.md` for the plan.
