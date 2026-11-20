# Developing

## Layout

```
src/ion_chambers/    the package
tests/               unit and integration tests
scripts/             parameter studies (write JSON to results/, not tracked)
data/experimental/   digitised measurements
docs/                this documentation
```

Inside the package, each module holds one kind of object:

| Module | Content |
|---|---|
| `geometry.py` | `Geometry`, `CoaxialChamber` |
| `gas.py` | `Gas`, `ChargeCarrier`, `Reaction` |
| `library.py` | example gases and chambers |
| `beta.py`, `wall_loss.py` | beta spectrum, range, wall-loss Monte Carlo |
| `source.py` | `IonisationSource`, `TritiumInGas` |
| `field.py` | potential solvers |
| `stabilisation.py` | `DiffusionFloor` |
| `model.py` | `IonChamberModel` |
| `readout.py` | `Reading` and the functions that compute it |
| `studies.py` | mesh extrapolation and response curves |
| `mobility_electromigration.py`, `cellwise_diffusion_material.py` | subclasses of FESTIM classes |

## Tests

```bash
python -m pytest                      # everything, about 10 s
python -m pytest -m "not slow"        # unit tests only
python -m pytest tests/test_gas.py -k pressure
```

Tests marked `slow` solve a full model. `tests/test_model.py` includes a regression
against the currents of the original prototype.

## Style

```bash
ruff format .
ruff check src tests scripts
```

The conventions follow FESTIM's: annotated signatures, Google-style docstrings without
types, input validated when an object is created with an error naming the attribute
and the type received, and small tests with one assertion each.

## Documentation

```bash
sphinx-build -b html docs docs/_build/html
```

Pages are Markdown (MyST). The API reference is generated from the docstrings.
