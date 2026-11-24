# Developing

## Layout

```
src/ion_chambers/    the package
tests/               unit and integration tests
examples/            parameter studies, their post-processing, and the example air and chamber
docs/                this documentation
```

The package mirrors FESTIM: problem classes that are set up with a mesh, subdomains,
boundary conditions and exports, then initialised and run.

| Module | Content |
|---|---|
| `electrostatic_problem.py` | `ElectrostaticProblem`, `FixedPotentialBC`, `PotentialExport` |
| `charge_transport_problem.py` | `ChargeTransportProblem`, the stabilisation (`DiffusionFloor`), and three children of FESTIM classes: `MobilityElectromigration`, `AbsorbingElectrodeBC` and `CellwiseDiffusionMaterial` |
| `exports.py` | the chamber quantities, children of FESTIM's `DerivedQuantity` |
| `gas.py` | `Gas`, `ChargeCarrier` (a `festim.Species`), `CarrierReaction` (a `festim.GenericReaction`) |
| `geometry.py` | `Geometry`, `CoaxialChamber` |
| `source.py` | `IonisationSource`, `TritiumInGas` |
| `beta.py` | beta spectrum, range, wall-loss Monte Carlo |
| `constants.py` | physical constants |

## Tests

```bash
python -m pytest                      # everything, about a minute
python -m pytest -m "not slow"        # unit tests only
python -m pytest tests/test_gas.py -k pressure
```

Tests marked `slow` solve full problems. `tests/test_charge_transport_problem.py`
includes a regression against the currents of the original prototype.

## Style

```bash
ruff format src tests examples docs
ruff check src tests examples
```

The conventions follow FESTIM's: annotated signatures, Google-style docstrings without
types, input validated when an object is created with an error naming the attribute
and the type received, and small tests with one assertion each and a one-sentence
docstring.

Clarity comes before brevity: prefer plain loops and named intermediate values to
compact constructs.

## Documentation

```bash
sphinx-build -W -b html docs docs/_build/html
```

Pages are Markdown (MyST). The API reference is generated from the docstrings.
