# Example: air and the reference chamber

The package does not ship any gas data. You describe the gas yourself, so that every
number in a result is one you chose and can check.

This page shows a complete example: room air, the 1 L chamber used for the
[prototype results](results.md), and a function that sets up the two problems for
them. It is the file `examples/reference_case.py`, which the scripts in that folder
import. Copy it as a starting point for your own gas.

:::{warning}
The values below are reasonable starting points, not validated reference data. The
way they are scaled with pressure and temperature is deliberately simple. Check both
against your conditions before relying on a result.
:::

```{literalinclude} ../examples/reference_case.py
:language: python
```

## Values at 1 atm and 293.15 K

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

## How the example scales with pressure and temperature

With `n` the gas density relative to 1 atm and 293.15 K:

- the mass density scales as `n`
- the mobilities scale as `1 / n`
- the attachment frequency scales as `n²`, because attachment to O₂ is a three-body
  process
- the recombination coefficients and the electron characteristic energy are kept
  constant, which does not hold over a wide pressure range

## Using it

From a script in the `examples/` folder:

```python
from reference_case import air, reference_chamber

import ion_chambers as ic

chamber = reference_chamber()
gas = air(pressure=50_000.0)  # Pa
```
