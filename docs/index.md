# FESTIM - Ion chambers

A companion package for [FESTIM](https://github.com/festim-dev/FESTIM) to model
tritium ionisation chambers, the in-line process monitors that read a tritium
concentration as an electric current.

Tritium betas ionise the gas in the chamber. A bias between two electrodes sweeps the
charges out before they recombine, and the current that flows is the reading. This
package solves for the charge densities in the chamber and predicts that current from
the geometry, the gas, the bias and the tritium activity.

The [theory](theory.md) page shows how a chamber works, one mechanism at a time, and
the [walkthrough](walkthrough.md) builds a complete model step by step.

It extends FESTIM's classes and is used like FESTIM: an `ElectrostaticProblem` for the electric
potential and a `ChargeTransportProblem` for the charge carriers, each with a mesh,
subdomains, boundary conditions and exports.

:::{admonition} Status: prototype
:class: warning
The model runs for a generic coaxial chamber in air. It has **not been validated
against measurements**. See the [roadmap](roadmap.md).
:::

## What it can do today

- Model a closed coaxial chamber in 2D axisymmetry, with any dimensions.
- Describe the gas as a set of charge carriers and reactions that you define, with
  coefficients that may depend on the reduced electric field E/N. No gas data is
  built in.
- Account for the beta energy lost to the walls with a Monte Carlo over the tritium
  spectrum.
- Predict the chamber current, the collection efficiency and the wall-loss efficiency,
  as FESTIM-style exports.
- Report how much the space charge would distort the field, so you know when the
  fixed-field assumption stops holding.
- Examples for response curves (current against activity), saturation curves
  (current against bias) and shape studies, with mesh extrapolation.

```{toctree}
:caption: Getting started
:maxdepth: 1

installation
quickstart
walkthrough
example_air
```

```{toctree}
:caption: User guide
:maxdepth: 1

guide/geometry
guide/electrostatics
guide/gas
guide/sources
guide/transport
guide/exports
guide/studies
```

```{toctree}
:caption: Background
:maxdepth: 1

theory
results
roadmap
```

```{toctree}
:caption: Reference
:maxdepth: 1

api
developing
```
