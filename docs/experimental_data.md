# Experimental data

The folder `data/experimental/` holds the measurements the model will be validated
against. They were digitised from the figures of

> N. J. Sovare, W. T. Shmayda, *Modelling the Performance of Process Monitors from
> First Principles*, 50th Tritium Focus Group, Aiken SC, 19 October 2026.

:::{warning}
These points were read off raster plots. They are not the original measurements and
should be replaced by the raw data when it becomes available.
:::

| File | Device and gas | Content |
|---|---|---|
| `calibration_1L_N2_slide12.csv` | 1 L monitor, N₂ | current against concentration, for several flows and biases |
| `calibration_1L_N2_slide20.csv` | 1 L monitor, N₂ | the same calibration, plotted differently |
| `voltage_dependence_slide12.csv` | 1 L monitor (assumed) | detected activity against bias at three concentrations |
| `Ar_pressure_dependence_slide12.csv` | 1 L monitor (assumed), Ar | detected activity against pressure at 9.45 mCi/m³ |
| `wire_cage_20cc_pressures_slide14.csv` | 20 cc wire cage | current against concentration at ten pressures |
| `wire_cage_20cc_Ar_slide21.csv` | 20 cc wire cage, Ar | current against concentration |

Every row carries a `quality` flag (`ok`, `occluded`, `split` or `clipped`) that says
how well the marker could be read. `data/experimental/README.md` lists the caveats of
each figure: read it before comparing the model with a data set.

## Why there is no comparison yet

None of these data sets can be compared with the model as it stands:

- The measurements are in N₂ or Ar. The model's only gas is air, where electrons
  attach and the physics differs.
- The real 1 L monitor has a much larger anode than the reference chamber, and its
  drawings are not in hand.
- Some data sets vary the pressure, which needs gas data that scale properly with it.
- The 20 cc wire cage is an open cage, not a coaxial can.

Closing these gaps is the first part of the [roadmap](roadmap.md).
