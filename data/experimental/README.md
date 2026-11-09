# Experimental data (digitised)

Experimental points digitised from the figures in

> N. J. Sovare, W. T. Shmayda, *Modelling the Performance of Process Monitors from
> First Principles*, 50th Tritium Focus Group, Aiken SC, 19 October 2026.

These are **read off raster plots**, not the original measurements. Replace them with
the raw data from Tritium Solutions when it becomes available. Only the experimental
markers were extracted. Model and theory curves were left out.

## Files

| CSV | Slide | Device / gas | Columns |
|---|---|---|---|
| `calibration_1L_N2_slide12.csv` | 12 | 1 L monitor, N₂ | flow (LPM), bias (V), concentration (µCi/m³), current (fA) |
| `calibration_1L_N2_slide20.csv` | 20 | 1 L monitor, N₂ | same, y = indicated activity (µCi/m³) |
| `voltage_dependence_slide12.csv` | 12 | 1 L monitor (assumed), gas not stated | series (9.7 mCi/m³, 105 mCi/m³, 1.74 Ci/m³), bias magnitude (V), detected activity (Ci/m³) |
| `Ar_pressure_dependence_slide12.csv` | 12 | 1 L monitor (assumed), Ar, 9.45 mCi/m³ | pressure (Torr), detected activity (Ci/m³) |
| `wire_cage_20cc_pressures_slide14.csv` | 14 | 20 cc wire cage, gas not stated | pressure (Torr), concentration (Ci/m³), current (A) |
| `wire_cage_20cc_Ar_slide21.csv` | 21 | 20 cc wire cage, Ar | concentration (Ci/m³), current (A) |

Every row has a `quality` flag:

- `ok`: the marker is fully or nearly fully visible.
- `occluded`: less than 60% of the marker is visible.
- `split`: the marker was separated from touching markers of the same colour.
- `clipped`: the marker is cut by the plot frame.

## Caveats

- **Slide 12 axis label.** The y-axis reads "Current (fA) = Detected Activity (µCi/m³)", so the calibration is 1 fA per µCi/m³.
- **Slides 12 and 20** show the same calibration data with different plotting, and agree within about 0.1 µCi/m³ on the points they share. In slide 20, a 7.5 LPM / −25 V point at 40 µCi/m³ covers the 5 LPM / −100 V point that slide 12 shows there. Treat slide 12 as the primary copy.
- **Crowded region, 2–10 µCi/m³.** Markers overlap heavily here, so some points are probably hidden entirely and the `split` positions are uncertain to about ±0.5 µCi/m³.
- **Voltage dependence.** The nominal voltages (1, 3, 10, …, 500 V) are inferred from the evenly placed markers. The digitised values are kept alongside.
- **Slide 14.** The ten pressure series (100–1000 Torr) sit almost exactly on top of each other, with 1000 Torr drawn last. Only markers with at least 20 visible pixels are kept: all of 1000 Torr, 100 Torr where it shows below the stack, and 200 and 300 Torr at 15 Ci/m³. The other pressures lie within one marker diameter (about ±10% in current) of the 1000 Torr points.
- **Slide 21.** The point at 54 Ci/m³ is almost entirely hidden under the model marker. The point at 15 Ci/m³ is cut by the bottom axis.
- **Accuracy.** About ±1 px against axes calibrated from the tick marks. That's roughly ±0.1 µCi/m³ on the linear plots and ±1–2% on the log plots, which is better than the marker size.

## Reproducing

```bash
conda activate ion-chambers
python data/experimental/digitise_slides.py
```

The source images are in `figures/`, unmodified from the deck. The script writes the CSVs plus `figures/checks/`, which overlays the detected markers on each figure at 2× scale (magenta = ok, orange = flagged). Look at these before relying on a point.
