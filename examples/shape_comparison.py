"""Shape comparison: sensitivity and linear range of chambers of different shapes.

    python examples/aspect_ratio_sweep.py   # first, writes results/aspect_ratio.json
    python examples/shape_comparison.py

For each chamber of the shape sweep, this script works out the sensitivity and the
activity up to which the response is linear, as calibration.py does for one chamber.
The two pull in opposite directions, which is the trade-off in choosing a shape.

Edit the parameters below to change the case.
"""

import json

import matplotlib as mpl
import numpy as np

import ion_chambers as ic

mpl.use("Agg")  # the figure is saved to a file, not shown in a window
import matplotlib.pyplot as plt

INPUT_FILE = "results/aspect_ratio.json"
FIGURE_FILE = "results/shape_comparison.png"
ALLOWED_DEVIATION = 0.05  # the response counts as linear while it is within 5%

CURIE = ic.constants.BQ_PER_CI  # Bq

# the file holds one entry per chamber shape, each with one row per activity
with open(INPUT_FILE) as f:
    results = json.load(f)

plt.figure(figsize=(6, 4))

print("L/D     sensitivity (A per Bq/m3)   wall efficiency   linear up to (Ci/m3)")
for entry in results:
    activities = []
    currents = []
    for row in entry["rows"]:
        activities.append(row["activity"])
        currents.append(row["current_extrapolated"])
    activities = np.array(activities)  # Bq/m3
    currents = np.array(currents)  # A

    # at the lowest activity nothing recombines: the current follows the activity
    sensitivity = currents[0] / activities[0]
    wall_efficiency = entry["rows"][0]["wall_efficiency"]

    # The linear range ends where the current falls ALLOWED_DEVIATION below the
    # linear prediction, interpolated on a logarithmic activity axis. np.interp
    # needs increasing values, and the linearity decreases with the activity.
    linearity = currents / (sensitivity * activities)
    log_activity_at_limit = np.interp(
        1 - ALLOWED_DEVIATION, linearity[::-1], np.log10(activities)[::-1]
    )
    linear_limit = 10**log_activity_at_limit  # Bq/m3

    aspect_ratio = entry["aspect_ratio"]
    print(
        f"{aspect_ratio:<7g} {sensitivity:<27.4e} {wall_efficiency:<17.4f} "
        f"{linear_limit / CURIE:.3g}"
    )

    plt.semilogx(activities / CURIE, linearity, "o-", label=f"L/D = {aspect_ratio:g}")

plt.axhline(1 - ALLOWED_DEVIATION, color="k", linestyle="--")
plt.xlabel("activity (Ci/m³)")
plt.ylabel("current / linear calibration")
plt.legend()
plt.tight_layout()
plt.savefig(FIGURE_FILE, dpi=150)
print(f"saved {FIGURE_FILE}")
