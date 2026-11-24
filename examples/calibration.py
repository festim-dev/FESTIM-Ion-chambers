"""Calibration: turn a response curve into a way of reading an activity off a current.

    python examples/response_curve.py   # first, writes results/response.json
    python examples/calibration.py

A monitor is used the other way round from the model: the current is measured and the
activity is wanted. This script takes the response curve of the model and works out

1. the sensitivity, the current per unit of activity where the response is linear;
2. up to which activity that linear calibration holds;
3. the activity that corresponds to a measured current, with and without the
   correction for recombination.

Edit the parameters below to change the case.
"""

import json

import matplotlib as mpl
import numpy as np

import ion_chambers as ic

mpl.use("Agg")  # the figure is saved to a file, not shown in a window
import matplotlib.pyplot as plt

INPUT_FILE = "results/response.json"
FIGURE_FILE = "results/calibration.png"
ALLOWED_DEVIATION = 0.05  # the response counts as linear while it is within 5%
SPACE_CHARGE_LIMIT = 0.1  # the model is indicative only above this
MEASURED_CURRENTS = [1e-12, 1e-9, 1e-8, 5e-8]  # A, readings to convert

CURIE = ic.constants.BQ_PER_CI  # Bq

# ---- 1. read the response curve --------------------------------------------------

with open(INPUT_FILE) as f:
    results = json.load(f)

activities = []
currents = []
space_charge_fields = []
for row in results["rows"]:
    activities.append(row["activity"])
    currents.append(row["current_extrapolated"])
    space_charge_fields.append(row["space_charge_field"])
activities = np.array(activities)  # Bq/m3
currents = np.array(currents)  # A
space_charge_fields = np.array(space_charge_fields)

# ---- 2. the sensitivity ----------------------------------------------------------
# At the lowest activity nothing recombines, so the current is proportional to the
# activity. Their ratio is the calibration factor of the chamber.

sensitivity = currents[0] / activities[0]  # A per Bq/m3

print(f"sensitivity  {sensitivity:.3e} A per Bq/m3")
print(f"             {sensitivity * CURIE * 1e-6 * 1e15:.3f} fA per uCi/m3")

# ---- 3. how far the linear calibration holds -------------------------------------
# The ratio of the current to the linear prediction falls from 1 as recombination
# sets in. The linear range ends where it crosses 1 - ALLOWED_DEVIATION. The crossing
# is interpolated between the two points around it, on a logarithmic activity axis.

linearity = currents / (sensitivity * activities)

# np.interp needs increasing values, and the linearity decreases with the activity
log_activity_at_limit = np.interp(
    1 - ALLOWED_DEVIATION, linearity[::-1], np.log10(activities)[::-1]
)
linear_limit = 10**log_activity_at_limit  # Bq/m3

log_activity_at_space_charge = np.interp(
    np.log10(SPACE_CHARGE_LIMIT), np.log10(space_charge_fields), np.log10(activities)
)
space_charge_limit = 10**log_activity_at_space_charge  # Bq/m3

print(
    f"response {ALLOWED_DEVIATION:.0%} below linear at  "
    f"{linear_limit:.2e} Bq/m3 ({linear_limit / CURIE:.1f} Ci/m3)"
)
print(
    f"space-charge field reaches {SPACE_CHARGE_LIMIT} at  "
    f"{space_charge_limit:.2e} Bq/m3 ({space_charge_limit / CURIE:.1f} Ci/m3)"
)

# ---- 4. from a measured current to an activity -----------------------------------
# The linear calibration divides the current by the sensitivity. Above the linear
# range it reads too low, because part of the charge recombined. Reading the
# activity off the full response curve corrects for that.

print()
print("current (A)   activity, linear (Ci/m3)   activity, full curve (Ci/m3)")
for measured_current in MEASURED_CURRENTS:
    linear_activity = measured_current / sensitivity
    log_corrected_activity = np.interp(
        np.log10(measured_current), np.log10(currents), np.log10(activities)
    )
    corrected_activity = 10**log_corrected_activity
    print(
        f"{measured_current:.1e}       {linear_activity / CURIE:<24.4g} "
        f"{corrected_activity / CURIE:.4g}"
    )

# ---- 5. the figure ---------------------------------------------------------------

figure, axes = plt.subplots(1, 2, figsize=(11, 4))

axes[0].loglog(activities / CURIE, currents, "o-", label="model")
axes[0].loglog(
    activities / CURIE, sensitivity * activities, "k--", label="linear calibration"
)
axes[0].set_ylabel("chamber current (A)")
axes[0].legend()

axes[1].semilogx(activities / CURIE, linearity, "o-")
axes[1].axhline(1 - ALLOWED_DEVIATION, color="k", linestyle="--")
axes[1].axvline(linear_limit / CURIE, color="k", linestyle=":")
axes[1].set_ylabel("current / linear calibration")

for axis in axes:
    axis.set_xlabel("activity (Ci/m³)")
figure.tight_layout()
figure.savefig(FIGURE_FILE, dpi=150)
print()
print(f"saved {FIGURE_FILE}")
