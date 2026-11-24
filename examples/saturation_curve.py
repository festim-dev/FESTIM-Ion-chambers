"""Saturation curve: analyse the current against bias as a measured one would be.

    python examples/bias_sweep.py   # first, writes results/bias.json
    python examples/saturation_curve.py

Raising the bias collects the charges faster, so fewer recombine and the current
rises towards its saturation value, the current with no recombination at all. This
script takes the bias sweep of the model and works out, for each activity,

1. the bias needed to collect a given share of the charge;
2. the saturation current, by the extrapolation used on measured curves: a straight
   line through 1/I against 1/V^2, read at 1/V^2 = 0. The model also knows the true
   saturation current, so this tells how well the extrapolation does.

Edit the parameters below to change the case.
"""

import json

import matplotlib as mpl
import numpy as np

import ion_chambers as ic

mpl.use("Agg")  # the figure is saved to a file, not shown in a window
import matplotlib.pyplot as plt

INPUT_FILE = "results/bias.json"
FIGURE_FILE = "results/saturation_curve.png"
TARGET_EFFICIENCY = 0.99  # the share of the charge to collect
FIT_FROM_BIAS = 500.0  # V, the straight line is fitted to the points from here up

CURIE = ic.constants.BQ_PER_CI  # Bq

# ---- 1. read the bias sweep ------------------------------------------------------
# the file holds one entry per bias, each with one row per activity

with open(INPUT_FILE) as f:
    results = json.load(f)

biases = []
for entry in results:
    biases.append(entry["bias"])
biases = np.array(biases)  # V

activities = []
for row in results[0]["rows"]:
    activities.append(row["activity"])

figure, axes = plt.subplots(1, 2, figsize=(11, 4))

for i, activity in enumerate(activities):
    currents = []
    efficiencies = []
    for entry in results:
        row = entry["rows"][i]
        currents.append(row["current_extrapolated"])
        efficiencies.append(row["collection_efficiency"])
    currents = np.array(currents)  # A
    efficiencies = np.array(efficiencies)

    # the current if nothing recombined: every pair created in the gas is collected
    first_row = results[0]["rows"][i]
    true_saturation = first_row["ideal_current"] * first_row["wall_efficiency"]

    print(f"activity {activity / CURIE:g} Ci/m3")

    # ---- 2. the bias for the target collection efficiency ------------------------
    # interpolated between the two biases around it, on a logarithmic bias axis

    if efficiencies[-1] < TARGET_EFFICIENCY:
        print(f"  {TARGET_EFFICIENCY:.0%} collection needs more than {biases[-1]:g} V")
    elif efficiencies[0] > TARGET_EFFICIENCY:
        print(f"  {TARGET_EFFICIENCY:.0%} collection needs less than {biases[0]:g} V")
    else:
        log_bias = np.interp(TARGET_EFFICIENCY, efficiencies, np.log10(biases))
        print(f"  {TARGET_EFFICIENCY:.0%} collection at {10**log_bias:.0f} V")

    # ---- 3. the saturation current from 1/I against 1/V^2 ------------------------
    # Near saturation the loss to recombination goes as 1/V^2, so 1/I is a straight
    # line in 1/V^2 and its value at 1/V^2 = 0 is 1/I_saturation.

    in_fit = biases >= FIT_FROM_BIAS
    slope, intercept = np.polyfit(1 / biases[in_fit] ** 2, 1 / currents[in_fit], 1)
    fitted_saturation = 1 / intercept
    error = fitted_saturation / true_saturation - 1

    print(f"  saturation current, fitted  {fitted_saturation:.4e} A")
    print(f"  saturation current, true    {true_saturation:.4e} A")
    print(f"  error of the fit            {error:+.2%}")

    # ---- 4. the figure -----------------------------------------------------------

    label = f"{activity / CURIE:g} Ci/m³"
    axes[0].semilogx(biases, efficiencies, "o-", label=label)

    # divided by the true value so that the three activities share one axis
    line = axes[1].plot(1 / biases**2, true_saturation / currents, "o", label=label)
    fit_axis = np.linspace(0, 1 / FIT_FROM_BIAS**2, 20)
    axes[1].plot(
        fit_axis,
        true_saturation * (intercept + slope * fit_axis),
        color=line[0].get_color(),
    )

axes[0].axhline(TARGET_EFFICIENCY, color="k", linestyle="--")
axes[0].set_xlabel("bias (V)")
axes[0].set_ylabel("collection efficiency")
axes[0].legend()

axes[1].set_xlabel("1 / bias² (1/V²)")
axes[1].set_ylabel("saturation current / current")
axes[1].legend()

figure.tight_layout()
figure.savefig(FIGURE_FILE, dpi=150)
print(f"saved {FIGURE_FILE}")
