"""Response curve: chamber current against tritium activity, reference 1 L chamber.

    python scripts/response_curve.py

The current computed on one mesh has an error proportional to the cell size, so
every point is solved on a coarse and on a fine mesh and extrapolated to zero cell
size.

Edit the parameters below to change the case.
"""

import json
import os

import numpy as np
from reference_case import build_transport_problem, reference_chamber

import ion_chambers as ic

BIAS = 500.0  # V, potential of the anode rod
COARSE_RESOLUTION = 150  # mesh cells across the electrode gap
FINE_RESOLUTION = 300
OUTPUT_FILE = "results/response.json"

# activity concentrations (Bq/m3): three points at trace level, then 1e9 to 1e14
trace_activities = np.logspace(3, 8, 3)
high_activities = np.logspace(9, 14, 11)
activities = np.concatenate([trace_activities, high_activities])

chamber = reference_chamber()

# ---- 1. the problem on the coarse mesh, which only gives a current ---------------

coarse_current = ic.ChamberCurrent(surface=chamber.anode)

coarse_problem = build_transport_problem(chamber, BIAS, COARSE_RESOLUTION)
coarse_problem.exports = [coarse_current]
coarse_problem.initialise()
coarse_tritium = coarse_problem.ionisation_sources[0]

# ---- 2. the problem on the fine mesh, which gives every quantity -----------------

fine_current = ic.ChamberCurrent(surface=chamber.anode)
ideal_current = ic.IdealCurrent()
wall_efficiency = ic.WallEfficiency()
collection_efficiency = ic.CollectionEfficiency()
space_charge_field = ic.SpaceChargeField()

fine_problem = build_transport_problem(chamber, BIAS, FINE_RESOLUTION)
fine_problem.exports = [
    fine_current,
    ideal_current,
    wall_efficiency,
    collection_efficiency,
    space_charge_field,
]
fine_problem.initialise()
fine_tritium = fine_problem.ionisation_sources[0]

# ---- 3. solve both problems at every activity ------------------------------------
# a problem is initialised once, then run again after its activity is changed

cell_size_ratio = FINE_RESOLUTION / COARSE_RESOLUTION

rows = []
for activity in activities:
    coarse_tritium.activity = activity
    coarse_problem.run()

    fine_tritium.activity = activity
    fine_problem.run()

    # Richardson extrapolation of a first-order error to zero cell size
    difference = fine_current.value - coarse_current.value
    extrapolated_current = fine_current.value + difference / (cell_size_ratio - 1)

    print(f"activity {activity:.2e} Bq/m3   current {extrapolated_current:.4e} A")

    row = {
        "activity": float(activity),
        "current_coarse": float(coarse_current.value),
        "current_fine": float(fine_current.value),
        "current_extrapolated": float(extrapolated_current),
        "ideal_current": float(ideal_current.value),
        "wall_efficiency": float(wall_efficiency.value),
        "collection_efficiency": float(collection_efficiency.value),
        "space_charge_field": float(space_charge_field.value),
    }
    rows.append(row)

# ---- 4. save the results ---------------------------------------------------------

results = {
    "rod_radius": chamber.rod_radius,
    "can_radius": chamber.can_radius,
    "length": chamber.length,
    "bias": BIAS,
    "rows": rows,
}

os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
with open(OUTPUT_FILE, "w") as f:
    json.dump(results, f, indent=1)
print(f"saved {OUTPUT_FILE}")
