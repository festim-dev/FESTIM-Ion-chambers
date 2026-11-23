"""Saturation curve: chamber current against bias, reference 1 L chamber.

    python scripts/bias_sweep.py

Every bias needs new problems, because the electric potential changes. Every point
is solved on a coarse and on a fine mesh and extrapolated to zero cell size, as in
response_curve.py.

Edit the parameters below to change the case.
"""

import json
import os

from reference_case import build_transport_problem, reference_chamber

import ion_chambers as ic

BIASES = [250.0, 500.0, 1000.0]  # V, potential of the anode rod
COARSE_RESOLUTION = 150  # mesh cells across the electrode gap
FINE_RESOLUTION = 300
OUTPUT_FILE = "results/bias.json"

CURIE = ic.constants.BQ_PER_CI  # Bq
ACTIVITIES = [1 * CURIE, 10 * CURIE, 100 * CURIE]  # Bq/m3

cell_size_ratio = FINE_RESOLUTION / COARSE_RESOLUTION
chamber = reference_chamber()

results = []
for bias in BIASES:
    print(f"bias {bias:g} V")

    # the problem on the coarse mesh, which only gives a current
    coarse_current = ic.ChamberCurrent(surface=chamber.anode)

    coarse_problem = build_transport_problem(chamber, bias, COARSE_RESOLUTION)
    coarse_problem.exports = [coarse_current]
    coarse_problem.initialise()
    coarse_tritium = coarse_problem.ionisation_sources[0]

    # the problem on the fine mesh, which gives every quantity
    fine_current = ic.ChamberCurrent(surface=chamber.anode)
    ideal_current = ic.IdealCurrent()
    wall_efficiency = ic.WallEfficiency()
    collection_efficiency = ic.CollectionEfficiency()
    space_charge_field = ic.SpaceChargeField()

    fine_problem = build_transport_problem(chamber, bias, FINE_RESOLUTION)
    fine_problem.exports = [
        fine_current,
        ideal_current,
        wall_efficiency,
        collection_efficiency,
        space_charge_field,
    ]
    fine_problem.initialise()
    fine_tritium = fine_problem.ionisation_sources[0]

    # solve both problems at every activity
    rows = []
    for activity in ACTIVITIES:
        coarse_tritium.activity = activity
        coarse_problem.run()

        fine_tritium.activity = activity
        fine_problem.run()

        # Richardson extrapolation of a first-order error to zero cell size
        difference = fine_current.value - coarse_current.value
        extrapolated_current = fine_current.value + difference / (cell_size_ratio - 1)

        print(f"  activity {activity:.2e} Bq/m3   current {extrapolated_current:.4e} A")

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

    results.append({"bias": bias, "rows": rows})
    # saved after every bias, so a failure later on does not lose the work done
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        json.dump(results, f, indent=1)

print(f"saved {OUTPUT_FILE}")
