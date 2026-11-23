"""Shape study: response curves of 1 L chambers with different length/diameter.

    python scripts/aspect_ratio_sweep.py

Every shape needs new problems, because the mesh changes. Every point is solved on a
coarse and on a fine mesh and extrapolated to zero cell size, as in
response_curve.py.

Edit the parameters below to change the case.
"""

import json
import os

import numpy as np
from reference_case import build_transport_problem

import ion_chambers as ic

ASPECT_RATIOS = [0.4, 1.27, 7.0]  # length / diameter of the can
# mesh cells from one end cap to the other on the coarse mesh, one per aspect ratio:
# the long chamber needs more
COARSE_AXIAL_CELLS = [90, 90, 150]
COARSE_RADIAL_CELLS = 150  # mesh cells from the rod to the can
REFINEMENT = 2  # the fine mesh has this many times more cells in both directions
VOLUME = 1e-3  # m3
ROD_RADIUS = 1.5e-3  # m
BIAS = 500.0  # V, potential of the anode rod
OUTPUT_FILE = "results/aspect_ratio.json"

# activity concentrations (Bq/m3): 1 mCi/m3, then 1e8 to 1e14
ACTIVITIES = np.concatenate([[3.7e7], np.logspace(8, 14, 13)])


results = []
for aspect_ratio, coarse_axial_cells in zip(
    ASPECT_RATIOS, COARSE_AXIAL_CELLS, strict=True
):
    chamber = ic.CoaxialChamber.from_aspect_ratio(
        volume=VOLUME, rod_radius=ROD_RADIUS, aspect_ratio=aspect_ratio
    )
    print(
        f"L/D = {aspect_ratio}: can radius {chamber.can_radius * 100:.2f} cm, "
        f"length {chamber.length * 100:.2f} cm"
    )

    # the problem on the coarse mesh, which only gives a current
    coarse_current = ic.ChamberCurrent(surface=chamber.anode)

    coarse_problem = build_transport_problem(
        chamber,
        BIAS,
        radial_cells=COARSE_RADIAL_CELLS,
        axial_cells=coarse_axial_cells,
    )
    coarse_problem.exports = [coarse_current]
    coarse_problem.initialise()
    coarse_tritium = coarse_problem.ionisation_sources[0]

    # the problem on the fine mesh, which gives every quantity
    fine_current = ic.ChamberCurrent(surface=chamber.anode)
    ideal_current = ic.IdealCurrent()
    wall_efficiency = ic.WallEfficiency()
    collection_efficiency = ic.CollectionEfficiency()
    space_charge_field = ic.SpaceChargeField()

    fine_problem = build_transport_problem(
        chamber,
        BIAS,
        radial_cells=REFINEMENT * COARSE_RADIAL_CELLS,
        axial_cells=REFINEMENT * coarse_axial_cells,
    )
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
        extrapolated_current = fine_current.value + difference / (REFINEMENT - 1)

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

    result = {
        "aspect_ratio": aspect_ratio,
        "rod_radius": chamber.rod_radius,
        "can_radius": chamber.can_radius,
        "length": chamber.length,
        "bias": BIAS,
        "rows": rows,
    }
    results.append(result)
    # saved after every shape, so a failure later on does not lose the work done
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    with open(OUTPUT_FILE, "w") as f:
        json.dump(results, f, indent=1)

print(f"saved {OUTPUT_FILE}")
