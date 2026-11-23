"""Shape study: response curves of 1 L chambers with different length/diameter.

    python scripts/aspect_ratio_sweep.py

Edit the parameters below to change the case.
"""

import numpy as np
from reference_case import build_transport_problem

import ion_chambers as ic

ASPECT_RATIOS = [0.4, 1.27, 7.0]  # length / diameter of the can
VOLUME = 1e-3  # m3
ROD_RADIUS = 1.5e-3  # m
BIAS = 500.0  # V, potential of the anode rod
COARSE_RESOLUTION = 150  # mesh cells across the electrode gap
FINE_RESOLUTION = 300
OUTPUT_FILE = "results/aspect_ratio.json"

# activity concentrations (Bq/m3): 1 mCi/m3, then 1e8 to 1e14
activities = np.concatenate([[3.7e7], np.logspace(8, 14, 13)])

results = []
for aspect_ratio in ASPECT_RATIOS:
    chamber = ic.CoaxialChamber.from_aspect_ratio(
        volume=VOLUME, rod_radius=ROD_RADIUS, aspect_ratio=aspect_ratio
    )
    print(
        f"L/D = {aspect_ratio}: can radius {chamber.can_radius * 100:.2f} cm, "
        f"length {chamber.length * 100:.2f} cm"
    )

    coarse_problem = build_transport_problem(chamber, BIAS, COARSE_RESOLUTION)
    fine_problem = build_transport_problem(chamber, BIAS, FINE_RESOLUTION)

    rows = ic.response_curve(
        coarse_problem=coarse_problem,
        fine_problem=fine_problem,
        activities=activities,
        cell_size_ratio=FINE_RESOLUTION / COARSE_RESOLUTION,
        label=f"L/D={aspect_ratio}",
    )
    results.append(
        {
            "aspect_ratio": aspect_ratio,
            "rod_radius": chamber.rod_radius,
            "can_radius": chamber.can_radius,
            "length": chamber.length,
            "bias": BIAS,
            "rows": rows,
        }
    )

    # saved after every shape, so a failure later on does not lose the work done
    ic.save_json(results, OUTPUT_FILE)

print(f"saved {OUTPUT_FILE}")
