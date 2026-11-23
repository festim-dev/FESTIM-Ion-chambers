"""Saturation curve: chamber current against bias, reference 1 L chamber.

    python scripts/bias_sweep.py

Edit the parameters below to change the case.
"""

from reference_case import build_transport_problem, reference_chamber

import ion_chambers as ic

BIASES = [250.0, 500.0, 1000.0]  # V, potential of the anode rod
COARSE_RESOLUTION = 150  # mesh cells across the electrode gap
FINE_RESOLUTION = 300
OUTPUT_FILE = "results/bias.json"

CURIE = ic.constants.BQ_PER_CI  # Bq
ACTIVITIES = [1 * CURIE, 10 * CURIE, 100 * CURIE]  # Bq/m3

chamber = reference_chamber()

results = []
for bias in BIASES:
    coarse_problem = build_transport_problem(chamber, bias, COARSE_RESOLUTION)
    fine_problem = build_transport_problem(chamber, bias, FINE_RESOLUTION)

    rows = ic.response_curve(
        coarse_problem=coarse_problem,
        fine_problem=fine_problem,
        activities=ACTIVITIES,
        cell_size_ratio=FINE_RESOLUTION / COARSE_RESOLUTION,
        label=f"{bias:g} V",
    )
    results.append({"bias": bias, "rows": rows})

    # saved after every bias, so a failure later on does not lose the work done
    ic.save_json(results, OUTPUT_FILE)

print(f"saved {OUTPUT_FILE}")
