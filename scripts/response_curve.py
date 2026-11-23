"""Response curve: chamber current against tritium activity, reference 1 L chamber.

    python scripts/response_curve.py

Edit the parameters below to change the case.
"""

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
coarse_problem = build_transport_problem(chamber, BIAS, COARSE_RESOLUTION)
fine_problem = build_transport_problem(chamber, BIAS, FINE_RESOLUTION)

rows = ic.response_curve(
    coarse_problem=coarse_problem,
    fine_problem=fine_problem,
    activities=activities,
    cell_size_ratio=FINE_RESOLUTION / COARSE_RESOLUTION,
    label="reference",
)

results = {
    "rod_radius": chamber.rod_radius,
    "can_radius": chamber.can_radius,
    "length": chamber.length,
    "bias": BIAS,
    "rows": rows,
}
ic.save_json(results, OUTPUT_FILE)
print(f"saved {OUTPUT_FILE}")
