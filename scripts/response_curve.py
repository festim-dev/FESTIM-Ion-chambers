"""Response curve: chamber current against tritium activity, reference 1 L chamber.

python scripts/response_curve.py
python scripts/response_curve.py --bias 250 --output results/response_250V.json
"""

import argparse

import numpy as np
from reference_case import air, reference_chamber

import ion_chambers as ic


def build_model(bias, resolution):
    """The reference chamber in air, on a mesh with the given resolution."""
    return ic.IonChamberModel(
        geometry=reference_chamber(),
        gas=air(),
        sources=ic.TritiumInGas(activity=0.0),
        bias=bias,
        resolution=resolution,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bias", type=float, default=500.0, help="anode bias (V)")
    parser.add_argument("--output", default="results/response.json")
    args = parser.parse_args()

    activities = np.r_[np.logspace(3, 8, 3), np.logspace(9, 14, 11)]  # Bq/m3

    rows = ic.response_curve(
        coarse_model=build_model(args.bias, resolution=150),
        fine_model=build_model(args.bias, resolution=300),
        activities=activities,
        label="reference",
    )

    geometry = reference_chamber()
    results = {
        "rod_radius": geometry.rod_radius,
        "can_radius": geometry.can_radius,
        "length": geometry.length,
        "bias": args.bias,
        "rows": rows,
    }
    ic.save_json(results, args.output)


if __name__ == "__main__":
    main()
