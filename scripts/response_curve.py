"""Response curve: chamber current against tritium activity, reference 1 L chamber.

python scripts/response_curve.py
python scripts/response_curve.py --bias 250 --output results/response_250V.json
"""

import argparse

import numpy as np

import ion_chambers as ic


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bias", type=float, default=500.0, help="anode bias (V)")
    parser.add_argument("--output", default="results/response.json")
    args = parser.parse_args()

    geometry = ic.library.reference_chamber()

    def build_model(resolution: int) -> ic.IonChamberModel:
        return ic.IonChamberModel(
            geometry=geometry,
            gas=ic.library.air(),
            sources=ic.TritiumInGas(activity=0.0),
            bias=args.bias,
            resolution=resolution,
        )

    activities = np.r_[np.logspace(3, 8, 3), np.logspace(9, 14, 11)]  # Bq/m3
    rows = ic.response_curve(build_model, activities, label="reference")
    ic.save_json(
        {
            "rod_radius": geometry.rod_radius,
            "can_radius": geometry.can_radius,
            "length": geometry.length,
            "bias": args.bias,
            "rows": rows,
        },
        args.output,
    )


if __name__ == "__main__":
    main()
