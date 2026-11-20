"""Shape study: response curves of 1 L chambers with different length/diameter.

python scripts/aspect_ratio_sweep.py 0.4 1.27 7
python scripts/aspect_ratio_sweep.py 1.27 --output results/aspect_test.json
"""

import argparse
import time

import numpy as np

import ion_chambers as ic

ROD_RADIUS = 1.5e-3  # m
VOLUME = 1e-3  # m3


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("aspect_ratios", type=float, nargs="+", help="length/diameter")
    parser.add_argument("--bias", type=float, default=500.0, help="anode bias (V)")
    parser.add_argument("--output", default="results/aspect_ratio.json")
    args = parser.parse_args()

    activities = np.r_[3.7e7, np.logspace(8, 14, 13)]  # Bq/m3
    results = []
    for aspect_ratio in args.aspect_ratios:
        geometry = ic.CoaxialChamber.from_aspect_ratio(
            volume=VOLUME, rod_radius=ROD_RADIUS, aspect_ratio=aspect_ratio
        )

        def build_model(resolution: int, geometry=geometry) -> ic.IonChamberModel:
            return ic.IonChamberModel(
                geometry=geometry,
                gas=ic.library.air(),
                sources=ic.TritiumInGas(activity=0.0),
                bias=args.bias,
                resolution=resolution,
            )

        start = time.perf_counter()
        rows = ic.response_curve(build_model, activities, label=f"L/D={aspect_ratio}")
        results.append(
            {
                "aspect_ratio": aspect_ratio,
                "rod_radius": geometry.rod_radius,
                "can_radius": geometry.can_radius,
                "length": geometry.length,
                "bias": args.bias,
                "rows": rows,
            }
        )
        print(
            f"== L/D={aspect_ratio}: R={geometry.can_radius * 100:.2f} cm, "
            f"L={geometry.length * 100:.2f} cm, {time.perf_counter() - start:.0f} s",
            flush=True,
        )
        ic.save_json(results, args.output)  # keep what is done if a later case fails


if __name__ == "__main__":
    main()
