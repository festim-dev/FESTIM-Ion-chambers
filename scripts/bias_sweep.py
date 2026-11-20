"""Saturation curve: chamber current against bias, reference 1 L chamber.

python scripts/bias_sweep.py 250 500 1000
python scripts/bias_sweep.py 500 --output results/bias_test.json
"""

import argparse
import time

import ion_chambers as ic
from ion_chambers.constants import BQ_PER_CI

ACTIVITIES = [1 * BQ_PER_CI, 10 * BQ_PER_CI, 100 * BQ_PER_CI]  # Bq/m3


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("biases", type=float, nargs="+", help="anode bias (V)")
    parser.add_argument("--output", default="results/bias.json")
    args = parser.parse_args()

    geometry = ic.library.reference_chamber()
    results = []
    for bias in args.biases:

        def build_model(resolution: int, bias=bias) -> ic.IonChamberModel:
            return ic.IonChamberModel(
                geometry=geometry,
                gas=ic.library.air(),
                sources=ic.TritiumInGas(activity=0.0),
                bias=bias,
                resolution=resolution,
            )

        start = time.perf_counter()
        rows = ic.response_curve(build_model, ACTIVITIES, label=f"{bias:g} V")
        results.append({"bias": bias, "rows": rows})
        print(f"== {bias:g} V done in {time.perf_counter() - start:.0f} s", flush=True)
        ic.save_json(results, args.output)  # keep what is done if a later case fails


if __name__ == "__main__":
    main()
