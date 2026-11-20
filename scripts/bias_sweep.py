"""Saturation curve: chamber current against bias, reference 1 L chamber.

python scripts/bias_sweep.py 250 500 1000
python scripts/bias_sweep.py 500 --output results/bias_test.json
"""

import argparse
import time

from reference_case import air, reference_chamber

import ion_chambers as ic

BQ_PER_CI = ic.constants.BQ_PER_CI
ACTIVITIES = [1 * BQ_PER_CI, 10 * BQ_PER_CI, 100 * BQ_PER_CI]  # Bq/m3


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
    parser.add_argument("biases", type=float, nargs="+", help="anode bias (V)")
    parser.add_argument("--output", default="results/bias.json")
    args = parser.parse_args()

    results = []
    for bias in args.biases:
        start = time.perf_counter()
        rows = ic.response_curve(
            coarse_model=build_model(bias, resolution=150),
            fine_model=build_model(bias, resolution=300),
            activities=ACTIVITIES,
            label=f"{bias:g} V",
        )
        results.append({"bias": bias, "rows": rows})
        elapsed = time.perf_counter() - start
        print(f"== {bias:g} V done in {elapsed:.0f} s", flush=True)
        # saved after every bias, so a failure later on does not lose the work done
        ic.save_json(results, args.output)


if __name__ == "__main__":
    main()
