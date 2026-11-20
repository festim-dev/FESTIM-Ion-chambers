"""Quick end-to-end check: reference chamber on a coarse mesh, three activities.

python scripts/quick_check.py
"""

import time

from reference_case import air, reference_chamber

import ion_chambers as ic


def main():
    model = ic.IonChamberModel(
        geometry=reference_chamber(),
        gas=air(),
        sources=ic.TritiumInGas(activity=0.0),
        bias=500.0,
        resolution=100,
    )
    start = time.perf_counter()
    model.initialise()
    print(f"setup {time.perf_counter() - start:.1f} s")

    for activity in [3.7e4, 3.7e10, 3.7e13]:  # Bq/m3
        model.sources[0].activity = activity
        start = time.perf_counter()
        reading = model.run()
        print(f"\nactivity {activity:.2e} Bq/m3 ({time.perf_counter() - start:.1f} s)")
        for name, value in reading.as_dict().items():
            print(f"  {name:26s} {value:.4g}")


if __name__ == "__main__":
    main()
