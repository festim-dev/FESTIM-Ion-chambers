"""Helpers for parameter studies: mesh extrapolation, response curves, saving."""

import json
from collections.abc import Callable, Sequence
from pathlib import Path

from ion_chambers.model import IonChamberModel


def richardson(coarse: float, fine: float, order: int = 1, ratio: float = 2.0) -> float:
    """Extrapolates a quantity computed on two meshes to zero cell size.

    Args:
        coarse: the value on the coarse mesh
        fine: the value on the fine mesh
        order: the order of the discretisation error. The diffusion floor makes the
            currents first order.
        ratio: the ratio of the coarse to the fine cell size

    Returns:
        the extrapolated value
    """
    return fine + (fine - coarse) / (ratio**order - 1)


def response_curve(
    build_model: Callable[[int], IonChamberModel],
    activities: Sequence[float],
    resolutions: tuple[int, int] = (150, 300),
    label: str = "",
) -> list[dict[str, float]]:
    """The chamber current against tritium activity, extrapolated from two meshes.

    For each resolution the model is built and initialised once, then solved for every
    activity by changing the activity of its first source.

    Args:
        build_model: called with a resolution, returns a model that is not yet
            initialised
        activities: the activity concentrations (Bq/m3)
        resolutions: the coarse and the fine mesh resolution
        label: a prefix for the progress lines

    Returns:
        one row per activity: the fine-mesh reading as a dictionary, plus ``activity``,
        ``current_coarse``, ``current_fine`` and ``current_extrapolated``
    """
    readings = {}
    for resolution in resolutions:
        model = build_model(resolution)
        model.initialise()
        readings[resolution] = []
        for activity in activities:
            model.sources[0].activity = activity
            reading = model.run()
            readings[resolution].append(reading)
            print(
                f"[{label}] resolution={resolution} activity={activity:.2e} Bq/m3 "
                f"I/I_ideal={reading.current / reading.ideal_current:.4f} "
                f"collection={reading.collection_efficiency:.4f} "
                f"space charge={reading.space_charge_field:.3g}",
                flush=True,
            )

    coarse_resolution, fine_resolution = resolutions
    ratio = fine_resolution / coarse_resolution
    rows = []
    for activity, coarse, fine in zip(
        activities, readings[coarse_resolution], readings[fine_resolution]
    ):
        row = {"activity": float(activity), **fine.as_dict()}
        row["current_coarse"] = coarse.current
        row["current_fine"] = fine.current
        row["current_extrapolated"] = richardson(
            coarse.current, fine.current, ratio=ratio
        )
        rows.append(row)
    return rows


def save_json(data: object, path: str | Path) -> None:
    """Writes results to a JSON file, creating the folder if needed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=1, default=float)
