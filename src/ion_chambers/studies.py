"""Helpers for parameter studies: mesh extrapolation, response curves, saving."""

import json
from collections.abc import Sequence
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
    coarse_model: IonChamberModel,
    fine_model: IonChamberModel,
    activities: Sequence[float],
    label: str = "",
) -> list[dict[str, float]]:
    """The chamber current against tritium activity, extrapolated from two meshes.

    The two models must describe the same chamber and differ only in their
    ``resolution``. Each is initialised once, then solved for every activity by
    changing the activity of its first source.

    Args:
        coarse_model: the model on the coarse mesh, not yet initialised
        fine_model: the model on the fine mesh, not yet initialised
        activities: the activity concentrations (Bq/m3)
        label: a prefix for the progress lines

    Returns:
        one row per activity: the fine-mesh reading as a dictionary, plus ``activity``,
        ``current_coarse``, ``current_fine`` and ``current_extrapolated``
    """
    coarse_readings = _solve_for_each_activity(coarse_model, activities, label)
    fine_readings = _solve_for_each_activity(fine_model, activities, label)
    cell_size_ratio = fine_model.resolution / coarse_model.resolution

    rows = []
    for i, activity in enumerate(activities):
        coarse = coarse_readings[i]
        fine = fine_readings[i]
        row = fine.as_dict()
        row["activity"] = float(activity)
        row["current_coarse"] = coarse.current
        row["current_fine"] = fine.current
        row["current_extrapolated"] = richardson(
            coarse.current, fine.current, ratio=cell_size_ratio
        )
        rows.append(row)
    return rows


def _solve_for_each_activity(
    model: IonChamberModel, activities: Sequence[float], label: str
) -> list:
    """Initialises a model and returns its reading at each activity."""
    model.initialise()
    readings = []
    for activity in activities:
        model.sources[0].activity = activity
        reading = model.run()
        readings.append(reading)
        print(
            f"[{label}] resolution={model.resolution} activity={activity:.2e} Bq/m3 "
            f"I/I_ideal={reading.current / reading.ideal_current:.4f} "
            f"collection={reading.collection_efficiency:.4f} "
            f"space charge={reading.space_charge_field:.3g}",
            flush=True,
        )
    return readings


def save_json(data: object, path: str | Path) -> None:
    """Writes results to a JSON file, creating the folder if needed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=1, default=float)
