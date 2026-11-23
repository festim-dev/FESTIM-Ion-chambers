"""Helpers for parameter studies: mesh extrapolation, response curves, saving."""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path

from ion_chambers.charge_transport_problem import ChargeTransportProblem
from ion_chambers.exports import ChamberCurrent, ChamberQuantity


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


def find_chamber_current(problem: ChargeTransportProblem) -> ChamberCurrent:
    """Returns the ChamberCurrent among the exports of a problem."""
    for export in problem.exports:
        if isinstance(export, ChamberCurrent):
            return export
    raise ValueError("the problem needs an ion_chambers.ChamberCurrent in its exports")


def response_curve(
    coarse_problem: ChargeTransportProblem,
    fine_problem: ChargeTransportProblem,
    activities: Sequence[float],
    cell_size_ratio: float = 2.0,
    label: str = "",
) -> list[dict[str, float]]:
    """The chamber current against tritium activity, extrapolated from two meshes.

    The two problems must describe the same chamber on a coarse and on a fine mesh,
    and each must have a :class:`ion_chambers.ChamberCurrent` among its exports. Each
    is initialised once, then solved for every activity by changing the activity of
    its first ionisation source.

    Args:
        coarse_problem: the problem on the coarse mesh, not yet initialised
        fine_problem: the problem on the fine mesh, not yet initialised
        activities: the activity concentrations (Bq/m3)
        cell_size_ratio: the cell size of the coarse mesh divided by that of the fine
            mesh
        label: a prefix for the progress lines

    Returns:
        one row per activity, with ``activity``, ``current_coarse``, ``current_fine``,
        ``current_extrapolated`` and the value of every chamber quantity exported by
        the fine problem, keyed by its title
    """
    _solve_for_each_activity(coarse_problem, activities, f"{label} coarse")
    _solve_for_each_activity(fine_problem, activities, f"{label} fine")

    coarse_current = find_chamber_current(coarse_problem)
    fine_current = find_chamber_current(fine_problem)

    rows = []
    for i, activity in enumerate(activities):
        row = {"activity": float(activity)}
        row["current_coarse"] = coarse_current.data[i]
        row["current_fine"] = fine_current.data[i]
        row["current_extrapolated"] = richardson(
            coarse_current.data[i], fine_current.data[i], ratio=cell_size_ratio
        )
        for export in fine_problem.exports:
            if isinstance(export, ChamberQuantity):
                row[export.title] = export.data[i]
        rows.append(row)
    return rows


def _solve_for_each_activity(
    problem: ChargeTransportProblem, activities: Sequence[float], label: str
) -> None:
    """Initialises a problem and runs it at each activity."""
    current = find_chamber_current(problem)
    problem.initialise()
    for activity in activities:
        problem.ionisation_sources[0].activity = activity
        problem.run()
        print(
            f"[{label}] activity={activity:.2e} Bq/m3 current={current.value:.4e} A",
            flush=True,
        )


def save_json(data: object, path: str | Path) -> None:
    """Writes results to a JSON file, creating the folder if needed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=1, default=float)
