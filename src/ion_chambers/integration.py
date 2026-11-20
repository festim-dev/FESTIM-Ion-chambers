"""Integration helpers that hide the coordinate system and the MPI reduction."""

from __future__ import annotations

from mpi4py import MPI

import numpy as np
import ufl
from dolfinx import fem
from dolfinx import mesh as dmesh


def volume_weight(mesh: dmesh.Mesh, coordinate_system: str):
    """The factor turning an integral over the mesh into one over the real volume.

    An axisymmetric (r, z) mesh stands for a body of revolution, so its integrals
    carry ``2 pi r``. A cartesian mesh needs nothing.

    Args:
        mesh: the mesh
        coordinate_system: ``"cylindrical"`` or ``"cartesian"``

    Returns:
        the weight, to multiply the integrand with
    """
    if coordinate_system == "cylindrical":
        return 2 * np.pi * ufl.SpatialCoordinate(mesh)[0]
    if coordinate_system == "cartesian":
        return 1.0
    raise ValueError(
        f"coordinate_system must be 'cylindrical' or 'cartesian', "
        f"not '{coordinate_system}'"
    )


def assemble(form: ufl.Form) -> float:
    """Assembles a scalar form and sums it over the MPI ranks."""
    compiled = fem.form(form)
    return compiled.mesh.comm.allreduce(fem.assemble_scalar(compiled), op=MPI.SUM)
