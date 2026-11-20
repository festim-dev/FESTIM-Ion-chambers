"""Stabilisation of drift-dominated transport."""

from __future__ import annotations

import numpy as np
import ufl
from dolfinx import fem
from dolfinx import mesh as dmesh


def cell_field_length(
    mesh: dmesh.Mesh, potential: fem.Function
) -> tuple[np.ndarray, np.ndarray]:
    """The field strength in each cell and the length of the cell along the field.

    Only for triangular meshes.

    Args:
        mesh: a triangular mesh
        potential: the electric potential (V)

    Returns:
        ``|E|`` at the cell centroids (V/m) and the cell length along E (m)
    """
    vertices = mesh.geometry.x[:, :2][mesh.geometry.dofmaps[0]]  # (cells, 3, 2)
    centroids = vertices.mean(axis=1)
    points = np.c_[centroids, np.zeros(len(centroids))]
    cells = np.arange(len(centroids), dtype=np.int32)

    V = fem.functionspace(mesh, ("DG", 1, (2,)))
    field = fem.Function(V)
    field.interpolate(
        fem.Expression(-ufl.grad(potential), V.element.interpolation_points)
    )
    field_at_centroids = field.eval(points, cells)
    strength = np.linalg.norm(field_at_centroids, axis=1)
    direction = field_at_centroids / np.maximum(strength[:, None], 1e-300)

    edges = np.stack(
        [
            vertices[:, 1] - vertices[:, 0],
            vertices[:, 2] - vertices[:, 1],
            vertices[:, 0] - vertices[:, 2],
        ],
        axis=1,
    )
    length = np.max(np.abs(np.einsum("cej,cj->ce", edges, direction)), axis=1)
    return strength, length


class DiffusionFloor:
    r"""A cell-local lower bound on the diffusion coefficient.

    FESTIM does not stabilise drift-dominated transport, and the cell Péclet number of
    the ions reaches about 100. Each carrier therefore gets

    .. math::

        D = \max\left(D_\mathrm{phys},\ \tfrac{1}{2}\, \mathrm{factor}\, \mu |E| h_E
        \right)

    with :math:`h_E` the cell length along the field. This is first-order upwinding
    written as diffusion. The floor vanishes where the field is weak, so it cannot
    fake losses in low-field corners, and the error it leaves is first order in the
    cell size (see :func:`ion_chambers.richardson`).

    Args:
        factor: the multiplier of the floor. 1 gives a cell Péclet number of 2; 0
            switches the stabilisation off.
    """

    def __init__(self, factor: float = 1.0) -> None:
        if factor < 0:
            raise ValueError(f"factor must not be negative, not {factor}")
        self.factor = factor

    def apply(
        self,
        diffusivity: np.ndarray,
        mobility: np.ndarray,
        field_strength: np.ndarray,
        cell_length: np.ndarray,
    ) -> np.ndarray:
        """Returns the stabilised diffusion coefficient in each cell.

        Args:
            diffusivity: the physical diffusion coefficient in each cell (m2/s)
            mobility: the mobility in each cell (m2/V/s)
            field_strength: ``|E|`` in each cell (V/m)
            cell_length: the cell length along the field (m)
        """
        floor = self.factor * 0.5 * mobility * field_strength * cell_length
        return np.maximum(diffusivity, floor)
