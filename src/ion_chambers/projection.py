"""Conservative projection of a point-wise function onto P1."""

from collections.abc import Callable

import numpy as np
from dolfinx import fem
from dolfinx import mesh as dmesh


def lumped_projection(
    mesh: dmesh.Mesh,
    V1: fem.FunctionSpace,
    func: Callable[[np.ndarray, np.ndarray], np.ndarray],
    coordinate_system: str = "cylindrical",
    k: int = 6,
) -> np.ndarray:
    """Nodal values ``f_i = int f phi_i dV / int phi_i dV`` of a function.

    ``f`` is sampled on ``k**2`` points per cell, so the projection keeps the integral
    of a layer thinner than the cells, where nodal interpolation would smear it into
    the neighbouring cells. ``dV`` carries the factor r on a cylindrical mesh.

    Only for P1 on triangles, on one MPI rank.

    Args:
        mesh: a triangular mesh
        V1: the P1 space on ``mesh``
        func: the function, called with the two coordinates as flat arrays
        coordinate_system: ``"cylindrical"`` or ``"cartesian"``
        k: the number of subdivisions of each cell edge

    Returns:
        the nodal values, ordered like the dofs of ``V1``
    """
    vertices = mesh.geometry.x[:, :2][mesh.geometry.dofmaps[0]]  # (cells, 3, 2)
    dofs = V1.dofmap.list  # same local vertex order as the P1 geometry

    # barycentric coordinates (l1, l2) of the centroids of the k**2 sub-triangles
    barycentric = []
    for i in range(k):
        for j in range(k - i):
            barycentric.append(((i + 1 / 3) / k, (j + 1 / 3) / k))  # pointing up
            if i + j < k - 1:
                barycentric.append(((i + 2 / 3) / k, (j + 2 / 3) / k))  # pointing down
    barycentric = np.array(barycentric)
    basis = np.c_[1 - barycentric.sum(1), barycentric]  # phi_i at the sample points
    samples = np.einsum("qv,cvd->cqd", basis, vertices)  # (cells, samples, 2)

    edge_1, edge_2 = vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]
    area = 0.5 * np.abs(edge_1[:, 0] * edge_2[:, 1] - edge_1[:, 1] * edge_2[:, 0])
    weights = np.repeat(area[:, None] / len(barycentric), len(barycentric), axis=1)
    if coordinate_system == "cylindrical":
        weights = weights * samples[..., 0]

    values = func(samples[..., 0].ravel(), samples[..., 1].ravel())
    values = values.reshape(samples.shape[:2])

    index_map = V1.dofmap.index_map
    numerator = np.zeros(index_map.size_local + index_map.num_ghosts)
    denominator = np.zeros_like(numerator)
    np.add.at(numerator, dofs, np.einsum("cq,qv->cv", weights * values, basis))
    np.add.at(denominator, dofs, np.einsum("cq,qv->cv", weights, basis))
    return numerator / denominator
