"""Ionisation sources: where and how fast ion pairs are created."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable

import festim as F
import numpy as np
from dolfinx import fem
from scipy.interpolate import RegularGridInterpolator

from ion_chambers.beta import efficiency_grid, ion_pairs_per_decay
from ion_chambers.gas import Gas
from ion_chambers.geometry import Geometry


class IonisationSource(ABC):
    """Base class of the ionisation sources.

    A transport problem sums the contributions of all its sources, so tritium in the
    gas and tritium held on the walls can be described side by side.
    """

    @abstractmethod
    def create_pair_production(self, mesh: F.Mesh, gas: Gas):
        """Creates the source term on the mesh.

        Args:
            mesh: the mesh of the gas volume
            gas: the gas

        Returns:
            the production rate of ion pairs (1/m3/s), a ufl expression
        """

    @abstractmethod
    def ideal_pair_rate(self, gas: Gas, gas_volume: float) -> float:
        """The pair production rate (1/s) if no beta energy were lost to the walls.

        Args:
            gas: the gas
            gas_volume: the volume of gas in the chamber (m3)
        """


class TritiumInGas(IonisationSource):
    r"""Tritium mixed uniformly into the gas.

    The local production rate of ion pairs is

    .. math::

        S(r, z) = A \, \frac{\bar{E}}{W} \, f(r, z)

    with :math:`A` the activity concentration, :math:`\bar{E}` the mean beta energy,
    :math:`W` the energy per ion pair in the gas and :math:`f` the wall-loss
    efficiency (see :mod:`ion_chambers.beta`).

    The activity can be changed after the transport problem is initialised; the next
    ``run()`` uses the new value without rebuilding anything.

    Args:
        activity: the activity concentration (Bq/m3)
        geometry: the chamber geometry, used to compute how much beta energy is lost
            to the walls. ``None`` ignores wall losses (f = 1).
        n_rays: the number of betas followed per point in the wall-loss Monte Carlo

    Attributes:
        activity: the activity concentration (Bq/m3)
        geometry: the chamber geometry, or ``None``
        n_rays: the number of betas followed per point in the wall-loss Monte Carlo
        efficiency: the wall-loss efficiency on the mesh, once the source is created
    """

    def __init__(
        self,
        activity: float,
        geometry: Geometry | None = None,
        n_rays: int = 3000,
    ) -> None:
        self._rate = None
        self._pairs_per_decay = None
        self.activity = activity
        self.geometry = geometry
        self.n_rays = n_rays
        self.efficiency = None

    @property
    def activity(self) -> float:
        return self._activity

    @activity.setter
    def activity(self, value: float) -> None:
        is_number = isinstance(value, int | float | np.number)
        if not is_number or isinstance(value, bool):
            raise TypeError(f"activity must be a float, not {type(value).__name__}")
        if value < 0:
            raise ValueError(f"activity must not be negative, not {value}")
        self._activity = float(value)
        # if the source is already part of a problem, update its fenics constant
        if self._rate is not None:
            self._rate.value = self._activity * self._pairs_per_decay

    def create_pair_production(self, mesh: F.Mesh, gas: Gas):
        self._pairs_per_decay = ion_pairs_per_decay(gas.w_value)
        self._rate = fem.Constant(mesh.mesh, self.activity * self._pairs_per_decay)

        V1 = fem.functionspace(mesh.mesh, ("Lagrange", 1))
        self.efficiency = fem.Function(V1, name="wall_loss_efficiency")
        if self.geometry is None:
            self.efficiency.x.array[:] = 1.0
        else:
            self.efficiency.x.array[:] = self._efficiency_at_nodes(mesh, V1, gas)

        return self._rate * self.efficiency

    def _efficiency_at_nodes(
        self, mesh: F.Mesh, V1: fem.FunctionSpace, gas: Gas
    ) -> np.ndarray:
        """The wall-loss efficiency at the nodes of the mesh."""
        r, z, efficiency = efficiency_grid(
            self.geometry, gas.density, n_rays=self.n_rays
        )
        interpolator = RegularGridInterpolator(
            (r, z), efficiency, bounds_error=False, fill_value=None
        )

        def efficiency_at(r_points: np.ndarray, z_points: np.ndarray) -> np.ndarray:
            # points a rounding error outside the grid are moved onto its edge
            r_inside = np.clip(r_points, r[0], r[-1])
            z_inside = np.clip(z_points, z[0], z[-1])
            values = interpolator(np.column_stack([r_inside, z_inside]))
            return np.clip(values, 0, 1)

        is_cylindrical = str(mesh.coordinate_system) == "cylindrical"
        return lumped_projection(V1, efficiency_at, weight_by_radius=is_cylindrical)

    def ideal_pair_rate(self, gas: Gas, gas_volume: float) -> float:
        return self.activity * ion_pairs_per_decay(gas.w_value) * gas_volume


def lumped_projection(
    V1: fem.FunctionSpace,
    function: Callable[[np.ndarray, np.ndarray], np.ndarray],
    weight_by_radius: bool = True,
    subdivisions: int = 6,
) -> np.ndarray:
    """Nodal values ``f_i = int f phi_i dV / int phi_i dV`` of a function.

    ``f`` is sampled on ``subdivisions**2`` points per cell, so the projection keeps
    the integral of a layer thinner than the cells, where nodal interpolation would
    smear it into the neighbouring cells.

    Only for P1 on triangles, on one MPI rank.

    Args:
        V1: the P1 function space
        function: the function, called with the two coordinates as flat arrays
        weight_by_radius: whether ``dV`` carries the factor r of a cylindrical mesh
        subdivisions: the number of subdivisions of each cell edge

    Returns:
        the nodal values, ordered like the dofs of ``V1``
    """
    mesh = V1.mesh
    vertices = mesh.geometry.x[:, :2][mesh.geometry.dofmaps[0]]  # (cells, 3, 2)
    dofs = V1.dofmap.list  # same local vertex order as the P1 geometry

    # Split each triangle into subdivisions**2 small triangles and take their
    # centroids as sample points, in barycentric coordinates (l1, l2).
    k = subdivisions
    barycentric = []
    for i in range(k):
        for j in range(k - i):
            barycentric.append(((i + 1 / 3) / k, (j + 1 / 3) / k))  # pointing up
            if i + j < k - 1:
                barycentric.append(((i + 2 / 3) / k, (j + 2 / 3) / k))  # pointing down
    barycentric = np.array(barycentric)
    n_samples = len(barycentric)

    # the three P1 basis functions at the sample points
    basis = np.column_stack([1 - barycentric.sum(axis=1), barycentric])
    samples = np.einsum("qv,cvd->cqd", basis, vertices)  # (cells, samples, 2)

    edge_1 = vertices[:, 1] - vertices[:, 0]
    edge_2 = vertices[:, 2] - vertices[:, 0]
    cell_area = 0.5 * np.abs(edge_1[:, 0] * edge_2[:, 1] - edge_1[:, 1] * edge_2[:, 0])
    weights = np.repeat(cell_area[:, None] / n_samples, n_samples, axis=1)
    if weight_by_radius:
        weights = weights * samples[..., 0]

    values = function(samples[..., 0].ravel(), samples[..., 1].ravel())
    values = values.reshape(samples.shape[:2])

    index_map = V1.dofmap.index_map
    numerator = np.zeros(index_map.size_local + index_map.num_ghosts)
    denominator = np.zeros_like(numerator)
    np.add.at(numerator, dofs, np.einsum("cq,qv->cv", weights * values, basis))
    np.add.at(denominator, dofs, np.einsum("cq,qv->cv", weights, basis))
    return numerator / denominator
