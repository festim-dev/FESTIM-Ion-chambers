"""Ionisation sources: where and how fast ion pairs are created."""

from abc import ABC, abstractmethod

import numpy as np
from dolfinx import fem
from dolfinx import mesh as dmesh
from scipy.interpolate import RegularGridInterpolator

from ion_chambers.beta import ion_pairs_per_decay
from ion_chambers.gas import Gas
from ion_chambers.geometry import Geometry
from ion_chambers.projection import lumped_projection
from ion_chambers.wall_loss import efficiency_grid


class IonisationSource(ABC):
    """Base class of the ionisation sources.

    A model sums the contributions of all its sources, so tritium in the gas and
    tritium held on the walls can be described side by side.
    """

    @abstractmethod
    def build(self, mesh: dmesh.Mesh, geometry: Geometry, gas: Gas):
        """Creates the source term on the mesh.

        Args:
            mesh: the mesh of the gas volume
            geometry: the chamber geometry
            gas: the gas

        Returns:
            the production rate of ion pairs (1/m3/s), a ufl expression
        """

    @abstractmethod
    def ideal_pair_rate(self, geometry: Geometry, gas: Gas) -> float:
        """The pair production rate (1/s) if no beta energy were lost to the walls."""


class TritiumInGas(IonisationSource):
    r"""Tritium mixed uniformly into the gas.

    The local production rate of ion pairs is

    .. math::

        S(r, z) = A \, \frac{\bar{E}}{W} \, f(r, z)

    with :math:`A` the activity concentration, :math:`\bar{E}` the mean beta energy,
    :math:`W` the energy per ion pair in the gas and :math:`f` the wall-loss
    efficiency (see :mod:`ion_chambers.wall_loss`).

    The activity can be changed after the model is initialised; the next
    :meth:`ion_chambers.IonChamberModel.run` uses the new value without rebuilding
    anything.

    Args:
        activity: the activity concentration (Bq/m3)
        wall_loss: whether betas lose energy to the walls. ``False`` sets f = 1.
        n_rays: the number of betas followed per point in the wall-loss Monte Carlo

    Attributes:
        activity: the activity concentration (Bq/m3)
        wall_loss: whether betas lose energy to the walls
        n_rays: the number of betas followed per point in the wall-loss Monte Carlo
        efficiency: the wall-loss efficiency on the mesh, after ``build``
    """

    def __init__(
        self, activity: float, wall_loss: bool = True, n_rays: int = 3000
    ) -> None:
        self._rate: fem.Constant | None = None
        self._pairs_per_decay: float | None = None
        self.activity = activity
        self.wall_loss = wall_loss
        self.n_rays = n_rays
        self.efficiency: fem.Function | None = None

    @property
    def activity(self) -> float:
        return self._activity

    @activity.setter
    def activity(self, value: float) -> None:
        if isinstance(value, bool) or not isinstance(value, int | float | np.number):
            raise TypeError(f"activity must be a float, not {type(value).__name__}")
        if value < 0:
            raise ValueError(f"activity must not be negative, not {value}")
        self._activity = float(value)
        if self._rate is not None:
            self._rate.value = self._activity * self._pairs_per_decay

    def build(self, mesh: dmesh.Mesh, geometry: Geometry, gas: Gas):
        self._pairs_per_decay = ion_pairs_per_decay(gas.w_value)
        self._rate = fem.Constant(mesh, self.activity * self._pairs_per_decay)

        V1 = fem.functionspace(mesh, ("Lagrange", 1))
        self.efficiency = fem.Function(V1, name="wall_loss_efficiency")
        if self.wall_loss:
            self.efficiency.x.array[:] = self._project_efficiency(
                mesh, V1, geometry, gas
            )
        else:
            self.efficiency.x.array[:] = 1.0
        return self._rate * self.efficiency

    def _project_efficiency(
        self, mesh: dmesh.Mesh, V1: fem.FunctionSpace, geometry: Geometry, gas: Gas
    ) -> np.ndarray:
        r, z, f = efficiency_grid(geometry, gas.density, n_rays=self.n_rays)
        interpolate = RegularGridInterpolator(
            (r, z), f, bounds_error=False, fill_value=None
        )

        def efficiency(r_points: np.ndarray, z_points: np.ndarray) -> np.ndarray:
            inside = np.c_[
                np.clip(r_points, r[0], r[-1]), np.clip(z_points, z[0], z[-1])
            ]
            return np.clip(interpolate(inside), 0, 1)

        return lumped_projection(mesh, V1, efficiency, geometry.coordinate_system)

    def ideal_pair_rate(self, geometry: Geometry, gas: Gas) -> float:
        return self.activity * ion_pairs_per_decay(gas.w_value) * geometry.volume
