"""Chamber geometries.

A geometry knows four things: how to mesh the gas volume, where the electrodes are,
how large the gas volume is, and how far a beta travels before it hits a wall. The
rest of the package only talks to a geometry through :class:`Geometry`, so a new shape
is a new subclass.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from mpi4py import MPI

import numpy as np
from dolfinx import mesh as dmesh
from scipy.optimize import brentq


class Geometry(ABC):
    """Base class of the chamber geometries.

    Attributes:
        coordinate_system: ``"cylindrical"`` for an axisymmetric (r, z) mesh or
            ``"cartesian"``
    """

    coordinate_system: str = "cylindrical"

    @abstractmethod
    def build_mesh(self, resolution: int) -> dmesh.Mesh:
        """Creates the mesh of the gas volume.

        Args:
            resolution: the number of cells across the electrode gap
        """

    @abstractmethod
    def is_anode(self, x: np.ndarray) -> np.ndarray:
        """Marks the points on the anode, the electrode held at the bias.

        Args:
            x: the coordinates, one column per point (the DOLFINx locator convention)
        """

    @abstractmethod
    def is_cathode(self, x: np.ndarray) -> np.ndarray:
        """Marks the points on the cathode, the grounded electrode.

        Args:
            x: the coordinates, one column per point (the DOLFINx locator convention)
        """

    @property
    @abstractmethod
    def volume(self) -> float:
        """The gas volume (m3)."""

    @abstractmethod
    def distance_to_wall(self, points: np.ndarray, directions: np.ndarray):
        """The straight-line distance from each point to the nearest wall hit.

        Args:
            points: the starting points in mesh coordinates, shape (n, 2)
            directions: unit vectors in 3D, shape (m, 3)

        Returns:
            the distances (m), shape (n, m)
        """

    @abstractmethod
    def wall_clustered_grid(self, nr: int, nz: int) -> tuple[np.ndarray, np.ndarray]:
        """The axes of a tensor grid covering the gas, refined towards every wall.

        Used to tabulate the wall-loss efficiency, which varies over a few millimetres
        next to the walls and is flat elsewhere.
        """


@dataclass(frozen=True)
class CoaxialChamber(Geometry):
    """An anode rod on the axis of a closed cylindrical can.

    The can and its two end caps are the cathode. Where the rod passes through the
    caps, an insulating ring of width ``feedthrough_gap`` separates it from the cap::

            z = length  +--+=====================+
                        |  : insulator   cathode |
                        |a |                     |
                        |n |        gas          | cathode
                        |o |                     |
                        |d |                     |
                        |e : insulator   cathode |
            z = 0       +--+=====================+
                     r = rod_radius          r = can_radius

    Args:
        rod_radius: the radius of the anode rod (m)
        can_radius: the inner radius of the can (m)
        length: the inner length of the can (m)
        feedthrough_gap: the radial width of the insulating rings (m)
        axial_cells: the number of mesh cells along z. ``None`` derives it from the
            radial resolution.
        cap_refinement: how strongly the mesh is refined towards the end caps (0 for a
            uniform mesh)
    """

    rod_radius: float
    can_radius: float
    length: float
    feedthrough_gap: float = 3e-3
    axial_cells: int | None = field(default=None, compare=False)
    cap_refinement: float = field(default=1.6, compare=False)

    def __post_init__(self) -> None:
        for attribute in ("rod_radius", "can_radius", "length", "feedthrough_gap"):
            value = getattr(self, attribute)
            if isinstance(value, bool) or not isinstance(value, int | float):
                raise TypeError(
                    f"{attribute} must be a float, not {type(value).__name__}"
                )
            if value <= 0:
                raise ValueError(f"{attribute} must be positive, not {value}")
        if self.rod_radius >= self.can_radius:
            raise ValueError(
                f"rod_radius ({self.rod_radius}) must be smaller than can_radius "
                f"({self.can_radius})"
            )
        if self.rod_radius + self.feedthrough_gap >= self.can_radius:
            raise ValueError("feedthrough_gap leaves no cathode on the end caps")

    @classmethod
    def from_volume(
        cls, volume: float, rod_radius: float, can_radius: float, **kwargs
    ) -> CoaxialChamber:
        """A chamber of given gas volume and radii; the length follows.

        Args:
            volume: the gas volume (m3)
            rod_radius: the radius of the anode rod (m)
            can_radius: the inner radius of the can (m)
            **kwargs: passed on to :class:`CoaxialChamber`
        """
        length = volume / (np.pi * (can_radius**2 - rod_radius**2))
        return cls(rod_radius, can_radius, length, **kwargs)

    @classmethod
    def from_aspect_ratio(
        cls, volume: float, rod_radius: float, aspect_ratio: float, **kwargs
    ) -> CoaxialChamber:
        """A chamber of given gas volume and length-to-diameter ratio.

        Args:
            volume: the gas volume (m3)
            rod_radius: the radius of the anode rod (m)
            aspect_ratio: the length divided by the diameter of the can
            **kwargs: passed on to :class:`CoaxialChamber`
        """

        def excess_volume(can_radius: float) -> float:
            length = 2 * can_radius * aspect_ratio
            return np.pi * (can_radius**2 - rod_radius**2) * length - volume

        can_radius = brentq(excess_volume, rod_radius * (1 + 1e-9), 1.0)
        return cls(rod_radius, can_radius, 2 * can_radius * aspect_ratio, **kwargs)

    @property
    def volume(self) -> float:
        return np.pi * (self.can_radius**2 - self.rod_radius**2) * self.length

    @property
    def _tolerance(self) -> float:
        return 1e-9 * self.can_radius

    def default_axial_cells(self, radial_cells: int) -> int:
        """The number of cells along z that matches a radial resolution.

        The spacing at mid-height is about three times the outer radial spacing,
        bounded to between 0.6 and 1 times the radial count.
        """
        a, R = self.rod_radius, self.can_radius
        outer_spacing = R * np.log(R / a) / radial_cells
        cells = self.length / outer_spacing * 0.35
        return int(np.clip(cells, 0.6 * radial_cells, 1.0 * radial_cells))

    def build_mesh(self, resolution: int) -> dmesh.Mesh:
        """Creates a triangular (r, z) mesh, graded towards the rod and the caps.

        The radial spacing grows in proportion to r, which keeps the product of drift
        speed and cell size uniform in the 1/r field.

        Args:
            resolution: the number of cells in the radial direction
        """
        a, R, L = self.rod_radius, self.can_radius, self.length
        axial_cells = self.axial_cells or self.default_axial_cells(resolution)
        mesh = dmesh.create_rectangle(
            MPI.COMM_WORLD,
            [[0, -1], [1, 1]],
            [resolution, axial_cells],
            dmesh.CellType.triangle,
            diagonal=dmesh.DiagonalType.crossed,
        )
        x = mesh.geometry.x
        s, t = x[:, 0].copy(), x[:, 1].copy()
        x[:, 0] = a * (R / a) ** s
        gamma = self.cap_refinement
        stretched = np.tanh(gamma * t) / np.tanh(gamma) if gamma > 0 else t
        x[:, 1] = 0.5 * L * (1 + stretched)
        return mesh

    def is_anode(self, x: np.ndarray) -> np.ndarray:
        return np.isclose(x[0], self.rod_radius, atol=self._tolerance)

    def is_cathode(self, x: np.ndarray) -> np.ndarray:
        tol = self._tolerance
        on_cap = np.isclose(x[1], 0, atol=tol) | np.isclose(x[1], self.length, atol=tol)
        beyond_insulator = x[0] >= self.rod_radius + self.feedthrough_gap - tol
        return np.isclose(x[0], self.can_radius, atol=tol) | (on_cap & beyond_insulator)

    def distance_to_wall(self, points: np.ndarray, directions: np.ndarray):
        r, z = points[:, 0], points[:, 1]
        ux, uy, uz = directions[:, 0], directions[:, 1], directions[:, 2]
        to_can = _distance_to_cylinder(r, ux, uy, self.can_radius, inside=True)
        to_rod = _distance_to_cylinder(r, ux, uy, self.rod_radius, inside=False)
        with np.errstate(divide="ignore"):
            to_caps = np.where(
                uz[None, :] > 0,
                (self.length - z[:, None]) / uz[None, :],
                np.where(uz[None, :] < 0, -z[:, None] / uz[None, :], np.inf),
            )
        return np.minimum(np.minimum(to_can, to_rod), to_caps)

    def wall_clustered_grid(self, nr: int, nz: int) -> tuple[np.ndarray, np.ndarray]:
        a, R, L = self.rod_radius, self.can_radius, self.length
        s = np.linspace(-1, 1, nr)
        r = a + 0.5 * (R - a) * (1 + np.sin(0.5 * np.pi * s))
        t = np.linspace(-1, 1, nz)
        z = 0.5 * L * (1 + np.sin(0.5 * np.pi * t))
        return r, z


def _distance_to_cylinder(
    r0: np.ndarray, ux: np.ndarray, uy: np.ndarray, radius: float, inside: bool
) -> np.ndarray:
    """Distance along rays from the points (r0, 0) to a cylinder about the z axis.

    Args:
        r0: the radial positions of the starting points, shape (n,)
        ux: the x components of the ray directions, shape (m,)
        uy: the y components of the ray directions, shape (m,)
        radius: the radius of the cylinder
        inside: True when the points are inside the cylinder (the ray leaves through
            it), False when they are outside a solid rod (the ray may miss it)

    Returns:
        the distances, shape (n, m), infinite for a ray that misses
    """
    A = ux**2 + uy**2
    B = 2 * r0[:, None] * ux[None, :]
    C = (r0**2 - radius**2)[:, None]
    discriminant = B**2 - 4 * A[None, :] * C
    with np.errstate(invalid="ignore", divide="ignore"):
        root = np.sqrt(np.where(discriminant > 0, discriminant, np.nan))
        if inside:
            s = (-B + root) / (2 * A[None, :])
        else:
            s = (-B - root) / (2 * A[None, :])
            s = np.where(s > 0, s, np.inf)
    return np.where(np.isnan(s), np.inf, s)
