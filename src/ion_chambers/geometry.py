"""Chamber geometries.

A geometry knows where the walls of a chamber are, which the wall-loss calculation
needs. :class:`CoaxialChamber` also creates its own mesh and subdomains. A mesh made
with another tool is used as in FESTIM, without a geometry: see the documentation.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

from mpi4py import MPI

import festim as F
import numpy as np
from dolfinx import mesh as dmesh
from scipy.optimize import brentq


class Geometry(ABC):
    """Base class of the chamber geometries: what the wall-loss calculation needs.

    Attributes:
        dimensions: the numbers that define the shape, as a tuple. Two geometries with
            the same dimensions share the cached wall-loss calculation.
    """

    dimensions: tuple

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

    The model is axisymmetric: the mesh covers the (r, z) half section of the gas.

    Args:
        rod_radius: the radius of the anode rod (m)
        can_radius: the inner radius of the can (m)
        length: the inner length of the can (m)
        feedthrough_gap: the radial width of the insulating rings (m)

    Attributes:
        rod_radius: the radius of the anode rod (m)
        can_radius: the inner radius of the can (m)
        length: the inner length of the can (m)
        feedthrough_gap: the radial width of the insulating rings (m)
        gas: the gas volume, a :class:`festim.VolumeSubdomain` with id 1
        anode: the surface of the rod, a :class:`festim.SurfaceSubdomain` with id 1
        cathode: the can and the end caps beyond the insulating rings, a
            :class:`festim.SurfaceSubdomain` with id 2
    """

    def __init__(
        self,
        rod_radius: float,
        can_radius: float,
        length: float,
        feedthrough_gap: float = 3e-3,
    ) -> None:
        dimensions = {
            "rod_radius": rod_radius,
            "can_radius": can_radius,
            "length": length,
            "feedthrough_gap": feedthrough_gap,
        }
        for name, value in dimensions.items():
            is_number = isinstance(value, int | float) and not isinstance(value, bool)
            if not is_number:
                raise TypeError(f"{name} must be a float, not {type(value).__name__}")
            if value <= 0:
                raise ValueError(f"{name} must be positive, not {value}")
        if rod_radius >= can_radius:
            raise ValueError(
                f"rod_radius ({rod_radius}) must be smaller than can_radius "
                f"({can_radius})"
            )
        if rod_radius + feedthrough_gap >= can_radius:
            raise ValueError("feedthrough_gap leaves no cathode on the end caps")

        self.rod_radius = rod_radius
        self.can_radius = can_radius
        self.length = length
        self.feedthrough_gap = feedthrough_gap

        self.dimensions = (rod_radius, can_radius, length)
        # The gas has no material here: the transport problem defines the transport
        # properties of the gas once the electric field is known.
        self.gas = F.VolumeSubdomain(id=1, material=None)
        self.anode = F.SurfaceSubdomain(id=1, locator=self.is_anode)
        self.cathode = F.SurfaceSubdomain(id=2, locator=self.is_cathode)

    @classmethod
    def from_volume(
        cls,
        volume: float,
        rod_radius: float,
        can_radius: float,
        feedthrough_gap: float = 3e-3,
    ) -> CoaxialChamber:
        """A chamber of given gas volume and radii; the length follows.

        Args:
            volume: the gas volume (m3)
            rod_radius: the radius of the anode rod (m)
            can_radius: the inner radius of the can (m)
            feedthrough_gap: the radial width of the insulating rings (m)
        """
        length = volume / (np.pi * (can_radius**2 - rod_radius**2))
        return cls(
            rod_radius=rod_radius,
            can_radius=can_radius,
            length=length,
            feedthrough_gap=feedthrough_gap,
        )

    @classmethod
    def from_aspect_ratio(
        cls,
        volume: float,
        rod_radius: float,
        aspect_ratio: float,
        feedthrough_gap: float = 3e-3,
    ) -> CoaxialChamber:
        """A chamber of given gas volume and length-to-diameter ratio.

        Args:
            volume: the gas volume (m3)
            rod_radius: the radius of the anode rod (m)
            aspect_ratio: the length divided by the diameter of the can
            feedthrough_gap: the radial width of the insulating rings (m)
        """

        def excess_volume(can_radius: float) -> float:
            length = 2 * can_radius * aspect_ratio
            return np.pi * (can_radius**2 - rod_radius**2) * length - volume

        can_radius = brentq(excess_volume, rod_radius * (1 + 1e-9), 1.0)
        length = 2 * can_radius * aspect_ratio
        return cls(
            rod_radius=rod_radius,
            can_radius=can_radius,
            length=length,
            feedthrough_gap=feedthrough_gap,
        )

    @property
    def volume(self) -> float:
        """The gas volume (m3)."""
        return np.pi * (self.can_radius**2 - self.rod_radius**2) * self.length

    def create_mesh(self, radial_cells: int, axial_cells: int) -> F.Mesh:
        """Creates a triangular (r, z) mesh of the gas, with cells of equal size.

        Args:
            radial_cells: the number of cells from the rod to the can
            axial_cells: the number of cells from one end cap to the other

        Returns:
            the mesh, in cylindrical coordinates
        """
        for name, value in {
            "radial_cells": radial_cells,
            "axial_cells": axial_cells,
        }.items():
            is_integer = isinstance(value, int) and not isinstance(value, bool)
            if not is_integer:
                raise TypeError(f"{name} must be an int, not {type(value).__name__}")
            if value <= 0:
                raise ValueError(f"{name} must be positive, not {value}")

        mesh = dmesh.create_rectangle(
            MPI.COMM_WORLD,
            [[self.rod_radius, 0.0], [self.can_radius, self.length]],
            [radial_cells, axial_cells],
            dmesh.CellType.triangle,
            diagonal=dmesh.DiagonalType.crossed,
        )
        return F.Mesh(mesh, coordinate_system="cylindrical")

    def is_anode(self, x: np.ndarray) -> np.ndarray:
        """Marks the points on the rod.

        Args:
            x: the coordinates, one column per point (the DOLFINx locator convention)
        """
        tolerance = 1e-9 * self.can_radius
        return np.isclose(x[0], self.rod_radius, atol=tolerance)

    def is_cathode(self, x: np.ndarray) -> np.ndarray:
        """Marks the points on the can and on the end caps beyond the insulating rings.

        Args:
            x: the coordinates, one column per point (the DOLFINx locator convention)
        """
        tolerance = 1e-9 * self.can_radius
        on_can = np.isclose(x[0], self.can_radius, atol=tolerance)
        on_bottom_cap = np.isclose(x[1], 0, atol=tolerance)
        on_top_cap = np.isclose(x[1], self.length, atol=tolerance)
        beyond_insulator = x[0] >= self.rod_radius + self.feedthrough_gap - tolerance
        return on_can | ((on_bottom_cap | on_top_cap) & beyond_insulator)

    def distance_to_wall(self, points: np.ndarray, directions: np.ndarray):
        r = points[:, 0]
        z = points[:, 1]
        ux = directions[:, 0]
        uy = directions[:, 1]
        uz = directions[:, 2]

        to_can = _distance_to_cylinder(r, ux, uy, self.can_radius, inside=True)
        to_rod = _distance_to_cylinder(r, ux, uy, self.rod_radius, inside=False)

        with np.errstate(divide="ignore"):
            to_top_cap = (self.length - z[:, None]) / uz[None, :]
            to_bottom_cap = -z[:, None] / uz[None, :]
        to_caps = np.full_like(to_can, np.inf)
        going_up = np.broadcast_to(uz[None, :] > 0, to_caps.shape)
        going_down = np.broadcast_to(uz[None, :] < 0, to_caps.shape)
        to_caps[going_up] = to_top_cap[going_up]
        to_caps[going_down] = to_bottom_cap[going_down]

        return np.minimum(np.minimum(to_can, to_rod), to_caps)

    def wall_clustered_grid(self, nr: int, nz: int) -> tuple[np.ndarray, np.ndarray]:
        gap = self.can_radius - self.rod_radius
        s = np.linspace(-1, 1, nr)
        r = self.rod_radius + 0.5 * gap * (1 + np.sin(0.5 * np.pi * s))
        t = np.linspace(-1, 1, nz)
        z = 0.5 * self.length * (1 + np.sin(0.5 * np.pi * t))
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
    # the ray hits the cylinder where A s^2 + B s + C = 0
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
