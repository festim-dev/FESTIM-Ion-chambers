"""Electric potential solvers.

The transport model takes the potential as a given field. Which equation that field
solves is a choice, made by passing a :class:`PotentialSolver` to the model.
"""

from abc import ABC, abstractmethod

import ufl
from dolfinx import fem
from dolfinx import mesh as dmesh
from dolfinx.fem.petsc import LinearProblem

from ion_chambers.constants import EPS0
from ion_chambers.geometry import Geometry
from ion_chambers.integration import volume_weight

_DIRECT_SOLVER = {"ksp_type": "preonly", "pc_type": "lu"}


def _electrode_dofs(V: fem.FunctionSpace, locator) -> list:
    fdim = V.mesh.topology.dim - 1
    facets = dmesh.locate_entities_boundary(V.mesh, fdim, locator)
    return fem.locate_dofs_topological(V, fdim, facets)


class PotentialSolver(ABC):
    """Base class of the potential solvers."""

    @abstractmethod
    def solve(self, mesh: dmesh.Mesh, geometry: Geometry, bias: float) -> fem.Function:
        """Computes the electric potential.

        Args:
            mesh: the mesh of the gas volume
            geometry: the chamber geometry, which locates the electrodes
            bias: the potential of the anode (V); the cathode is at 0 V

        Returns:
            the potential (V)
        """


class LaplaceField(PotentialSolver):
    """The potential of the empty chamber: Laplace's equation, no space charge.

    Valid while the charge in the gas is too small to distort the applied field. The
    reading reports how far that holds (see
    :attr:`ion_chambers.Reading.space_charge_field`).

    Args:
        degree: the polynomial degree of the potential
    """

    def __init__(self, degree: int = 2) -> None:
        self.degree = degree

    def solve(self, mesh: dmesh.Mesh, geometry: Geometry, bias: float) -> fem.Function:
        V = fem.functionspace(mesh, ("Lagrange", self.degree))
        bcs = [
            fem.dirichletbc(
                fem.Constant(mesh, float(bias)),
                _electrode_dofs(V, geometry.is_anode),
                V,
            ),
            fem.dirichletbc(
                fem.Constant(mesh, 0.0), _electrode_dofs(V, geometry.is_cathode), V
            ),
        ]
        u, w = ufl.TrialFunction(V), ufl.TestFunction(V)
        weight = volume_weight(mesh, geometry.coordinate_system)
        problem = LinearProblem(
            weight * ufl.dot(ufl.grad(u), ufl.grad(w)) * ufl.dx,
            fem.Constant(mesh, 0.0) * w * ufl.dx,
            bcs=bcs,
            petsc_options_prefix="potential_",
            petsc_options=_DIRECT_SOLVER,
        )
        return problem.solve()


class SpaceChargePerturbation:
    r"""The potential created by the charge in the gas, with grounded electrodes.

    Solves :math:`-\nabla \cdot (\varepsilon_0 \nabla \varphi_1) = \rho`. Comparing
    :math:`\nabla \varphi_1` with the applied field tells how much the space charge,
    which :class:`LaplaceField` ignores, would distort it.

    Args:
        mesh: the mesh of the gas volume
        geometry: the chamber geometry, which locates the electrodes
        degree: the polynomial degree of the potential
    """

    def __init__(self, mesh: dmesh.Mesh, geometry: Geometry, degree: int = 2) -> None:
        V = fem.functionspace(mesh, ("Lagrange", degree))
        self.charge_density = fem.Function(fem.functionspace(mesh, ("Lagrange", 1)))

        def on_electrode(x):
            return geometry.is_anode(x) | geometry.is_cathode(x)

        bc = fem.dirichletbc(
            fem.Constant(mesh, 0.0), _electrode_dofs(V, on_electrode), V
        )
        u, w = ufl.TrialFunction(V), ufl.TestFunction(V)
        weight = volume_weight(mesh, geometry.coordinate_system)
        self._problem = LinearProblem(
            weight * ufl.dot(ufl.grad(u), ufl.grad(w)) * ufl.dx,
            weight * self.charge_density / EPS0 * w * ufl.dx,
            bcs=[bc],
            petsc_options_prefix="space_charge_",
            petsc_options=_DIRECT_SOLVER,
        )

    def solve(self, charge_density) -> fem.Function:
        """Computes the perturbation potential.

        Args:
            charge_density: the charge density (C/m3), a ufl expression

        Returns:
            the perturbation potential (V)
        """
        points = self.charge_density.function_space.element.interpolation_points
        self.charge_density.interpolate(fem.Expression(charge_density, points))
        return self._problem.solve()
