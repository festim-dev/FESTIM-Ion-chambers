"""The electric potential in the chamber.

:class:`ElectrostaticProblem` plays the role :class:`festim.HeatTransferProblem` plays
for hydrogen transport: it is solved on its own, and its solution is then given to the
transport problem.
"""

from __future__ import annotations

import festim as F
import ufl
from dolfinx import fem

from ion_chambers.constants import EPS0


class FixedPotentialBC(F.FixedTemperatureBC):
    """The electric potential of an electrode.

    It behaves exactly like :class:`festim.FixedTemperatureBC`, applied to the
    potential instead of the temperature.

    Args:
        subdomain: the surface of the electrode
        value: the potential (V)
    """


class PotentialExport(F.exports.FieldExportBase):
    """Export the electric potential to a file, for instance for ParaView.

    Args:
        filename: the name of the output file
        format: the output format, see :class:`festim.exports.FieldExportBase`
    """


class ElectrostaticProblem(F.ProblemBase):
    r"""The electric potential between the electrodes.

    Solves Poisson's equation

    .. math::

        -\nabla \cdot (\varepsilon_0 \nabla \varphi) = \rho

    with the potential fixed on the electrodes and no normal field on the other
    boundaries (insulators). Without a charge density this is Laplace's equation: the
    field of the empty chamber.

    It is used like a FESTIM problem: set the attributes, then call
    :meth:`initialise` and :meth:`run`. The solution is ``potential``.

    Args:
        mesh: the mesh
        subdomains: the gas volume and the surfaces of the electrodes
        boundary_conditions: the potential of each electrode, a list of
            :class:`FixedPotentialBC`
        charge_density: the charge density in the gas (C/m3), a ufl expression or a
            DOLFINx function. ``None`` for no charge.
        exports: the exports, a list of :class:`PotentialExport`
        settings: the solver settings. The default suits this linear problem.
        element_degree: the polynomial degree of the potential

    Attributes:
        mesh: the mesh
        subdomains: the gas volume and the surfaces of the electrodes
        boundary_conditions: the potential of each electrode
        charge_density: the charge density in the gas (C/m3), or ``None``
        exports: the exports
        settings: the solver settings
        element_degree: the polynomial degree of the potential
        u: the potential (V), after ``initialise``. Also available as ``potential``.

    Examples:
        .. code-block:: python

            import ion_chambers as ic

            chamber = ic.CoaxialChamber(
                rod_radius=1.5e-3, can_radius=0.05, length=0.127
            )

            electrostatics = ic.ElectrostaticProblem()
            electrostatics.mesh = chamber.create_mesh(
                radial_cells=120, axial_cells=72
            )
            electrostatics.subdomains = [chamber.gas, chamber.anode, chamber.cathode]
            electrostatics.boundary_conditions = [
                ic.FixedPotentialBC(subdomain=chamber.anode, value=500.0),
                ic.FixedPotentialBC(subdomain=chamber.cathode, value=0.0),
            ]
            electrostatics.initialise()
            electrostatics.run()
    """

    def __init__(
        self,
        mesh: F.Mesh | None = None,
        subdomains: list | None = None,
        boundary_conditions: list | None = None,
        charge_density=None,
        exports: list | None = None,
        settings: F.Settings | None = None,
        element_degree: int = 2,
    ) -> None:
        if settings is None:
            settings = F.Settings(atol=1e-10, rtol=1e-10, transient=False)
        super().__init__(
            mesh=mesh,
            subdomains=subdomains,
            boundary_conditions=boundary_conditions,
            exports=exports,
            settings=settings,
        )
        self.charge_density = charge_density
        self.element_degree = element_degree
        self.u = None

    @property
    def potential(self) -> fem.Function | None:
        """The electric potential (V). ``None`` before ``initialise``."""
        return self.u

    @property
    def boundary_conditions(self) -> list:
        return self._boundary_conditions

    @boundary_conditions.setter
    def boundary_conditions(self, value: list) -> None:
        for bc in value:
            if not isinstance(bc, FixedPotentialBC):
                raise TypeError(
                    f"boundary_conditions must be a list of "
                    f"ion_chambers.FixedPotentialBC, not {type(bc).__name__}"
                )
        self._boundary_conditions = value

    def initialise(self) -> None:
        """Creates the function space, the boundary conditions, the weak form and the
        solver."""
        self.define_function_space()
        self.define_meshtags_and_measures()
        self.t = fem.Constant(self.mesh.mesh, 0.0)
        self.define_boundary_conditions()
        self.create_formulation()
        self.create_solver()
        self.initialise_exports()

    def define_function_space(self) -> None:
        """Creates the function space, the potential ``u`` and the test function."""
        self.function_space = fem.functionspace(
            self.mesh.mesh, ("Lagrange", self.element_degree)
        )
        self.u = fem.Function(self.function_space, name="potential")
        self.test_function = ufl.TestFunction(self.function_space)

    def create_dirichletbc_form(self, bc: FixedPotentialBC):
        """Creates the DOLFINx boundary condition of an electrode.

        Args:
            bc: the boundary condition

        Returns:
            the DOLFINx boundary condition
        """
        bc.create_value(function_space=self.function_space, t=self.t)
        dofs = bc.define_surface_subdomain_dofs(
            facet_meshtags=self.facet_meshtags, function_space=self.function_space
        )
        if isinstance(bc.value_fenics, fem.Function):
            return fem.dirichletbc(value=bc.value_fenics, dofs=dofs)
        return fem.dirichletbc(value=bc.value_fenics, dofs=dofs, V=self.function_space)

    def create_formulation(self) -> None:
        """Creates the weak form of Poisson's equation, divided by epsilon_0."""
        # the factor 2 pi r of a cylindrical mesh, or 1 on a cartesian mesh
        weight = self.mesh.coordinate_system.integration_weight(self.mesh.mesh)

        self.formulation = (
            weight * ufl.dot(ufl.grad(self.u), ufl.grad(self.test_function)) * self.dx
        )
        if self.charge_density is None:
            # Laplace's equation. A zero source keeps the form well defined.
            zero = fem.Constant(self.mesh.mesh, 0.0)
            self.formulation -= weight * zero * self.test_function * self.dx
        else:
            self.formulation -= (
                weight * self.charge_density / EPS0 * self.test_function * self.dx
            )

    def initialise_exports(self) -> None:
        """Opens the file of each export."""
        for export in self.exports:
            if not isinstance(export, PotentialExport):
                raise TypeError(
                    f"exports must be a list of ion_chambers.PotentialExport, "
                    f"not {type(export).__name__}"
                )
            export.define_writer([self.u], ["potential"], self.mesh.mesh)

    def post_processing(self) -> None:
        """Writes the exports."""
        for export in self.exports:
            export.update()
            export.write(float(self.t))
