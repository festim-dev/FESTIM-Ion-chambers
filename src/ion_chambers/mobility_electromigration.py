"""Electromigration with the drift set by a mobility."""

from __future__ import annotations

import festim as F
import ufl
from dolfinx import fem


class MobilityElectromigration(F.ElectromigrationTerm):
    """Electromigration with the drift velocity ``v = -z mu grad(phi)``.

    :class:`festim.ElectromigrationTerm` derives the mobility from the diffusion
    coefficient with the Einstein relation. That is wrong here for two reasons: the
    diffusion coefficient carries a numerical floor (see
    :class:`ion_chambers.DiffusionFloor`), and electrons heated by the field do not
    follow the Einstein relation at the gas temperature.

    Args:
        species: the species the drift acts on
        charge: the charge number of the species
        mobility: the mobility (m2/V/s), a constant or a field
        potential: the electric potential (V)
        subdomain: the volume subdomain where the drift applies
    """

    def __init__(
        self,
        species: F.Species,
        charge: int,
        mobility: float | fem.Function,
        potential: fem.Function,
        subdomain: F.VolumeSubdomain,
    ) -> None:
        super().__init__(
            species=species, charge=charge, potential=potential, subdomain=subdomain
        )
        self.mobility = mobility

    def drift_velocity(self, D, temperature):
        return -self.charge * self.mobility * ufl.grad(self.potential.fenics_object)
