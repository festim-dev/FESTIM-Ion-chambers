"""A material with one diffusivity field per species."""

from __future__ import annotations

import festim as F
from dolfinx import fem


class CellwiseDiffusionMaterial(F.Material):
    """Material whose diffusion coefficient is a separate field for each species.

    :class:`festim.Material` accepts a field through ``D``, but uses that one field for
    every species. Here each charge carrier needs its own.

    Args:
        D_fields: the diffusion coefficient of each species (m2/s), keyed by species
            name
        name: the name of the material
    """

    def __init__(self, D_fields: dict[str, fem.Function], name: str | None = None):
        super().__init__(
            D_0=dict.fromkeys(D_fields, 1.0),
            E_D=dict.fromkeys(D_fields, 0.0),
            name=name,
        )
        self.D_fields = D_fields

    def get_diffusion_coefficient(self, mesh=None, temperature=None, species=None):
        return self.D_fields[species.name]
