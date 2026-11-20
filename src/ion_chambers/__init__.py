"""FESTIM model of a tritium ionisation chamber (process monitor)."""

# mpi4py first: it has to initialise MPI before FESTIM pulls in DOLFINx
from mpi4py import MPI as _MPI  # noqa: F401

import festim as _festim

from ion_chambers._festim_check import check_festim as _check_festim

_check_festim(_festim)

from ion_chambers import beta, constants, library, wall_loss
from ion_chambers.cellwise_diffusion_material import CellwiseDiffusionMaterial
from ion_chambers.field import LaplaceField, PotentialSolver, SpaceChargePerturbation
from ion_chambers.gas import ChargeCarrier, Gas, Reaction
from ion_chambers.geometry import CoaxialChamber, Geometry
from ion_chambers.mobility_electromigration import MobilityElectromigration
from ion_chambers.model import IonChamberModel
from ion_chambers.readout import Reading
from ion_chambers.source import IonisationSource, TritiumInGas
from ion_chambers.stabilisation import DiffusionFloor
from ion_chambers.studies import response_curve, richardson, save_json

__all__ = [
    "CellwiseDiffusionMaterial",
    "ChargeCarrier",
    "CoaxialChamber",
    "DiffusionFloor",
    "Gas",
    "Geometry",
    "IonChamberModel",
    "IonisationSource",
    "LaplaceField",
    "MobilityElectromigration",
    "PotentialSolver",
    "Reaction",
    "Reading",
    "SpaceChargePerturbation",
    "TritiumInGas",
    "beta",
    "constants",
    "library",
    "response_curve",
    "richardson",
    "save_json",
    "wall_loss",
]
