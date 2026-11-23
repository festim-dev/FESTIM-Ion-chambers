"""Quick end-to-end check: the reference chamber on a coarse mesh, three activities.

    python scripts/quick_check.py

This script spells out every step, so it also serves as a first example.
"""

from reference_case import air, reference_chamber

import ion_chambers as ic

BIAS = 500.0  # V, potential of the anode rod
RADIAL_CELLS = 100  # mesh cells from the rod to the can
AXIAL_CELLS = 60  # mesh cells from one end cap to the other
ACTIVITIES = [3.7e4, 3.7e10, 3.7e13]  # Bq/m3

chamber = reference_chamber()

# ---- 1. the electric potential ---------------------------------------------------

electrostatics = ic.ElectrostaticProblem()
electrostatics.mesh = chamber.create_mesh(
    radial_cells=RADIAL_CELLS, axial_cells=AXIAL_CELLS
)
electrostatics.subdomains = [chamber.gas, chamber.anode, chamber.cathode]
electrostatics.boundary_conditions = [
    ic.FixedPotentialBC(subdomain=chamber.anode, value=BIAS),
    ic.FixedPotentialBC(subdomain=chamber.cathode, value=0.0),
]
electrostatics.initialise()
electrostatics.run()

# ---- 2. the transport of the charge carriers ---------------------------------------

tritium = ic.TritiumInGas(activity=0.0, geometry=chamber)

current = ic.ChamberCurrent(surface=chamber.anode)
anode_current = ic.ElectrodeCurrent(surface=chamber.anode)
ideal_current = ic.IdealCurrent()
wall_efficiency = ic.WallEfficiency()
collection_efficiency = ic.CollectionEfficiency()
space_charge_field = ic.SpaceChargeField()

transport = ic.ChargeTransportProblem()
transport.electrostatic_problem = electrostatics
transport.gas = air()
transport.ionisation_sources = [tritium]
transport.exports = [
    current,
    anode_current,
    ideal_current,
    wall_efficiency,
    collection_efficiency,
    space_charge_field,
]
transport.initialise()

# ---- 3. solve for each activity ---------------------------------------------------

for activity in ACTIVITIES:
    tritium.activity = activity
    transport.run()

    print(f"\nactivity {activity:.2e} Bq/m3")
    print(f"  chamber current        {current.value:.4e} A")
    print(f"  current through anode  {anode_current.value:.4e} A")
    print(f"  ideal current          {ideal_current.value:.4e} A")
    print(f"  wall efficiency        {wall_efficiency.value:.4f}")
    print(f"  collection efficiency  {collection_efficiency.value:.4f}")
    print(f"  space-charge field     {space_charge_field.value:.3g}")
