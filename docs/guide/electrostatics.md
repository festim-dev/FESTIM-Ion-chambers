# Electric potential

The charge carriers drift in an electric potential. That potential is the solution of
an {class}`~ion_chambers.ElectrostaticProblem`, which is solved first and then given to
the transport problem.

It plays the role that `festim.HeatTransferProblem` plays for hydrogen transport, and
is used the same way.

```python
import ion_chambers as ic

chamber = ic.CoaxialChamber(rod_radius=1.5e-3, can_radius=0.05, length=0.127)

electrostatics = ic.ElectrostaticProblem()
electrostatics.mesh = chamber.create_mesh(radial_cells=120, axial_cells=72)
electrostatics.subdomains = [chamber.gas, chamber.anode, chamber.cathode]
electrostatics.boundary_conditions = [
    ic.FixedPotentialBC(subdomain=chamber.anode, value=500.0),
    ic.FixedPotentialBC(subdomain=chamber.cathode, value=0.0),
]
electrostatics.initialise()
electrostatics.run()

potential = electrostatics.potential  # a DOLFINx function, in volts
```

## What it solves

Poisson's equation, with the potential fixed on each electrode and no normal field on
the other boundaries, which are the insulators:

$$
-\nabla \cdot (\varepsilon_0 \nabla \varphi) = \rho
$$

By default there is no charge density, so this is Laplace's equation: the field of the
empty chamber.

## Electrodes

Each electrode is a surface subdomain with a {class}`~ion_chambers.FixedPotentialBC`.
`FixedPotentialBC` is a child of `festim.FixedTemperatureBC` and takes the same
arguments.

A negative potential on the anode rod makes it collect the positive ions instead:

```python
ic.FixedPotentialBC(subdomain=chamber.anode, value=-500.0)
```

The transport problem reads the electrodes and their potentials from the electrostatic
problem. It works for two electrodes.

## Space charge

The transport problem takes the potential as fixed. It ignores the field created by the
charges in the gas themselves. That is a good approximation at low activity and breaks
down at high activity, where slow ions accumulate.

The export {class}`~ion_chambers.SpaceChargeField` tells you where you stand. It solves
a second electrostatic problem with the charge density of the solution as source and
grounded electrodes, and reports the size of that field relative to the applied one.
Treat results as indicative once it exceeds about 0.1. For the reference chamber at
500 V that happens near 37 Ci/m³.

You can also give a charge density yourself:

```python
electrostatics.charge_density = my_charge_density  # C/m3, a ufl expression
```

Iterating between the two problems, so that the space charge feeds back into the
potential, is on the [roadmap](../roadmap.md).

## Writing the potential to a file

```python
electrostatics.exports = [ic.PotentialExport(filename="results/potential.bp")]
```

The file opens in ParaView. `PotentialExport` is a child of FESTIM's field export and
accepts the same `format` argument.
