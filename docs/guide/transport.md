# Charge transport

{class}`~ion_chambers.ChargeTransportProblem` solves for the density of every charge
carrier of the gas. It is a child of `festim.HydrogenTransportProblem`.

```python
transport = ic.ChargeTransportProblem()
transport.electrostatic_problem = electrostatics  # already run
transport.gas = gas
transport.ionisation_sources = [tritium]
transport.exports = [current]
transport.initialise()
transport.run()
```

## What you give it

`electrostatic_problem`
: An {class}`~ion_chambers.ElectrostaticProblem` that has been run. The transport
  problem takes its potential, its mesh and its electrodes from it.

`gas`
: The [gas](gas.md): its charge carriers and their reactions.

`ionisation_sources`
: The [sources](sources.md) of ion pairs.

`exports`
: The [quantities](exports.md) to compute after each run.

## What it sets up for you

During `initialise()`, the problem builds the FESTIM model from those inputs:

- the carriers become the species, and the gas reactions become the reactions
- each carrier gets a drift term in the electric potential
- each ionisation source feeds the two ionisation products of the gas
- each electrode absorbs the carriers the field drives into it and emits none; the
  insulators let nothing through

It also evaluates the mobilities, diffusivities and reaction rates in the local
electric field, for the coefficients that depend on E/N.

A problem is initialised once. To change the gas, the potential or the mesh, create new
problems. Only the activity of a source can be changed between two calls to `run()`.

## Adding your own FESTIM objects

`sources`, `drift_terms` and `boundary_conditions` can be set as on any FESTIM problem.
`initialise()` adds its own to what you set:

```python
# a flow of gas that carries the positive ions along
flow = F.AdvectionTerm(
    velocity=gas_velocity, subdomain=chamber.gas, species=positive_ions
)
transport.drift_terms = [flow]

# a different condition for the electrons on the anode
transport.boundary_conditions = [
    F.FixedConcentrationBC(subdomain=chamber.anode, value=0.0, species=electrons)
]
```

A boundary condition that you set for a carrier on an electrode replaces the automatic
one for that carrier and electrode. Sources and drift terms are given on the volume
subdomain of the electrostatic problem.

`mesh`, `subdomains`, `species`, `reactions` and `temperature` must not be set on the
transport problem: they come from the electrostatic problem and the gas, and
`initialise()` raises a `ValueError` if they are set.

## Electrode conditions

By default each electrode gets, for each carrier, one of two conditions, chosen by
comparing the potential of the electrode with the mean of the electrode potentials:
outflow for the carriers it attracts, zero density for those it repels.

{class}`~ion_chambers.AbsorbingElectrodeBC` is an alternative that needs no such rule.
It looks at the force on the carrier at each point of the electrode: where the field
drives the carrier into the surface it leaves at its drift velocity, and where the
field drives it away nothing crosses. The same condition is given to every carrier
on every electrode:

```python
for electrode in [chamber.anode, chamber.cathode]:
    for carrier in gas.carriers:
        transport.boundary_conditions.append(
            ic.AbsorbingElectrodeBC(subdomain=electrode, species=carrier)
        )
```

The two differ on a single mesh. With the default, the zero density on the repelling
electrode combines with the [stabilisation](#stabilisation) into a sink, and the
current comes out too low by an amount that falls with the cell size (about 6% at 30
radial and 18 axial cells for the reference chamber at trace level). With absorbing
electrodes the same mesh gives the current of all the pairs created to about 0.02%. Both converge to the same current.
The absorbing condition neglects the carriers that diffuse against the field to a
repelling electrode, a share of about the thermal voltage over the bias. It is not
the default yet: the pinned results of this package were computed with the first.

## Stabilisation

Ions in a chamber are carried by the field far more than they diffuse. Standard finite
elements produce oscillations in that regime, and FESTIM does not stabilise it. The
problem therefore raises the diffusion coefficient in each cell to at least

$$
D = \max\left(D_\mathrm{phys},\ \tfrac{1}{2}\,\mu\,|E|\,h_E\right)
$$

where $h_E$ is the size of the cell along the field. This is the classic upwind
scheme, written as extra diffusion. Two properties matter:

- The extra diffusion vanishes where the field is weak, so it does not create
  artificial losses in the low-field corners of a chamber.
- The error it leaves shrinks in proportion to the cell size. This is why currents are
  [extrapolated from two meshes](studies.md#mesh-extrapolation).

The strength can be changed, mainly to test how sensitive a result is:

```python
transport.stabilisation = ic.DiffusionFloor(factor=0.5)
```

## Inspecting the problem

After `initialise()`:

| Attribute | Content |
|---|---|
| `potential` | the electric potential (V) |
| `reduced_field` | E/N in every cell (Td) |
| `mobilities` | the mobility of each carrier, keyed by name |
| `diffusivities` | the stabilised diffusion coefficient of each carrier, keyed by name |
| `rates` | the rate coefficient of each reaction |
| `volume_of_gas` | the gas volume (m³) |

`reduced_field` helps to check that field-dependent gas data covers the range the
chamber actually sees.
