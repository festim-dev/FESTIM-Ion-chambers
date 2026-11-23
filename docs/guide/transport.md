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
