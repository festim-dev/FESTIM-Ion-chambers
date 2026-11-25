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

## What it solves

For each charge carrier $s$, a steady-state balance of its density $n_s$ (m⁻³):

$$
\nabla \cdot J_s = S_s + R_s, \qquad
J_s = -D_s \nabla n_s + n_s v_s, \qquad
v_s = -z_s \, \mu_s \nabla \varphi
$$

| Symbol | Meaning | Unit |
|---|---|---|
| $J_s$ | flux of the carrier: diffusion plus drift | m⁻² s⁻¹ |
| $D_s$ | diffusion coefficient | m²/s |
| $v_s$ | drift velocity in the electric field | m/s |
| $z_s$ | charge number, +1 or −1 | |
| $\mu_s$ | mobility | m²/V/s |
| $\varphi$ | electric potential, from the electrostatic problem | V |
| $S_s$ | ionisation [source](sources.md), for the two ionisation products | m⁻³ s⁻¹ |
| $R_s$ | net production by the [reactions](gas.md#reactions) | m⁻³ s⁻¹ |

### Boundary conditions

With $n$ the normal pointing out of the gas:

| Condition | Equation | Class |
|---|---|---|
| Outflow | $J_s \cdot n = n_s \, (v_s \cdot n)$ | `festim.OutflowBC` |
| Zero density | $n_s = 0$ | `festim.FixedConcentrationBC` with `value=0.0` |
| Absorbing electrode | $J_s \cdot n = n_s \, \max(v_s \cdot n,\ 0)$ | {class}`~ion_chambers.AbsorbingElectrodeBC` |
| No flux | $J_s \cdot n = 0$ | nothing: any surface without a condition |

Outflow
: The carrier leaves at its drift velocity and the diffusive flux through the surface
  is zero. Where the field points the other way ($v_s \cdot n < 0$) the same condition
  lets the carrier in.

Zero density
: The surface removes every carrier that touches it.

Absorbing electrode
: Outflow where the field drives the carrier into the surface, no flux where it drives
  the carrier away.

No flux
: Nothing crosses. This is what an insulator gets.

[Electrode conditions](#electrode-conditions) says which of these each electrode gets
by default.

### Weak formulation

The densities are solved together with continuous linear finite elements, by FESTIM.
With $w_s$ the test function of carrier $s$:

$$
\int_\Omega \left( D_s \nabla n_s - n_s v_s \right) \cdot \nabla w_s \, dV
+ \int_{\Gamma_\mathrm{out}} n_s \, (v_s \cdot n) \, w_s \, dS
= \int_\Omega \left( S_s + R_s \right) w_s \, dV
$$

$\Gamma_\mathrm{out}$ is where the carrier has an outflow condition. For an absorbing
electrode, $v_s \cdot n$ in the surface term becomes $\max(v_s \cdot n,\ 0)$. A
surface without a condition contributes no term, which makes it a no-flux surface.
On a cylindrical mesh $dV = 2\pi r \, dr \, dz$ and $dS = 2\pi r \, dl$.

$D_s$ here is the diffusion coefficient after [stabilisation](#stabilisation).

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

| | Electrode above the mean potential | Electrode below the mean potential |
|---|---|---|
| Negative carrier | outflow | zero density |
| Positive carrier | zero density | outflow |

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
D = \max\left(D_\mathrm{phys},\ \tfrac{1}{2}\, a \,\mu\,|E|\,h_E\right)
$$

where $D_\mathrm{phys}$ is the diffusion coefficient of the gas data, $|E|$ the field
strength, $h_E$ the size of the cell along the field and $a$ the `factor` of the
stabilisation, 1 by default. This is the classic upwind
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
