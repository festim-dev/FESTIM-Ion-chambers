# Gas

The gas is described by what carries charge in it and how those carriers react. All of
it is input that you control.

The package ships no gas data. Transport and reaction coefficients depend on the gas,
its purity, the pressure and the field, and a value that looks official is easy to
trust too much. So the gas is always something you define.
The [example case](../example_case.md#the-script) has a complete definition to
copy and adapt.

## Defining a gas

A {class}`~ion_chambers.Gas` is built from {class}`~ion_chambers.ChargeCarrier` and
{class}`~ion_chambers.CarrierReaction` objects. This example is a gas in which
electrons do not attach, so there are no negative ions at all:

```python
import ion_chambers as ic

electrons = ic.ChargeCarrier(name="e", charge=-1, mobility=0.3, diffusivity=0.3)
positive_ions = ic.ChargeCarrier(name="P", charge=1, mobility=1.5e-4)

recombination = ic.CarrierReaction(reactant=[positive_ions, electrons], rate=1e-13)

gas = ic.Gas(
    name="my gas",
    w_value=26.4,  # eV per ion pair
    density=1.66,  # kg/m3, sets the range of the betas
    pressure=101325.0,  # Pa
    temperature=293.15,  # K
    carriers=[electrons, positive_ions],
    ionisation_products=(electrons, positive_ions),
    reactions=[recombination],
)
```

The numbers above only illustrate the syntax.

Carriers
: Each has a name, a charge number and a mobility in m²/V/s. If you leave out
  `diffusivity`, the Einstein relation `D = mobility × k_B T / e` is used. That is
  right for ions, but not for electrons, which the field heats well above the gas
  temperature: give their diffusivity explicitly.

Reactions
: Each has one or more reactants, a rate coefficient and, optionally, products. The
  rate of the reaction is the coefficient times the densities of the reactants, so
  the coefficient is in 1/s for one reactant and m³/s for two. A reaction with no
  product is a recombination and counts as a loss of signal.

Ionisation products
: The two carriers that each ionisation creates, one negative and one positive.

## Reactions

A reaction turns its reactants into its product at a rate set by a coefficient and
the densities of the reactants. With one reactant:

$$
\mathrm{A} \xrightarrow{\;k\;} \mathrm{C}, \qquad r = k \, n_\mathrm{A}
$$

and with two:

$$
\mathrm{A} + \mathrm{B} \xrightarrow{\;k\;} \mathrm{C}, \qquad
r = k \, n_\mathrm{A} \, n_\mathrm{B}
$$

$r$ is the number of reactions per unit volume and time (m⁻³ s⁻¹), and $n$ a
density (m⁻³). Each reactant loses $r$ and the product gains $r$. Without a product
the charges neutralise and are gone: a recombination.

The three reactions of air, with electrons $\mathrm{e^-}$, positive ions
$\mathrm{P^+}$ and negative ions $\mathrm{N^-}$:

| | Reaction | Rate | Coefficient |
|---|---|---|---|
| Attachment | $\mathrm{e^-} + \mathrm{O_2} \rightarrow \mathrm{N^-}$ | $\nu_a \, n_e$ | $\nu_a$ (1/s) |
| Ion–ion recombination | $\mathrm{P^+} + \mathrm{N^-} \rightarrow \varnothing$ | $\alpha \, n_P \, n_N$ | $\alpha$ (m³/s) |
| Electron–ion recombination | $\mathrm{P^+} + \mathrm{e^-} \rightarrow \varnothing$ | $\beta \, n_P \, n_e$ | $\beta$ (m³/s) |

```python
attachment = ic.CarrierReaction(reactant=electrons, product=negative_ions, rate=1.0e8)
ion_ion = ic.CarrierReaction(reactant=[positive_ions, negative_ions], rate=2.0e-12)
electron_ion = ic.CarrierReaction(reactant=[positive_ions, electrons], rate=2.0e-13)
```

The oxygen is not a carrier: it is part of the gas, and its density is folded into
the attachment frequency $\nu_a$.

Summed over the reactions, these rates give the net production $R_s$ of each carrier
in the [transport equation](transport.md#what-it-solves). For air:

$$
\begin{aligned}
R_e &= - \nu_a n_e - \beta \, n_P n_e \\
R_P &= - \alpha \, n_P n_N - \beta \, n_P n_e \\
R_N &= \nu_a n_e - \alpha \, n_P n_N
\end{aligned}
$$

## Relation to FESTIM

`ChargeCarrier` is a child of `festim.Species`, and `CarrierReaction` is a child of
`festim.GenericReaction`. The transport problem passes them to FESTIM as they are, and
their arguments (`reactant`, `product`) are named as in FESTIM.

They add what a chamber needs:

- A `ChargeCarrier` has a `charge`, a `mobility` and a `diffusivity`.
- A `CarrierReaction` is created without a volume subdomain. The transport problem
  assigns its own gas volume when it is initialised. Its `rate` may depend on the reduced
  electric field, which FESTIM's rate coefficients do not support.

Like any FESTIM species, a carrier holds the solution of the problem it belongs to. Do
not share one gas between two transport problems. Build a gas for each, for instance
with a function that returns a new `Gas`.

## Pressure and temperature

The coefficients you give are used as they are. `Gas` does not rescale them with
`pressure` or `temperature`. Those two only set the number density (for E/N) and the
thermal voltage (for the Einstein relation).

To study a range of pressures, write a function that takes the pressure and returns a
`Gas` with the coefficients scaled the way you judge right for your gas.
The [example case](../example_case.md#the-script) does exactly that.

## Coefficients that depend on the field

Any mobility, diffusivity or rate can be a function of the reduced electric field E/N
instead of a number. The function receives E/N in townsend (1 Td = 10⁻²¹ V m²) as a
numpy array and returns the coefficient in SI units:

```python
import numpy as np


def electron_mobility(reduced_field):
    """Illustration only: mobility falling with E/N."""
    return 0.3 / (1 + np.sqrt(reduced_field))


electron = ic.ChargeCarrier(
    name="e", charge=-1, mobility=electron_mobility, diffusivity=0.3
)
```

This is how tabulated swarm data (for example from a Boltzmann solver) enters the
model: wrap the table in a function that interpolates it.

After `transport.initialise()`, `transport.reduced_field` holds E/N in every cell,
which helps to check that your data covers the range the chamber actually sees.

## Impurities

There is no separate impurity object. An impurity that captures electrons, such as a
trace of O₂ or water in N₂, is described by adding a negative-ion carrier and an
attachment `CarrierReaction` whose rate reflects the impurity concentration.
