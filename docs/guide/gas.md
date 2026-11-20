# Gas

The gas is described by what carries charge in it and how those carriers react. All of
it is input that you control.

## Using the example gas

```python
import ion_chambers as ic

gas = ic.library.air()  # 1 atm, 293.15 K
gas = ic.library.air(pressure=50_000.0)  # Pa
gas = ic.library.air(temperature=320.0)  # K
```

{func}`~ion_chambers.library.air` has three carriers (electrons, positive ions and
O₂⁻ negative ions) and three reactions (attachment of electrons to O₂, ion–ion
recombination and electron–ion recombination). Its docstring lists every value and how
each one is scaled with pressure and temperature.

:::{admonition} The library holds examples, not reference data
:class: warning
The values in `ion_chambers.library` are reasonable starting points with simple,
documented scaling rules. Check them against your conditions before relying on a
result. Air is the only gas provided so far.
:::

## Defining your own gas

A {class}`~ion_chambers.Gas` is built from {class}`~ion_chambers.ChargeCarrier` and
{class}`~ion_chambers.Reaction` objects. This example is a gas in which electrons do
not attach, so there are no negative ions at all:

```python
gas = ic.Gas(
    name="my gas",
    w_value=26.4,  # eV per ion pair
    density=1.66,  # kg/m3, sets the range of the betas
    pressure=101325.0,  # Pa
    temperature=293.15,  # K
    carriers=[
        ic.ChargeCarrier(name="e", charge=-1, mobility=0.3, diffusivity=0.3),
        ic.ChargeCarrier(name="P", charge=1, mobility=1.5e-4),
    ],
    reactions=[
        ic.Reaction(reactants=["P", "e"], rate=1e-13),  # recombination, m3/s
    ],
    ionisation_products=("e", "P"),
)
```

The numbers above only illustrate the syntax.

Carriers
: Each has a name, a charge number and a mobility in m²/V/s. If you leave out
  `diffusivity`, the Einstein relation `D = mobility × k_B T / e` is used. That is
  right for ions, but not for electrons, which the field heats well above the gas
  temperature: give their diffusivity explicitly.

Reactions
: Each has a list of reactants, a rate coefficient and a list of products. The rate of
  the reaction is the coefficient times the densities of the reactants, so the
  coefficient is in 1/s for one reactant and m³/s for two. A reaction with no product
  is a recombination and counts as a loss of signal.

Ionisation products
: The two carriers that each ionisation creates.

The coefficients you give are used as they are. `Gas` does not rescale them with
`pressure` or `temperature`; those two set the number density (for E/N) and the
thermal voltage (for the Einstein relation).

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

After `model.initialise()`, `model.reduced_field` holds E/N in every cell, which helps
to check that your data covers the range the chamber actually sees.

## Impurities

There is no separate impurity object. An impurity that captures electrons, such as a
trace of O₂ or water in N₂, is described by adding a negative-ion carrier and an
attachment reaction whose rate reflects the impurity concentration.
