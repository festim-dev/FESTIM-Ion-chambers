# The reading

`model.run()` returns a {class}`~ion_chambers.Reading`.

```python
reading = model.run()
reading.current  # A
reading.as_dict()  # everything, for saving
```

## What the chamber reads

`current`
: The chamber current in amperes, from the Shockley–Ramo theorem. This is the number
  to compare with a measurement. Its sign follows the bias.

`ideal_current`
: The current if every beta stopped in the gas and every charge reached an electrode:
  `e × activity × volume × pairs per decay`.

The ratio of the two is the product of the two efficiencies below, up to
discretisation error.

`wall_efficiency`
: The share of the beta energy deposited in the gas rather than in the walls.

`collection_efficiency`
: The share of the ion pairs that escape recombination.

`generation`, `recombination`
: The underlying rates, in ion pairs per second over the whole gas volume.

## How far to trust it

`anode_current`, `cathode_current`, `balance_current`
: Three other estimates of the same current: from the carrier fluxes into each
  electrode, and from the volume balance of generation minus recombination. In exact
  arithmetic all four agree. Their spread shows the discretisation error of the mesh
  you used. The Ramo current and the anode flux typically agree to better than 0.01%.

`space_charge_field`
: The field the charges in the gas would create, relative to the applied field. The
  model assumes this is small. Above about 0.1 the result is indicative only.

`space_charge_potential`
: The same check expressed as the largest space-charge potential over the bias.

`mean_net_charge_density`
: The volume average of the absolute net charge density, in elementary charges per m³.

## Why the Ramo current

The current through an electrode can be computed directly from the flux of carriers
into it, but that needs the gradient of the densities at the boundary, where the
numerical error is largest. The Shockley–Ramo theorem gives the same current as an
integral over the whole gas volume, which is much less sensitive to the mesh.
