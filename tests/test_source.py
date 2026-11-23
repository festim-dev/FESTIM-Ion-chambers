import numpy as np
import pytest
from dolfinx import fem
from prototype_case import air

import ion_chambers as ic
from ion_chambers.source import lumped_projection


@pytest.mark.parametrize("activity", ["high", None, True])
def test_activity_rejects_wrong_types(activity):
    """An activity that is not a number raises a TypeError."""
    with pytest.raises(TypeError, match="activity must be a float"):
        ic.TritiumInGas(activity=activity)


def test_activity_rejects_negative_values():
    """A negative activity raises a ValueError."""
    with pytest.raises(ValueError, match="activity must not be negative"):
        ic.TritiumInGas(activity=-1.0)


def test_ideal_pair_rate():
    """The ideal pair rate is the activity times the gas volume times the ion pairs
    per decay."""
    source = ic.TritiumInGas(activity=3.7e10)
    expected = 3.7e10 * 1e-3 * ic.beta.ion_pairs_per_decay(33.97)
    assert np.isclose(source.ideal_pair_rate(air(), gas_volume=1e-3), expected)


def test_without_geometry_the_efficiency_is_one(chamber):
    """A source created without a geometry ignores wall losses: the efficiency is
    exactly one everywhere."""
    source = ic.TritiumInGas(activity=1.0)
    source.create_pair_production(chamber.create_mesh(6), air())
    assert np.all(source.efficiency.x.array == 1.0)


def test_wall_loss_efficiency_stays_between_zero_and_one(chamber):
    """The projected wall-loss efficiency stays within its physical bounds of 0 and
    1."""
    source = ic.TritiumInGas(activity=1.0, geometry=chamber, n_rays=200)
    source.create_pair_production(chamber.create_mesh(6), air())
    values = source.efficiency.x.array
    assert values.min() >= 0.0 and values.max() <= 1.0


def test_changing_the_activity_updates_the_created_source(chamber):
    """Setting a new activity on a source that is part of a problem updates its rate
    constant without a rebuild."""
    source = ic.TritiumInGas(activity=1.0)
    source.create_pair_production(chamber.create_mesh(6), air())
    source.activity = 10.0
    assert np.isclose(source._rate.value, 10.0 * ic.beta.ion_pairs_per_decay(33.97))


def test_lumped_projection_reproduces_a_constant(chamber):
    """Projecting a constant function gives that constant at every node."""
    V1 = fem.functionspace(chamber.create_mesh(8).mesh, ("Lagrange", 1))
    values = lumped_projection(V1, lambda r, z: np.full_like(r, 3.0))
    assert np.allclose(values, 3.0)


def test_lumped_projection_conserves_the_integral_of_a_thin_layer(chamber):
    """A layer thinner than the cells keeps its volume integral within 5% after
    projection."""
    import ufl
    from scifem import assemble_scalar

    mesh = chamber.create_mesh(8)
    V1 = fem.functionspace(mesh.mesh, ("Lagrange", 1))
    thickness = 5e-3
    R = chamber.can_radius

    def layer(r, z):
        return (r > R - thickness).astype(float)

    f = fem.Function(V1)
    f.x.array[:] = lumped_projection(V1, layer, subdivisions=40)
    weight = mesh.coordinate_system.integration_weight(mesh.mesh)
    integral = assemble_scalar(weight * f * ufl.dx)
    exact = np.pi * (R**2 - (R - thickness) ** 2) * chamber.length
    assert np.isclose(integral, exact, rtol=0.05)
