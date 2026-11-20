import json

import numpy as np
import pytest
from dolfinx import fem

import ion_chambers as ic


@pytest.mark.parametrize("activity", ["high", None, True])
def test_activity_rejects_wrong_types(activity):
    """An activity that is not a number raises a TypeError."""
    with pytest.raises(TypeError, match="activity must be a float"):
        ic.TritiumInGas(activity=activity)


def test_activity_rejects_negative_values():
    """A negative activity raises a ValueError."""
    with pytest.raises(ValueError, match="activity must not be negative"):
        ic.TritiumInGas(activity=-1.0)


def test_ideal_pair_rate(chamber):
    """The ideal pair rate is the activity times the gas volume times the ion pairs
    per decay."""
    source = ic.TritiumInGas(activity=3.7e10)
    expected = 3.7e10 * 1e-3 * ic.beta.ion_pairs_per_decay(33.97)
    assert np.isclose(source.ideal_pair_rate(chamber, ic.library.air()), expected)


def test_no_wall_loss_gives_unit_efficiency(chamber):
    """With wall losses switched off the efficiency is exactly one everywhere."""
    source = ic.TritiumInGas(activity=1.0, wall_loss=False)
    source.build(chamber.build_mesh(6), chamber, ic.library.air())
    assert np.all(source.efficiency.x.array == 1.0)


def test_wall_loss_efficiency_stays_between_zero_and_one(chamber):
    """The projected wall-loss efficiency stays within its physical bounds of 0 and
    1."""
    source = ic.TritiumInGas(activity=1.0, n_rays=200)
    source.build(chamber.build_mesh(6), chamber, ic.library.air())
    values = source.efficiency.x.array
    assert values.min() >= 0.0 and values.max() <= 1.0


def test_changing_the_activity_updates_the_built_source(chamber):
    """Setting a new activity on a built source updates its rate constant without a
    rebuild."""
    source = ic.TritiumInGas(activity=1.0, wall_loss=False)
    source.build(chamber.build_mesh(6), chamber, ic.library.air())
    source.activity = 10.0
    assert np.isclose(source._rate.value, 10.0 * ic.beta.ion_pairs_per_decay(33.97))


def test_built_source_is_a_ufl_expression(chamber):
    """Building a source returns an expression and stores the efficiency as a DOLFINx
    function."""
    source = ic.TritiumInGas(activity=1.0, wall_loss=False)
    expression = source.build(chamber.build_mesh(6), chamber, ic.library.air())
    assert isinstance(source.efficiency, fem.Function) and expression is not None


@pytest.mark.parametrize(
    "coarse, fine, order, ratio, expected",
    [
        (1.2, 1.1, 1, 2.0, 1.0),  # error halves: 0.2 -> 0.1 -> 0
        (1.4, 1.1, 2, 2.0, 1.0),  # error quarters: 0.4 -> 0.1 -> 0
        (1.3, 1.1, 1, 3.0, 1.0),  # error divided by three
        (5.0, 5.0, 1, 2.0, 5.0),  # already converged
    ],
)
def test_richardson(coarse, fine, order, ratio, expected):
    """Richardson extrapolation recovers the exact value from two results whose error
    follows the stated order."""
    assert np.isclose(ic.richardson(coarse, fine, order=order, ratio=ratio), expected)


def test_save_json_creates_the_folder(tmp_path):
    """save_json creates a missing folder and writes numpy floats as plain numbers."""
    path = tmp_path / "new_folder" / "out.json"
    ic.save_json({"current": np.float64(1.5)}, path)
    assert json.loads(path.read_text()) == {"current": 1.5}
