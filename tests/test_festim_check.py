from types import SimpleNamespace

import festim as F
import pytest

from ion_chambers._festim_check import REQUIRED_FEATURES, check_festim


def test_installed_festim_passes_the_check():
    """The FESTIM installed in the environment has every feature the package needs."""
    assert check_festim(F) is None


@pytest.mark.parametrize("feature", REQUIRED_FEATURES)
def test_festim_without_a_required_feature_is_rejected(feature):
    """A FESTIM missing one required class raises an ImportError that names it."""
    others = {name: object for name in REQUIRED_FEATURES if name != feature}
    old_festim = SimpleNamespace(__version__="2.1", **others)
    with pytest.raises(ImportError, match=f"found FESTIM 2.1, which lacks {feature}"):
        check_festim(old_festim)
