"""Check that the installed FESTIM is recent enough for this package."""

REQUIRED_FEATURES = ("ElectromigrationTerm", "OutflowBC", "GenericReaction")
"""FESTIM classes this package builds on. All first appeared in FESTIM 2.2."""


def check_festim(festim) -> None:
    """Raises a clear error if FESTIM lacks a feature this package needs.

    The check looks for the features rather than comparing version numbers, so
    development builds of FESTIM (versions like ``2.2.dev0``) are judged correctly.

    Args:
        festim: the imported ``festim`` module
    """
    missing = [name for name in REQUIRED_FEATURES if not hasattr(festim, name)]
    if missing:
        version = getattr(festim, "__version__", "unknown")
        raise ImportError(
            f"ion_chambers needs FESTIM 2.2 or newer, but found FESTIM {version}, "
            f"which lacks {', '.join(missing)}. Update FESTIM, for example with "
            f"'conda install -c conda-forge festim=2.2'."
        )
