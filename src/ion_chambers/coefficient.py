"""Transport and rate coefficients that may depend on the electric field."""

from collections.abc import Callable

import numpy as np

Coefficient = float | Callable[[np.ndarray], np.ndarray]
"""A transport or rate coefficient.

Either a constant, or a callable taking the reduced electric field E/N in townsend
(1 Td = 1e-21 V m2) as a numpy array and returning the coefficient in SI units.
"""


def check_coefficient(name: str, value: object) -> None:
    """Raises a TypeError unless ``value`` is a number or a callable.

    Args:
        name: the name of the attribute, used in the error message
        value: the value to check
    """
    is_number = isinstance(value, int | float) and not isinstance(value, bool)
    if not is_number and not callable(value):
        raise TypeError(
            f"{name} must be a float or a callable of the reduced field, "
            f"not {type(value).__name__}"
        )
