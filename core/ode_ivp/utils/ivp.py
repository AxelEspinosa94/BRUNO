import math
from typing import Any, Dict, Optional, Tuple

import numpy as np

from core.exceptions import ConstructionError

IVP_STEP_RULES = {
    "fixed": {
        "required": [],
        "exactly_one_of": ["n", "h"],
    },
    "adaptive": {
        "required": ["tol", "h0", "hmin", "hmax"],
        "exactly_one_of": [],
    },
}


def _validate_ivp_step(step: Optional[Dict[str, Any]] = None) -> None:
    """Checks that `step` has a valid mode and the fields that mode requires."""
    if not isinstance(step, dict):
        raise ConstructionError("Missing 'step' for IVP.")

    mode = step.get("mode")
    rules = IVP_STEP_RULES.get(mode)
    if rules is None:
        valid = ", ".join(f"'{m}'" for m in IVP_STEP_RULES)
        raise ConstructionError(
            f"Invalid mode '{mode}' for IVP. Must be one of: {valid}."
        )

    missing = [f for f in rules["required"] if step.get(f) is None]
    if missing:
        raise ConstructionError(f"Missing {missing} for IVP '{mode}' mode.")

    group = rules["exactly_one_of"]
    if group:
        given = [f for f in group if step.get(f) is not None]
        if len(given) != 1:
            raise ConstructionError(
                f"IVP '{mode}' mode requires exactly one of {group};"
                f" got {given or 'none'}."
            )


def _is_real(value: Any) -> bool:
    """
    True for finite int/float (numpy included), 
    False for bool, None, str, nan, inf.
    """
    if isinstance(value, bool) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        return False
    return math.isfinite(value)


def _normalize_fixed_step(
    a: float, b: float, step: Dict[str, Any], rtol: float = 1e-9
) -> Tuple[int, float]:
    """Returns (N, h) with h = (b - a) / N, from either step['n'] or step['h']."""
    n = step.get("n")
    if n is not None:
        if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 1:
            raise ConstructionError("n must be a positive integer.")
        return int(n), (b - a) / n

    h = step.get("h")
    if not _is_real(h) or h <= 0:
        raise ConstructionError("h must be a positive number.")

    if h > (b - a) * (1 + rtol):
        raise ConstructionError(f"h = {h} is larger than the interval length {b - a}.")

    ratio = (b - a) / h
    N = round(ratio)
    if abs(ratio - N) > rtol * max(1.0, ratio):
        raise ConstructionError(
            f"(b - a) / h = {ratio:.6g} is not an integer; "
            f"try h = {(b - a) / N:.6g} (n = {N})."
        )
    return N, (b - a) / N  # recompute h to remove rounding error


def _validate_adaptive_step(step: Dict[str, Any]) -> None:
    """tol, h0, hmin, hmax must be positive numbers with hmin <= h0 <= hmax."""
    for field in ("tol", "h0", "hmin", "hmax"):
        if not _is_real(step.get(field)) or step[field] <= 0:
            raise ConstructionError(f"'{field}' must be a positive number.")
    if not step["hmin"] <= step["h0"] <= step["hmax"]:
        raise ConstructionError("Adaptive step requires hmin <= h0 <= hmax.")


def _validate_y0(y0: Any, m: int) -> np.ndarray:
    """Returns y0 as a float array of length m (one initial value per equation)."""
    if y0 is None:
        raise ConstructionError("Initial condition 'y0' is required.")

    values = list(y0) if isinstance(y0, (list, tuple, np.ndarray)) else [y0]
    if not all(_is_real(v) for v in values):
        raise ConstructionError("Every value in 'y0' must be a finite number.")
    if len(values) != m:
        raise ConstructionError(
            f"Got {m} equation(s) but {len(values)} initial value(s) in 'y0'."
        )
    return np.asarray(values, dtype=float)
