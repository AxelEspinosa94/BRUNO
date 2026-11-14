import numpy as np

from core.exceptions import ExecutionError


def require_fixed_mesh(instance):
    """Fixed-step executors need the mesh built by the constructor."""
    if getattr(instance, "step_mode", None) != "fixed" or instance.x is None:
        raise ExecutionError(
            f"Method '{instance.calculation_mode}'"
            " requires a fixed step (mode='fixed')."
        )
    return np.asarray(instance.x, dtype=float), float(instance.h)


def check_finite(w, x, method):
    """
    Stops as soon as the solution overflows
    (typically h outside the stability region).
    """
    if not np.all(np.isfinite(w)):
        raise ExecutionError(
            f"{method}: the solution diverged at x = {x:.6g}. "
            "The step size may be outside the method's stability region;"
            " try a smaller h."
        )
