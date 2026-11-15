from app.utils.bounds_validator import validate_bounds
from core.exceptions import ValidationError
from core.ode_ivp.utils.ivp import _normalize_fixed_step


def number_of_equations(input_data):
    """1 for a scalar IVP ('function' is a str), m for a system (list of m str)."""
    funcs = input_data.get("function")
    return 1 if isinstance(funcs, str) else len(funcs)


def fixed_steps(input_data):
    """N for a fixed step, computed exactly as the constructor does (from n or h)."""
    a, b = validate_bounds(input_data.get("bounds"))[0]
    N, _ = _normalize_fixed_step(float(a), float(b), input_data["step"])
    return N


def validate_step_mode(input_data, params):
    """The step mode ('fixed' / 'adaptive') must be one the method supports."""
    mode = input_data["step"]["mode"]
    allowed = params.get("step_modes", ["fixed"])
    if mode not in allowed:
        method = input_data.get("calculation_mode")
        raise ValidationError(
            f"Method '{method}' does not support step mode '{mode}'. "
            f"Supported: {', '.join(allowed)}."
        )


def validate_system_support(input_data, params):
    """Systems (m > 1) are only allowed for methods released with systems=true."""
    m = number_of_equations(input_data)
    if m > 1 and not params.get("systems", False):
        method = input_data.get("calculation_mode")
        raise ValidationError(
            f"Method '{method}' does not support systems of equations (got {m})."
        )
