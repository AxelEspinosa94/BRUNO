from core.exceptions import ValidationError
from strategies.validators.ode_ivp.rules.common_rules import fixed_steps


def validate_min_steps(input_data, params):
    """
    A k-step method needs k - 1 starting values plus at least one step of its own,
    so the mesh must have N >= min_steps.
    """
    k = params["min_steps"]
    N = fixed_steps(input_data)
    if N < k:
        method = input_data.get("calculation_mode")
        raise ValidationError(
            f"Method '{method}' needs at least {k} steps; got N = {N}. "
            "Use a larger n or a smaller h."
        )
