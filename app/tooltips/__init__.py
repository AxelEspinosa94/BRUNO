from ._global import TOOLTIPS as GLOBAL
from .integration import TOOLTIPS as INTEGRATION
from .linear_algebra import TOOLTIPS as LINEAR_ALGEBRA
from .nonlinear import TOOLTIPS as NONLINEAR
from .numerical_derivative import TOOLTIPS as NUMERICAL_DERIVATIVE
from .ode import TOOLTIPS as ODE
from .ode_ivp import TOOLTIPS as ODE_IVP

ALL_TOOLTIPS = {
    **GLOBAL,
    **INTEGRATION,
    **LINEAR_ALGEBRA,
    **NONLINEAR,
    **ODE,
    **ODE_IVP,
    **NUMERICAL_DERIVATIVE,
}


def get_tooltip(key):
    return ALL_TOOLTIPS.get(key, "")
