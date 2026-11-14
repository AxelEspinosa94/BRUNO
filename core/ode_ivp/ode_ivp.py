import os

import numpy as np

from app.utils.bounds_validator import validate_bounds
from app.utils.build_function import build_ode_function
from app.utils.catalog_loader import load_catalog
from core.exceptions import ConstructionError
from core.ode_ivp.utils.ivp import (
    _normalize_fixed_step,
    _validate_adaptive_step,
    _validate_ivp_step,
    _validate_y0,
)


class IVP:
    """
    Constructor for initial value problems y' = F(x, y), y(a) = y0, x in [a, b].
    Performs structural validation and stores parsed fields.
    Method-level validation (e.g. N >= k for multistep) is done in the validator layer.

    Input:
        function : str for a scalar IVP (symbols x, y) or list[str] for a system
                   of m equations (symbols x, y0, ..., y{m-1}).
        bounds   : [a, b] with a < b.
        y0       : number (scalar) or list of m numbers (system).
        step     : {"mode": "fixed", "n": int} or {"mode": "fixed", "h": float}
                   or {"mode": "adaptive", "tol", "h0", "hmin", "hmax"}.
    """

    def __init__(self, input_data):
        base = os.path.dirname(__file__)
        path = os.path.join(base, "ode_ivp_constructor_catalog.json")
        self.catalog = load_catalog(path)

        self.input_data = input_data
        self.mode = self.input_data.get("mode")
        # IVP calculation modes: euler, heun, rk2, rk4, rk4_system,
        # adams_bashforth_2, adams_bashforth_3, adams_moulton_2.
        self.calculation_mode = self.input_data.get("calculation_mode")

        supported = self.catalog.get("SUPPORTED_METHODS", [])
        if self.calculation_mode not in supported:
            raise ConstructionError(
                f"Not supported method for IVP: {self.calculation_mode}. "
                f"Must be one of: {', '.join(supported)}."
            )

        if self.mode != "function":
            raise ConstructionError("ODE IVP only supports mode='function'.")

        # ------------------------------------------------------------
        # 1. Function(s): scalar IVP or system
        # ------------------------------------------------------------
        funcs = self.input_data.get("function")
        if isinstance(funcs, str):
            self.func_str = [funcs]
        elif isinstance(funcs, (list, tuple)) and len(funcs) > 0:
            self.func_str = list(funcs)
        else:
            raise ConstructionError(
                "IVP requires 'function' as a string or a non-empty list of strings."
            )

        self.m = len(self.func_str)
        self.is_system = self.m > 1
        self.F, self.symbols = build_ode_function(self.func_str)

        # ------------------------------------------------------------
        # 2. Interval (always 1D)
        # ------------------------------------------------------------
        self.bounds = validate_bounds(self.input_data.get("bounds"))
        if len(self.bounds) != 1:
            raise ConstructionError("IVP requires a single interval [a, b].")
        self.a, self.b = (float(v) for v in self.bounds[0])

        # ------------------------------------------------------------
        # 3. Initial condition: one value per equation
        # ------------------------------------------------------------
        self.y0 = _validate_y0(self.input_data.get("y0"), self.m)

        # ------------------------------------------------------------
        # 4. Step: normalized to a canonical form
        # ------------------------------------------------------------
        self.step = self.input_data.get("step")
        _validate_ivp_step(self.step)
        self.step_mode = self.step["mode"]

        if self.step_mode == "fixed":
            self.N, self.h = _normalize_fixed_step(self.a, self.b, self.step)
            self.x = np.linspace(self.a, self.b, self.N + 1)
        else:
            _validate_adaptive_step(self.step)
            self.N, self.h, self.x = None, None, None  # mesh is an output here
