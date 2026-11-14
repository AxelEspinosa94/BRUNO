import os

from dash import ALL, Input, Output, State, ctx, no_update

from app.layout.ode_ivp.components.ivp_input_area import equation_row
from app.utils.catalog_loader import load_catalog
from core.base_method import NumericalMethod
from core.contract import UIContract

_CATALOG = load_catalog(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "layout",
        "ode_ivp",
        "ode_ivp_layout_catalog.json",
    )
)
SYSTEM_METHODS = set(_CATALOG.get("system_methods", []))
MIN_EQ, MAX_EQ = 1, int(_CATALOG.get("max_equations", 6))

contract = UIContract()


# ============================================================
# Helpers (module level so they can be unit tested)
# ============================================================
def supports_systems(method) -> bool:
    return method in SYSTEM_METHODS


def build_rows(m, exprs=None, y0s=None):
    """Rebuilds the m equation rows keeping the values already typed."""
    exprs = list(exprs or [])
    y0s = list(y0s or [])
    return [
        equation_row(
            i,
            m,
            expr=exprs[i] if i < len(exprs) else None,
            y0=y0s[i] if i < len(y0s) else None,
        )
        for i in range(m)
    ]


def build_input_data(method, exprs, y0s, a, b, step_by, n, h):
    """
    Maps the form to the ode_ivp input_data contract:
    a single equation is sent as a scalar IVP (str / number), several as a system (lists).
    """
    exprs = list(exprs or [])
    y0s = list(y0s or [])
    scalar = len(exprs) == 1

    step = {"mode": "fixed"}
    if step_by == "h":
        step["h"] = h
    else:
        step["n"] = int(n) if isinstance(n, float) and n.is_integer() else n

    return {
        "mode": "function",
        "calculation_mode": method,
        "function": exprs[0] if scalar else exprs,
        "y0": y0s[0] if scalar else y0s,
        "bounds": [a, b],
        "step": step,
    }


def run_ivp(input_data):
    """Construct → validate → execute, always returning a UIContract view."""
    method = input_data.get("calculation_mode")
    try:
        nm = NumericalMethod("ode_ivp", input_data)
        nm.validate_input()
    except Exception as e:
        return contract.resolve(
            method,
            {
                "status": "error",
                "error_type": type(e).__name__,
                "message": str(e),
                "context": input_data,
            },
        )
    return contract.resolve(method, nm.execute())


# ============================================================
# Callbacks
# ============================================================
def register_ode_ivp_callbacks(app):

    # ------------------------------------------------------------
    # Callback 1: muestra el formulario y el botón al elegir método
    # ------------------------------------------------------------
    @app.callback(
        Output("ivp-mode-function", "hidden"),
        Output("ivp-btn-card", "hidden"),
        Input("ivp-method", "value"),
        Input("ivp-input-mode", "value"),
        prevent_initial_call=True,
    )
    def build_input_area(method, mode):
        if not method:
            return True, True
        return mode != "function", False

    # ------------------------------------------------------------
    # Callback 2: agrega / quita ecuaciones (solo métodos con sistemas)
    # ------------------------------------------------------------
    @app.callback(
        Output("ivp-equation-rows", "children"),
        Output("ivp-m", "data"),
        Output("ivp-add-eq", "disabled"),
        Output("ivp-remove-eq", "disabled"),
        Input("ivp-add-eq", "n_clicks"),
        Input("ivp-remove-eq", "n_clicks"),
        Input("ivp-method", "value"),
        State("ivp-m", "data"),
        State({"type": "ivp-fn", "index": ALL}, "value"),
        State({"type": "ivp-y0", "index": ALL}, "value"),
        prevent_initial_call=True,
    )
    def change_equations(_add, _remove, method, m, exprs, y0s):
        systems = supports_systems(method)
        new_m = m

        if ctx.triggered_id == "ivp-add-eq" and systems and m < MAX_EQ:
            new_m = m + 1
        elif ctx.triggered_id == "ivp-remove-eq" and m > MIN_EQ:
            new_m = m - 1
        elif ctx.triggered_id == "ivp-method" and not systems and m > MIN_EQ:
            new_m = MIN_EQ  # the new method only solves scalar IVPs

        add_disabled = not systems or new_m >= MAX_EQ
        remove_disabled = new_m <= MIN_EQ

        # Rows only need rebuilding when m changes (labels switch between y and y0, y1, ...)
        rows = build_rows(new_m, exprs, y0s) if new_m != m else no_update
        return rows, new_m, add_disabled, remove_disabled

    # ------------------------------------------------------------
    # Callback 3: paso definido por n o por h
    # ------------------------------------------------------------
    @app.callback(
        Output("ivp-n-wrap", "hidden"),
        Output("ivp-h-wrap", "hidden"),
        Input("ivp-step-by", "value"),
        prevent_initial_call=True,
    )
    def toggle_step_input(step_by):
        return step_by == "h", step_by != "h"

    # ------------------------------------------------------------
    # Callback 4: ejecuta el cálculo
    # ------------------------------------------------------------
    @app.callback(
        Output("ivp-result-area", "children"),
        Input("ivp-run-btn", "n_clicks"),
        State("ivp-method", "value"),
        State({"type": "ivp-fn", "index": ALL}, "value"),
        State({"type": "ivp-y0", "index": ALL}, "value"),
        State("ivp-a", "value"),
        State("ivp-b", "value"),
        State("ivp-step-by", "value"),
        State("ivp-n", "value"),
        State("ivp-h", "value"),
        prevent_initial_call=True,
    )
    def run_ode_ivp(n_clicks, method, exprs, y0s, a, b, step_by, n, h):
        if not method:
            return no_update
        input_data = build_input_data(method, exprs, y0s, a, b, step_by, n, h)
        return run_ivp(input_data)
