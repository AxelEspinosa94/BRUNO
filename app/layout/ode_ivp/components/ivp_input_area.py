from dash import dcc, html

from app.tooltips import get_tooltip
from core.ui.styled_components import styled_input, styled_radioitems
from core.ui.tooltips import Tooltip


# ───────────────────────────────────────────────
# Fila de ecuación: y_i' = f_i(x, y)   |   y_i(a) = valor
# ───────────────────────────────────────────────
def equation_row(i: int, m: int, expr=None, y0=None):
    """Row i of an m-equation IVP. Scalar IVPs (m = 1) use the symbol y."""
    name = "y" if m == 1 else f"y{i}"
    return html.Div(
        className="input-row",
        children=[
            html.Div(f"{name}' =", className="na-label bound-level"),
            styled_input(
                id={"type": "ivp-fn", "index": i},
                type="text",
                placeholder=(
                    "ex: y - x**2 + 1" if m == 1 else f"f{i}(x, y0, ..., y{m - 1})"
                ),
                value=expr,
            ),
            html.Div(f"{name}(a) =", className="na-label bound-level"),
            styled_input(
                id={"type": "ivp-y0", "index": i},
                type="number",
                placeholder="y0",
                value=y0,
            ),
        ],
    )


def _label(text, tooltip_key):
    return html.Div(
        className="label-with-tooltip",
        children=[
            html.Div(text, className="na-label"),
            Tooltip(get_tooltip(tooltip_key)).render(),
        ],
    )


# ───────────────────────────────────────────────
# Área dinámica de input (contenedores predefinidos)
# ───────────────────────────────────────────────
def ivp_input_area():
    return html.Div(
        id="ivp-input-area",
        className="module-card input-area",
        children=[
            html.Div(
                id="ivp-mode-function",
                hidden=True,
                children=[
                    # Ecuaciones + condiciones iniciales
                    _label("Equations and initial values", "ivp-fn"),
                    html.Div(id="ivp-equation-rows", children=[equation_row(0, 1)]),
                    html.Div(
                        className="input-row",
                        children=[
                            html.Button(
                                "+ equation",
                                className="btn btn-secondary",
                                id="ivp-add-eq",
                                n_clicks=0,
                                disabled=True,
                            ),
                            html.Button(
                                "− equation",
                                className="btn btn-secondary",
                                id="ivp-remove-eq",
                                n_clicks=0,
                                disabled=True,
                            ),
                        ],
                    ),
                    dcc.Store(id="ivp-m", data=1),
                    # Intervalo [a, b]
                    _label("Interval [a, b]", "ivp-interval"),
                    html.Div(
                        className="input-row",
                        children=[
                            html.Div("x", className="na-label bound-level"),
                            styled_input(id="ivp-a", type="number", placeholder="a"),
                            styled_input(id="ivp-b", type="number", placeholder="b"),
                        ],
                    ),
                    # Paso: por n o por h
                    _label("Step definition", "ivp-step-by"),
                    styled_radioitems(
                        id="ivp-step-by",
                        options=[
                            {"label": "Number of steps (n)", "value": "n"},
                            {"label": "Step size (h)", "value": "h"},
                        ],
                        value="n",
                    ),
                    html.Div(
                        id="ivp-n-wrap",
                        children=[
                            _label("Number of steps (n)", "ivp-n"),
                            styled_input(
                                id="ivp-n", type="number", placeholder="ej: 10"
                            ),
                        ],
                    ),
                    html.Div(
                        id="ivp-h-wrap",
                        hidden=True,
                        children=[
                            _label("Step size (h)", "ivp-h"),
                            styled_input(
                                id="ivp-h", type="number", placeholder="ej: 0.1"
                            ),
                        ],
                    ),
                ],
            )
        ],
    )
