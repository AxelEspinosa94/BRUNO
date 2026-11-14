from dash import html

from app.tooltips import get_tooltip
from core.ui.styled_components import styled_button
from core.ui.tooltips import Tooltip


# ───────────────────────────────────────────────
# Botón de ejecución
# ───────────────────────────────────────────────
def ivp_run_btn():
    return html.Div(
        className="module-card",
        id="ivp-btn-card",
        hidden=True,
        children=[
            Tooltip(get_tooltip("ivp-run-btn")).render(),
            styled_button("ivp-run-btn", "Calcular"),
        ],
    )
