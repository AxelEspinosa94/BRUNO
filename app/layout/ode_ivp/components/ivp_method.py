from typing import Dict, List

from dash import html

from app.tooltips import get_tooltip
from core.ui.styled_components import styled_dropdown
from core.ui.tooltips import Tooltip


# ───────────────────────────────────────────────
# Selector de método
# ───────────────────────────────────────────────
def ivp_method(options: List[Dict] = None):
    return html.Div(
        className="module-card",
        children=[
            html.Div(
                className="label-with-tooltip",
                children=[
                    html.Div("Method", className="na-label"),
                    Tooltip(get_tooltip("ivp-method")).render(),
                ],
            ),
            styled_dropdown(
                id="ivp-method",
                options=options,
                placeholder="Selecciona un método",
            ),
        ],
    )
