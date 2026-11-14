from typing import Dict, List

from dash import html

from app.tooltips import get_tooltip
from core.ui.styled_components import styled_radioitems
from core.ui.tooltips import Tooltip


# ───────────────────────────────────────────────
# Modo de entrada (siempre función)
# ───────────────────────────────────────────────
def ivp_input_mode(options: List[Dict] = None):
    return html.Div(
        className="module-card",
        id="ivp-mode-card",
        children=[
            html.Div(
                className="label-with-tooltip",
                children=[
                    html.Div("Input Type", className="na-label"),
                    Tooltip(get_tooltip("ivp-input-mode")).render(),
                ],
            ),
            styled_radioitems(
                id="ivp-input-mode",
                options=options,
                value="function",
            ),
        ],
    )
