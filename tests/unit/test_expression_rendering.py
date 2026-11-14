"""Expressions: LaTeX (expression_format="latex") goes through MathJax, plain text stays a code block."""

import pytest
from dash import dcc

from core.base_method import NumericalMethod
from core.contract import UIContract
from core.renderer import Renderer


def _markdown_blocks(component, found=None):
    found = [] if found is None else found
    if isinstance(component, dcc.Markdown):
        found.append(component)
    children = getattr(component, "children", None)
    if isinstance(children, (list, tuple)):
        for c in children:
            _markdown_blocks(c, found)
    elif children is not None and not isinstance(children, str):
        _markdown_blocks(children, found)
    return found


def _expression_md(view):
    return [
        m for m in _markdown_blocks(view) if "result-expression" in (m.className or "")
    ]


# ---------------------------------------------------------------- renderer
@pytest.mark.parametrize(
    "result, fmt",
    [
        ({"expression": "x**2 + 1"}, "text"),
        ({"expression": "$a^2$", "expression_format": "latex"}, "latex"),
        (
            {"markdown": "### Title", "expression_format": "latex"},
            "text",
        ),  # only applies to expression
    ],
)
def test_renderer_markdown_format(result, fmt):
    block = Renderer().render("any", result)["blocks"][0]
    assert block["type"] == "markdown" and block["format"] == fmt


# ---------------------------------------------------------------- contract
@pytest.mark.parametrize(
    "content", ["$w_{i+1} = w_i$", "$$w_{i+1} = w_i$$", "w_{i+1} = w_i"]
)
def test_contract_latex_uses_mathjax_display(content):
    block = {"type": "markdown", "content": content, "format": "latex"}
    md = _expression_md(UIContract()._block_expression(block))[0]
    assert md.mathjax is True
    assert md.children == "$$\nw_{i+1} = w_i\n$$"


def test_contract_text_keeps_code_block():
    block = {"type": "markdown", "content": "L(x) = 2*x**2 + 1"}
    md = _expression_md(UIContract()._block_expression(block))[0]
    assert md.children == "```\nL(x) = 2*x**2 + 1\n```"
    assert not getattr(md, "mathjax", False)


# ---------------------------------------------------------------- full pipeline
def _view(module, data):
    nm = NumericalMethod(module, data)
    nm.validate_input()
    return UIContract().resolve(data["calculation_mode"], nm.execute())


def test_ode_ivp_expression_is_rendered_as_latex():
    view = _view(
        "ode_ivp",
        {
            "mode": "function",
            "calculation_mode": "rk2",
            "function": "y - x**2 + 1",
            "bounds": [0, 2],
            "y0": 0.5,
            "step": {"mode": "fixed", "n": 10},
        },
    )
    md = _expression_md(view)[0]
    assert (
        md.mathjax is True
        and "\\tfrac{h}{2}" in md.children
        and "```" not in md.children
    )


def test_integration_expression_is_rendered_as_latex():
    view = _view(
        "integration",
        {
            "mode": "function",
            "calculation_mode": "simpson_1_3",
            "function": "x**2",
            "bounds": [[0, 1]],
            "n": 10,
        },
    )
    md = _expression_md(view)[0]
    assert md.mathjax is True
    assert md.children.startswith("$$\n\\int_{0}^{1} f(x)\\,dx \\approx 0.333333")
