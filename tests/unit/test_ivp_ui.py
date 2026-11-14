"""UI tests for the ode_ivp module: layout, form → input_data mapping and full pipeline."""

import json
from pathlib import Path

import pytest
from dash import html

from app.callbacks.ode_ivp.ode_ivp_callbacks import (
    SYSTEM_METHODS,
    build_input_data,
    build_rows,
    run_ivp,
    supports_systems,
)
from app.layout.ode_ivp.ode_ivp_layout import ode_ivp_section


def _ids(component, found=None):
    """Collects every component id in a Dash tree."""
    found = [] if found is None else found
    cid = getattr(component, "id", None)
    if cid is not None:
        found.append(json.dumps(cid, sort_keys=True) if isinstance(cid, dict) else cid)
    children = getattr(component, "children", None)
    if isinstance(children, (list, tuple)):
        for c in children:
            _ids(c, found)
    elif children is not None and hasattr(children, "children"):
        _ids(children, found)
    return found


# ---------------------------------------------------------------- layout
def test_layout_has_expected_ids_and_no_duplicates():
    ids = _ids(ode_ivp_section)
    for expected in (
        "ivp-container",
        "ivp-method",
        "ivp-input-mode",
        "ivp-mode-function",
        "ivp-equation-rows",
        "ivp-add-eq",
        "ivp-remove-eq",
        "ivp-m",
        "ivp-a",
        "ivp-b",
        "ivp-step-by",
        "ivp-n",
        "ivp-h",
        "ivp-btn-card",
        "ivp-run-btn",
        "ivp-result-area",
    ):
        assert expected in ids, expected
    assert len(ids) == len(set(ids)), "duplicated component ids"
    assert all(not isinstance(i, str) or not i.startswith("integr-") for i in ids)


def test_layout_methods_are_released_in_backend():
    layout = json.loads(
        Path("app/layout/ode_ivp/ode_ivp_layout_catalog.json").read_text(
            encoding="utf-8"
        )
    )
    ctor = json.loads(
        Path("core/ode_ivp/ode_ivp_constructor_catalog.json").read_text(
            encoding="utf-8"
        )
    )
    options = {o["value"] for o in layout["method_selector_options"]}
    assert options == set(layout["supported_methods"])
    assert options <= set(ctor["SUPPORTED_METHODS"])
    assert set(layout["system_methods"]) <= options


def test_rows_switch_labels_between_scalar_and_system():
    scalar = build_rows(1, ["y"], [1.0])
    system = build_rows(2, ["y1", "-y0"], [1.0, 0.0])
    assert scalar[0].children[0].children == "y' ="
    assert [r.children[0].children for r in system] == ["y0' =", "y1' ="]
    assert system[1].children[1].value == "-y0" and system[1].children[3].value == 0.0


def test_rows_keep_typed_values_when_growing():
    rows = build_rows(3, ["y1", "-y0"], [1.0, 0.0])
    assert rows[0].children[1].value == "y1" and rows[2].children[1].value is None


def test_system_support_comes_from_catalog():
    assert supports_systems("rk4_system") and not supports_systems("euler")
    assert SYSTEM_METHODS == {"rk4_system"}


# ---------------------------------------------------------------- form → input_data
def test_scalar_form_with_n():
    data = build_input_data("heun", ["y - x**2 + 1"], [0.5], 0, 2, "n", 10.0, None)
    assert data == {
        "mode": "function",
        "calculation_mode": "heun",
        "function": "y - x**2 + 1",
        "y0": 0.5,
        "bounds": [0, 2],
        "step": {"mode": "fixed", "n": 10},
    }
    assert isinstance(data["step"]["n"], int)


def test_system_form_with_h():
    data = build_input_data("rk4_system", ["y1", "-y0"], [1, 0], 0, 1, "h", None, 0.1)
    assert data["function"] == ["y1", "-y0"] and data["y0"] == [1, 0]
    assert data["step"] == {"mode": "fixed", "h": 0.1}


# ---------------------------------------------------------------- full pipeline
def _text(component):
    """Concatenates every string in a Dash tree (to look for messages)."""
    if isinstance(component, str):
        return component
    children = getattr(component, "children", None)
    if isinstance(children, (list, tuple)):
        return " ".join(_text(c) for c in children if c is not None)
    return _text(children) if children is not None else ""


def test_pipeline_success_scalar():
    view = run_ivp(
        build_input_data("rk4", ["y - x**2 + 1"], [0.5], 0, 2, "n", 10, None)
    )
    assert isinstance(view, html.Div) and view.className == "result-container"
    assert "5.30536" in _text(view)


def test_pipeline_success_system():
    view = run_ivp(
        build_input_data(
            "rk4_system",
            ["-4*y0 + 3*y1 + 6", "-2.4*y0 + 1.6*y1 + 3.6"],
            [0, 0],
            0,
            0.5,
            "h",
            None,
            0.1,
        )
    )
    assert view.className == "result-container"


@pytest.mark.parametrize(
    "args, message",
    [
        (("rk4", ["y + z"], [0.5], 0, 2, "n", 10, None), "unknown symbols"),
        (("rk4", ["1 if x < 1 else 0"], [0], 0, 2, "n", 10, None), "Piecewise"),
        (("rk4", ["y"], [None], 0, 2, "n", 10, None), "required"),
        (("rk4", ["y"], [1], 2, 0, "n", 10, None), "Invalid interval"),
        (("rk4", ["y"], [1], 0, 1, "h", None, 0.3), "not an integer"),
        (("euler", ["y1", "-y0"], [1, 0], 0, 1, "n", 10, None), "systems"),
        (("adams_bashforth_3", ["y"], [1], 0, 1, "n", 2, None), "at least 3"),
    ],
)
def test_pipeline_errors_are_shown(args, message):
    view = run_ivp(build_input_data(*args))
    assert message in _text(view)
