"""
Executor tests for ode_ivp.

Reference values come from the worked examples in docs/theory/ode_ivp/
(y' = y - x**2 + 1, y(0) = 0.5 on [0, 2], h = 0.2; exact y = (x+1)**2 - 0.5 e^x).
"""

import numpy as np
import pandas as pd
import pytest

from core.exceptions import ExecutionError
from core.ode_ivp.ode_ivp import IVP
from strategies.executors.ode_ivp.methods.tableaus import TABLEAUS
from strategies.executors.ode_ivp.ode_ivp_executors import IVPExecutor


def exact(x):
    return (x + 1) ** 2 - 0.5 * np.exp(x)


def make(method="euler", **overrides):
    data = {
        "mode": "function",
        "calculation_mode": method,
        "function": "y - x**2 + 1",
        "bounds": [0, 2],
        "y0": 0.5,
        "step": {"mode": "fixed", "n": 10},
    }
    data.update(overrides)
    return data


def run(method="euler", **overrides):
    return IVPExecutor().run(IVP(make(method, **overrides)))


# ---------------------------------------------------------------- theory examples
@pytest.mark.parametrize(
    "method, w_at_2",
    [
        ("euler", 4.8657845),
        ("heun", 5.2330546),
        ("rk2", 5.2903695),
        ("rk4", 5.3053630),
        ("adams_bashforth_2", 5.3992045),
        ("adams_bashforth_3", 5.3195640),
        ("adams_moulton_2", 5.3048287),
    ],
)
def test_matches_theory_examples(method, w_at_2):
    out = run(method)
    assert out["value"] == pytest.approx(w_at_2, abs=5e-7)
    assert len(out["x"]) == len(out["y"]) == 11


def test_first_steps_match_theory():
    assert run("euler")["y"][1] == pytest.approx(0.8)
    assert run("heun")["y"][1] == pytest.approx(0.826)
    assert run("rk2")["y"][1] == pytest.approx(0.828)
    assert run("adams_moulton_2")["y"][3] == pytest.approx(1.6489088, abs=5e-7)


def test_rk4_system_matches_theory_example():
    out = run(
        "rk4_system",
        function=["-4*y0 + 3*y1 + 6", "-2.4*y0 + 1.6*y1 + 3.6"],
        bounds=[0, 0.5],
        y0=[0, 0],
        step={"mode": "fixed", "h": 0.1},
    )
    assert np.allclose(out["y"][1], [0.53825520, 0.31962624], atol=5e-9)
    assert np.allclose(out["solution"], [1.79350749, 1.01440242], atol=5e-9)
    assert out["components"] == ["y0", "y1"] and out["plot_type"] == "multi_curve"


# ---------------------------------------------------------------- order of accuracy
@pytest.mark.parametrize(
    "method, order, ns",
    [
        ("euler", 1, (40, 80)),
        ("heun", 2, (40, 80)),
        ("rk2", 2, (40, 80)),
        ("rk4", 4, (20, 40)),
        ("adams_bashforth_2", 2, (160, 320)),
        ("adams_bashforth_3", 3, (160, 320)),
        ("adams_moulton_2", 3, (320, 640)),
    ],
)
def test_convergence_order(method, order, ns):
    errors = [
        abs(run(method, step={"mode": "fixed", "n": n})["value"] - exact(2.0))
        for n in ns
    ]
    observed = np.log2(errors[0] / errors[1])
    assert observed == pytest.approx(order, abs=0.2)


# ---------------------------------------------------------------- payload
def test_payload_shape_scalar():
    out = run("heun")
    for key in ("calculation_mode", "expression", "table", "x", "y", "value", "N", "h"):
        assert key in out
    table = out["table"]
    assert isinstance(table, pd.DataFrame)
    assert list(table.columns) == ["i", "x", "y", "k1", "k2"]
    assert table["k1"].iloc[0] is None or pd.isna(table["k1"].iloc[0])
    assert table["k1"].iloc[1] == pytest.approx(1.5)
    assert table["k2"].iloc[1] == pytest.approx(1.76)


def test_adams_tables_mark_starter_rows():
    table = run("adams_bashforth_3")["table"]
    assert list(table["source"][:4]) == ["y0", "rk4", "rk4", "ab3"]


def test_am2_reports_prediction_and_milne_estimate():
    table = run("adams_moulton_2")["table"]
    assert table["predicted"].iloc[3] == pytest.approx(1.6493272, abs=5e-7)
    assert table["milne_estimate"].iloc[3] == pytest.approx(4.18e-5, rel=0.01)


def test_constant_rhs_is_exact_for_every_method():
    for method in (
        "euler",
        "heun",
        "rk2",
        "rk4",
        "adams_bashforth_2",
        "adams_bashforth_3",
        "adams_moulton_2",
    ):
        out = run(method, function="2", y0=1.0)
        assert out["value"] == pytest.approx(5.0), method


# ---------------------------------------------------------------- errors
def test_divergence_is_reported():
    # y' = -100 y with h = 0.1: z = h*lambda = -10 is outside Euler's region (R = -9)
    with (
        pytest.raises(ExecutionError, match="diverged"),
        np.errstate(over="ignore", invalid="ignore"),
    ):
        run(
            "euler",
            function="-100*y",
            y0=1.0,
            bounds=[0, 40],
            step={"mode": "fixed", "n": 400},
        )


def test_adaptive_step_rejected_by_fixed_step_executors():
    ivp = IVP(
        make(
            "rk4",
            step={
                "mode": "adaptive",
                "tol": 1e-6,
                "h0": 0.1,
                "hmin": 1e-6,
                "hmax": 0.5,
            },
        )
    )
    with pytest.raises(ExecutionError, match="fixed step"):
        IVPExecutor().run(ivp)


def test_unsupported_method():
    ivp = IVP(make())
    ivp.calculation_mode = "shooting"
    with pytest.raises(ExecutionError, match="Not supported"):
        IVPExecutor().run(ivp)


# ---------------------------------------------------------------- catalog & tableaus
def test_every_executor_path_resolves():
    from app.utils.rule_loader import load_rule

    catalog = IVPExecutor().catalog
    for method in catalog["SUPPORTED_METHODS"]:
        assert callable(load_rule(catalog[method]["executor"])), method


@pytest.mark.parametrize("name", sorted(TABLEAUS))
def test_tableau_consistency(name):
    t = TABLEAUS[name]
    assert t.b.sum() == pytest.approx(1.0)
    assert np.allclose(np.triu(t.A), 0.0)  # explicit
    if t.order >= 2:
        assert t.b @ t.c == pytest.approx(0.5)


# ---------------------------------------------------------------- full pipeline
def test_full_pipeline_through_numerical_method():
    from core.base_method import NumericalMethod

    nm = NumericalMethod("ode_ivp", make("rk4"))
    assert nm.validate_input() is True
    out = nm.execute()
    assert out["status"] == "success"
    assert out["result"]["value"] == pytest.approx(5.3053630, abs=5e-7)
