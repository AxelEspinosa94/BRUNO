import numpy as np
import pytest

from core.exceptions import ConstructionError
from core.ode_ivp.ode_ivp import IVP


def make(**overrides):
    data = {
        "mode": "function",
        "calculation_mode": "euler",
        "function": "y - x**2 + 1",
        "bounds": [0, 2],
        "y0": 0.5,
        "step": {"mode": "fixed", "n": 10},
    }
    data.update(overrides)
    return data


# ---------------------------------------------------------------- valid inputs
def test_scalar_with_n():
    ivp = IVP(make())
    assert (ivp.a, ivp.b, ivp.m, ivp.N) == (0.0, 2.0, 1, 10)
    assert ivp.h == pytest.approx(0.2)
    assert ivp.x[0] == 0.0 and ivp.x[-1] == 2.0 and len(ivp.x) == 11
    assert np.allclose(ivp.y0, [0.5])
    assert not ivp.is_system
    assert np.allclose(ivp.F(0.0, ivp.y0), [1.5])


def test_scalar_with_h_gives_same_mesh():
    a = IVP(make())
    b = IVP(make(step={"mode": "fixed", "h": 0.2}))
    assert b.N == a.N and np.allclose(a.x, b.x)


def test_h_with_float_noise_is_accepted():
    ivp = IVP(make(bounds=[0, 1], step={"mode": "fixed", "h": 0.1}))
    assert ivp.N == 10 and ivp.x[-1] == 1.0


def test_system():
    ivp = IVP(
        make(
            calculation_mode="rk4_system",
            function=["-4*y0 + 3*y1 + 6", "-2.4*y0 + 1.6*y1 + 3.6"],
            bounds=[0, 0.5],
            y0=[0, 0],
            step={"mode": "fixed", "h": 0.1},
        )
    )
    assert ivp.is_system and ivp.m == 2 and ivp.N == 5
    assert ivp.symbols == ["x", "y0", "y1"]
    assert np.allclose(ivp.F(0.0, ivp.y0), [6.0, 3.6])


def test_single_equation_in_list_is_scalar():
    ivp = IVP(make(function=["y - x**2 + 1"], y0=[0.5]))
    assert ivp.m == 1 and ivp.symbols == ["x", "y"]


def test_adaptive_step_has_no_mesh():
    ivp = IVP(
        make(
            step={"mode": "adaptive", "tol": 1e-6, "h0": 0.1, "hmin": 1e-6, "hmax": 0.5}
        )
    )
    assert ivp.step_mode == "adaptive" and ivp.N is None and ivp.x is None


# ---------------------------------------------------------------- invalid inputs
@pytest.mark.parametrize(
    "overrides, match",
    [
        ({"mode": "table"}, "mode='function'"),
        ({"function": None}, "'function'"),
        ({"function": ""}, "non-empty"),
        ({"function": []}, "'function'"),
        ({"function": "y + z"}, "unknown symbols"),
        ({"function": "y0 + x"}, "unknown symbols"),
        ({"bounds": None}, "Bounds"),
        ({"bounds": [2, 0]}, "Invalid interval"),
        ({"bounds": [[0, 1], [0, 1]]}, "single interval"),
        ({"y0": None}, "required"),
        ({"y0": "0.5"}, "finite number"),
        ({"y0": [0.5, 1.0]}, "1 equation"),
        ({"y0": float("nan")}, "finite number"),
        ({"step": None}, "Missing 'step'"),
        ({"step": {"mode": "fixed"}}, "exactly one of"),
        ({"step": {"mode": "fixed", "n": 0}}, "positive integer"),
        ({"step": {"mode": "fixed", "n": 2.5}}, "positive integer"),
        ({"step": {"mode": "fixed", "n": True}}, "positive integer"),
        ({"step": {"mode": "fixed", "h": -0.1}}, "positive number"),
        ({"step": {"mode": "fixed", "h": 3}}, "larger than the interval"),
        ({"step": {"mode": "fixed", "h": 0.3}}, "not an integer"),
        (
            {
                "step": {
                    "mode": "adaptive",
                    "tol": 1e-6,
                    "h0": 1,
                    "hmin": 0.1,
                    "hmax": 0.5,
                }
            },
            "hmin <= h0 <= hmax",
        ),
        (
            {
                "step": {
                    "mode": "adaptive",
                    "tol": 0,
                    "h0": 0.1,
                    "hmin": 0.01,
                    "hmax": 0.5,
                }
            },
            "'tol'",
        ),
    ],
)
def test_invalid_inputs(overrides, match):
    with pytest.raises(ConstructionError, match=match):
        IVP(make(**overrides))


def test_system_y0_length_mismatch():
    with pytest.raises(ConstructionError, match="2 equation"):
        IVP(make(function=["y1", "-y0"], y0=[1.0]))
