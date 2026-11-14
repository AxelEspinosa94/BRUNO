import numpy as np
import pytest

from app.utils.build_function import build_function, build_ode_function, ode_symbols
from core.exceptions import ConstructionError


# ---------------------------------------------------------------- build_function
def test_univariate():
    f = build_function("x**2 + 1")
    assert np.allclose(f(np.array([0.0, 2.0])), [1.0, 5.0])


def test_multivariate_symbol_order():
    f = build_function("x0 - 2*x1", ["x0", "x1"])
    assert f(3.0, 1.0) == pytest.approx(1.0)


def test_constant_keeps_shape():
    f = build_function("5", ["x0", "x1"])
    out = f(np.zeros((2, 3)), np.zeros((2, 3)))
    assert out.shape == (2, 3) and np.all(out == 5.0)


def test_unused_symbols_are_allowed():
    f = build_function("x0", ["x0", "x1"])
    assert f(2.0, 7.0) == pytest.approx(2.0)


@pytest.mark.parametrize(
    "expr, symbols",
    [("x + z", ["x"]), ("x0 + x2", ["x0", "x1"]), ("x0", ["x"]), ("y", None)],
)
def test_unknown_symbols_rejected(expr, symbols):
    with pytest.raises(ConstructionError, match="unknown symbols"):
        build_function(expr, symbols)


def test_euler_number_hint():
    with pytest.raises(ConstructionError, match="exp"):
        build_function("e**x")


@pytest.mark.parametrize("expr", ["", "   ", None, 3])
def test_empty_or_non_string_rejected(expr):
    with pytest.raises(ConstructionError):
        build_function(expr)


@pytest.mark.parametrize("expr", ["x +* 2", "x >= 1"])
def test_unparseable_or_non_expression_rejected(expr):
    with pytest.raises(ConstructionError):
        build_function(expr)


# ---------------------------------------------------------------- ODE helpers
def test_ode_symbols():
    assert ode_symbols(1) == ["x", "y"]
    assert ode_symbols(3) == ["x", "y0", "y1", "y2"]
    with pytest.raises(ConstructionError):
        ode_symbols(0)


def test_ode_scalar():
    F, names = build_ode_function("y - x**2 + 1")
    assert names == ["x", "y"]
    assert np.allclose(F(0.0, 0.5), [1.5])


def test_ode_system_matches_rk4_theory_example():
    F, names = build_ode_function(["-4*y0 + 3*y1 + 6", "-2.4*y0 + 1.6*y1 + 3.6"])
    assert names == ["x", "y0", "y1"]
    assert np.allclose(F(0.05, [0.3, 0.18]), [5.34, 3.168])


def test_ode_constant_component():
    F, _ = build_ode_function(["1", "y0"])
    assert np.allclose(F(0.0, [2.0, 3.0]), [1.0, 2.0])


def test_ode_error_names_equation():
    with pytest.raises(ConstructionError, match="Equation 1"):
        build_ode_function(["y0", "y2"])


def test_python_ternary_suggests_piecewise():
    with pytest.raises(ConstructionError, match="Piecewise"):
        build_function("1 if x < 1 else 0")


def test_piecewise_is_supported():
    f = build_function("Piecewise((1, x < 1), (0, True))")
    assert f(0.5) == pytest.approx(1.0) and f(1.5) == pytest.approx(0.0)
    assert np.allclose(f(np.array([0.0, 0.99, 1.0, 2.0])), [1.0, 1.0, 0.0, 0.0])
