"""
Integration tests for the ode_ivp module (Stage 5).

Cover the full pipeline across layers:
UI form → input_data → NumericalMethod (constructor → validator → executor)
→ Renderer → UIContract, plus consistency between the per-layer catalogs.
"""

import json
from pathlib import Path

import numpy as np
import pytest
from dash import html

from app.callbacks.ode_ivp.ode_ivp_callbacks import build_input_data, run_ivp
from core.base_method import NumericalMethod
from core.contract import UIContract
from core.renderer import Renderer

contract = UIContract()

METHODS = [
    "euler",
    "heun",
    "rk2",
    "rk4",
    "rk4_system",
    "adams_bashforth_2",
    "adams_bashforth_3",
    "adams_moulton_2",
]
SCALAR_METHODS = [m for m in METHODS if m != "rk4_system"]

# Orden teórico de cada método (docs/theory/ode_ivp/)
ORDER = {
    "euler": 1,
    "heun": 2,
    "rk2": 2,
    "rk4": 4,
    "adams_bashforth_2": 2,
    "adams_bashforth_3": 3,
    "adams_moulton_2": 3,
}

CATALOGS = {
    "layout": "app/layout/ode_ivp/ode_ivp_layout_catalog.json",
    "constructor": "core/ode_ivp/ode_ivp_constructor_catalog.json",
    "validator": "strategies/validators/ode_ivp/ode_ivp_validation_catalog.json",
    "executor": "strategies/executors/ode_ivp/ode_ivp_executors_catalog.json",
}


def make_input(method, function="y - x**2 + 1", bounds=(0, 2), y0=0.5, **step):
    step = step or {"n": 10}
    return {
        "mode": "function",
        "calculation_mode": method,
        "function": function,
        "bounds": list(bounds),
        "y0": y0,
        "step": {"mode": "fixed", **step},
    }


def make_outcome(method, **kwargs):
    nm = NumericalMethod("ode_ivp", make_input(method, **kwargs))
    nm.validate_input()
    return nm.execute()


def exact(x):
    return (x + 1) ** 2 - 0.5 * np.exp(x)


# ────────────────────────────────────────────────────────────────
# Catálogos: todas las capas liberan los mismos métodos
# ────────────────────────────────────────────────────────────────


def _released(layer):
    data = json.loads(Path(CATALOGS[layer]).read_text(encoding="utf-8"))
    return set(data.get("SUPPORTED_METHODS") or data.get("supported_methods"))


@pytest.mark.parametrize("layer", ["layout", "constructor", "validator", "executor"])
def test_catalogos_liberan_mismos_metodos(layer):
    """Con el módulo completo, ningún método queda retenido en alguna capa."""
    assert _released(layer) == set(METHODS)


def test_registro_en_method_catalog():
    catalog = json.loads(Path("core/method_catalog.json").read_text(encoding="utf-8"))
    entry = catalog["ode_ivp"]
    assert entry["classConstructor"] == "core.ode_ivp.ode_ivp.IVP"
    assert entry["classInputValidator"].endswith("IVPValidator")
    assert entry["classExecutor"].endswith("IVPExecutor")


# ────────────────────────────────────────────────────────────────
# Pipeline: estructura del resultado
# ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", SCALAR_METHODS)
def test_estructura_resultado_escalar(method):
    outcome = make_outcome(method)
    assert outcome["status"] == "success"
    result = outcome["result"]
    for key in (
        "calculation_mode",
        "expression",
        "expression_format",
        "table",
        "x",
        "y",
        "value",
        "a",
        "b",
        "N",
        "h",
        "components",
        "plot_type",
    ):
        assert key in result, key
    assert result["calculation_mode"] == method
    assert result["N"] == 10 and result["h"] == pytest.approx(0.2)
    assert len(result["x"]) == len(result["y"]) == len(result["table"]) == 11
    assert result["x"][0] == 0.0 and result["x"][-1] == 2.0
    assert result["y"][0] == 0.5
    assert result["plot_type"] == "curve" and result["components"] == ["y"]


def test_estructura_resultado_sistema():
    outcome = make_outcome(
        "rk4_system", function=["y1", "-y0"], y0=[1.0, 0.0], bounds=(0, np.pi), n=50
    )
    assert outcome["status"] == "success"
    result = outcome["result"]
    assert result["components"] == ["y0", "y1"] and result["plot_type"] == "multi_curve"
    assert np.array(result["y"]).shape == (51, 2)
    assert len(result["solution"]) == 2
    assert "value" not in result
    assert list(result["table"].columns[:4]) == ["i", "x", "y0", "y1"]


# ────────────────────────────────────────────────────────────────
# Pipeline: exactitud contra soluciones analíticas
# ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "function, y0, bounds, solution",
    [
        ("y", 1.0, (0, 1), lambda x: np.exp(x)),  # crecimiento
        ("-2*y", 1.0, (0, 2), lambda x: np.exp(-2 * x)),  # decaimiento
        ("2*x", 0.0, (0, 3), lambda x: x**2),  # cuadratura pura
        ("cos(x)", 0.0, (0, np.pi), lambda x: np.sin(x)),  # trigonométrica
        ("-x*y", 1.0, (0, 2), lambda x: np.exp(-(x**2) / 2)),  # coeficiente variable
    ],
)
@pytest.mark.parametrize("method", SCALAR_METHODS)
def test_exactitud_soluciones_analiticas(method, function, y0, bounds, solution):
    outcome = make_outcome(method, function=function, y0=y0, bounds=bounds, n=200)
    assert outcome["status"] == "success"
    x = np.array(outcome["result"]["x"])
    y = np.array(outcome["result"]["y"])
    tol = {1: 5e-2, 2: 1e-3, 3: 1e-4, 4: 1e-6}[ORDER[method]]
    assert np.max(np.abs(y - solution(x))) < tol


@pytest.mark.parametrize("method", SCALAR_METHODS)
def test_orden_observado(method):
    """El error global en x = b debe escalar como h^p (ejemplo de la teoría)."""
    ns = (160, 320) if method.startswith("adams") else (40, 80)
    if method == "adams_moulton_2":
        ns = (320, 640)  # su régimen asintótico empieza tarde
    errs = [abs(make_outcome(method, n=n)["result"]["value"] - exact(2.0)) for n in ns]
    assert np.log2(errs[0] / errs[1]) == pytest.approx(ORDER[method], abs=0.2)


def test_sistema_oscilador_armonico_conserva_energia():
    """y0'' = -y0 → (y0, y1) = (cos x, -sin x); y0² + y1² = 1."""
    outcome = make_outcome(
        "rk4_system",
        function=["y1", "-y0"],
        y0=[1.0, 0.0],
        bounds=(0, 2 * np.pi),
        n=400,
    )
    Y = np.array(outcome["result"]["y"])
    x = np.array(outcome["result"]["x"])
    assert np.allclose(Y[:, 0], np.cos(x), atol=1e-7)
    assert np.allclose(Y[:, 1], -np.sin(x), atol=1e-7)
    assert np.max(np.abs(Y[:, 0] ** 2 + Y[:, 1] ** 2 - 1)) < 1e-7


def test_sistema_de_una_ecuacion_igual_a_escalar():
    """rk4_system con m = 1 debe coincidir con rk4."""
    a = make_outcome("rk4")["result"]["y"]
    b = make_outcome("rk4_system", function=["y - x**2 + 1"], y0=[0.5])["result"]["y"]
    assert np.allclose(a, b)


# ────────────────────────────────────────────────────────────────
# Pipeline: equivalencias de entrada
# ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", SCALAR_METHODS)
def test_n_y_h_equivalentes(method):
    by_n = make_outcome(method, n=20)["result"]
    by_h = make_outcome(method, h=0.1)["result"]
    assert np.array_equal(by_n["x"], by_h["x"]) and np.allclose(by_n["y"], by_h["y"])


def test_bounds_plano_o_anidado():
    flat = make_outcome("rk4", bounds=(0, 2))["result"]["y"]
    nested = NumericalMethod("ode_ivp", {**make_input("rk4"), "bounds": [[0, 2]]})
    nested.validate_input()
    assert np.allclose(flat, nested.execute()["result"]["y"])


# ────────────────────────────────────────────────────────────────
# Paridad con el módulo legacy (eliminar junto con core/ode/)
# ────────────────────────────────────────────────────────────────
# Solo Euler, Heun y RK4: el legacy implementa mal RK2 (orden 1), AB3 (arranque
# de orden 1) y AM2 (es Heun). h = 0.125 evita el desfase del while del legacy.


@pytest.mark.parametrize("method", ["euler", "heun", "rk4"])
def test_paridad_con_legacy(method):
    legacy = NumericalMethod(
        "ode",
        {
            "mode": "function",
            "function": "y - x**2 + 1",
            "x0": 0.0,
            "y0": 0.5,
            "x_end": 2.0,
            "h": 0.125,
            "calculation_mode": method,
        },
    )
    legacy.validate_input()
    old = legacy.execute()["result"]
    new = make_outcome(method, h=0.125)["result"]
    assert np.allclose(old["x"], new["x"]) and np.allclose(
        old["y"], new["y"], atol=1e-12
    )


# ────────────────────────────────────────────────────────────────
# Errores: cada capa responde con su tipo
# ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "overrides, error_type",
    [
        ({"function": "y + z"}, "ConstructionError"),
        ({"calculation_mode": "shooting"}, "ConstructionError"),
        ({"bounds": [2, 0]}, "ConstructionError"),
        ({"step": {"mode": "fixed", "h": 0.3}}, "ConstructionError"),
        (
            {"calculation_mode": "euler", "function": ["y1", "-y0"], "y0": [1, 0]},
            "ValidationError",
        ),
        (
            {
                "calculation_mode": "adams_bashforth_3",
                "step": {"mode": "fixed", "n": 2},
            },
            "ValidationError",
        ),
        (
            {
                "step": {
                    "mode": "adaptive",
                    "tol": 1e-6,
                    "h0": 0.1,
                    "hmin": 1e-4,
                    "hmax": 0.5,
                }
            },
            "ValidationError",
        ),
    ],
)
def test_errores_por_capa(overrides, error_type):
    data = {**make_input("rk4"), **overrides}
    with pytest.raises(Exception) as exc:
        nm = NumericalMethod("ode_ivp", data)
        nm.validate_input()
    assert type(exc.value).__name__ == error_type


def test_divergencia_es_error_controlado():
    """Inestabilidad numérica → ExecutionError normalizado, no excepción."""
    nm = NumericalMethod(
        "ode_ivp", make_input("euler", function="-100*y", y0=1.0, bounds=(0, 40), n=400)
    )
    nm.validate_input()
    with np.errstate(over="ignore", invalid="ignore"):
        outcome = nm.execute()
    assert outcome["status"] == "error" and outcome["error_type"] == "ExecutionError"
    assert "diverged" in outcome["message"]


# ────────────────────────────────────────────────────────────────
# Renderer + UIContract
# ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", SCALAR_METHODS)
def test_renderer_bloques_escalar(method):
    payload = Renderer().render(method, make_outcome(method)["result"])
    types = [b["type"] for b in payload["blocks"]]
    assert types == ["scalar", "plot", "markdown"]
    assert payload["blocks"][2]["format"] == "latex"


def test_renderer_bloques_sistema():
    result = make_outcome("rk4_system", function=["y1", "-y0"], y0=[1, 0], n=20)[
        "result"
    ]
    types = [b["type"] for b in Renderer().render("rk4_system", result)["blocks"]]
    assert types == ["plot", "markdown", "vector"]


@pytest.mark.parametrize("method", METHODS)
def test_contract_devuelve_div(method):
    if method == "rk4_system":
        outcome = make_outcome(method, function=["y1", "-y0"], y0=[1, 0])
    else:
        outcome = make_outcome(method)
    view = contract.resolve(method, outcome)
    assert isinstance(view, html.Div) and view.className == "result-container"
    assert len(view.children) == 3


@pytest.mark.parametrize("method", METHODS)
def test_contract_error_devuelve_div(method):
    outcome = {
        "status": "error",
        "error_type": "ValidationError",
        "message": "Test error",
        "context": {},
    }
    view = contract.resolve(method, outcome)
    assert isinstance(view, html.Div)


# ────────────────────────────────────────────────────────────────
# UI → backend: el formulario produce el mismo resultado que la API
# ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", SCALAR_METHODS)
def test_formulario_ui_igual_a_api(method):
    form = build_input_data(method, ["y - x**2 + 1"], [0.5], 0, 2, "n", 10, None)
    assert form == make_input(method)
    view = run_ivp(form)
    assert isinstance(view, html.Div) and view.className == "result-container"


def test_formulario_ui_sistema_con_h():
    form = build_input_data("rk4_system", ["y1", "-y0"], [1, 0], 0, 1, "h", None, 0.1)
    outcome = NumericalMethod("ode_ivp", form)
    outcome.validate_input()
    assert outcome.execute()["result"]["N"] == 10
