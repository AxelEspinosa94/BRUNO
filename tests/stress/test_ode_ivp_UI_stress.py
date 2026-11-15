# tests/stress/test_ode_ivp_UI_stress.py
"""
Stress tests for ode_ivp through the same path as the UI
(build_input_data → NumericalMethod → UIContract).

They do not check accuracy; they check that extreme or malformed inputs end
in a result or in a controlled error view, never in an uncaught exception.
"""

import numpy as np
import pytest
from dash import html

from app.callbacks.ode_ivp.ode_ivp_callbacks import MAX_EQ, build_input_data, run_ivp
from core.base_method import NumericalMethod
from core.contract import UIContract

contract = UIContract()

IVP_MODES = [
    "euler",
    "heun",
    "rk2",
    "rk4",
    "adams_bashforth_2",
    "adams_bashforth_3",
    "adams_moulton_2",
]


def make_outcome(method, function, y0, bounds, **step):
    """Como el callback: los errores de construcción también se vuelven outcome."""
    data = {
        "mode": "function",
        "calculation_mode": method,
        "function": function,
        "bounds": list(bounds),
        "y0": y0,
        "step": {"mode": "fixed", **step},
    }
    try:
        nm = NumericalMethod("ode_ivp", data)
        nm.validate_input()
    except Exception as e:
        return {
            "status": "error",
            "error_type": type(e).__name__,
            "message": str(e),
            "context": data,
        }
    with np.errstate(all="ignore"):
        return nm.execute()


def assert_view(method, outcome):
    view = contract.resolve(method, outcome)
    assert isinstance(view, html.Div)
    return view


# ═══════════════════════════════════════════════════════════════
# 1. INTERVALOS GRANDES
# ═══════════════════════════════════════════════════════════════


@pytest.mark.parametrize("mode", IVP_MODES)
def test_intervalo_grande(mode):
    """[0, 500] con h = 0.5 en una EDO de crecimiento lento."""
    outcome = make_outcome(mode, "sin(x) - 0.01*y", 1.0, (0, 500), h=0.5)
    assert outcome["status"] == "success"
    assert_view(mode, outcome)


@pytest.mark.parametrize("mode", IVP_MODES)
def test_intervalo_negativo_y_desplazado(mode):
    outcome = make_outcome(mode, "x + y", 1.0, (-50, -40), n=100)
    assert outcome["status"] == "success"
    assert outcome["result"]["x"][0] == -50.0 and outcome["result"]["x"][-1] == -40.0


# ═══════════════════════════════════════════════════════════════
# 2. PASOS MUY PEQUEÑOS / MUY POCOS
# ═══════════════════════════════════════════════════════════════


@pytest.mark.parametrize("mode", IVP_MODES)
def test_h_muy_pequeno(mode):
    outcome = make_outcome(mode, "sin(x) + y", 0.0, (0, 1e-3), h=1e-7)
    assert outcome["status"] == "success"
    assert outcome["result"]["N"] == 10_000


@pytest.mark.parametrize("mode", IVP_MODES)
def test_n_minimo_por_metodo(mode):
    """N = k pasos: lo mínimo que acepta cada método (k = 1, 2 o 3)."""
    k = {"adams_bashforth_2": 2, "adams_bashforth_3": 3, "adams_moulton_2": 3}.get(
        mode, 1
    )
    outcome = make_outcome(mode, "y", 1.0, (0, 1), n=k)
    assert outcome["status"] == "success" and len(outcome["result"]["y"]) == k + 1
    rejected = make_outcome(mode, "y", 1.0, (0, 1), n=k - 1) if k > 1 else None
    if rejected:
        assert (
            rejected["status"] == "error"
            and rejected["error_type"] == "ValidationError"
        )


# ═══════════════════════════════════════════════════════════════
# 3. VOLUMEN DE PASOS
# ═══════════════════════════════════════════════════════════════


@pytest.mark.parametrize("mode", IVP_MODES)
def test_volumen_pasos(mode):
    outcome = make_outcome(mode, "-y + x", 1.0, (0, 10), n=50_000)
    assert outcome["status"] == "success"
    assert len(outcome["result"]["table"]) == 50_001
    assert_view(mode, outcome)


# ═══════════════════════════════════════════════════════════════
# 4. INESTABILIDAD Y DESBORDAMIENTO → ERROR CONTROLADO
# ═══════════════════════════════════════════════════════════════


@pytest.mark.parametrize("mode", IVP_MODES)
def test_problema_rigido_fuera_de_region_de_estabilidad(mode):
    """y' = -1000 y con h = 0.1: z = -100, fuera de la región de todos los métodos."""
    outcome = make_outcome(mode, "-1000*y", 1.0, (0, 50), h=0.1)
    assert outcome["status"] == "error" and outcome["error_type"] == "ExecutionError"
    assert "diverged" in outcome["message"]
    assert_view(mode, outcome)


@pytest.mark.parametrize("mode", IVP_MODES)
def test_explosion_en_tiempo_finito(mode):
    """y' = y², y(0) = 1 explota en x = 1: el intervalo [0, 2] la cruza."""
    outcome = make_outcome(mode, "y**2", 1.0, (0, 2), n=2_000)
    assert outcome["status"] in ("success", "error")
    if outcome["status"] == "error":
        assert outcome["error_type"] == "ExecutionError"
    assert_view(mode, outcome)


@pytest.mark.parametrize("mode", IVP_MODES)
def test_dominio_invalido_de_f(mode):
    """sqrt de un valor negativo produce NaN: debe reportarse como divergencia."""
    outcome = make_outcome(mode, "sqrt(y)", -1.0, (0, 1), n=10)
    assert outcome["status"] == "error" and outcome["error_type"] == "ExecutionError"
    assert_view(mode, outcome)


# ═══════════════════════════════════════════════════════════════
# 5. FUNCIONES DIFÍCILES
# ═══════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    "function",
    [
        "Piecewise((1, x < 1), (0, True))",  # discontinua
        "Abs(x - 1)",  # no derivable
        "sin(100*x)",  # altamente oscilatoria
        "exp(-x**2)*cos(y)",  # sin primitiva elemental
        "1/(1 + y**2)",  # no lineal acotada
        "0",  # constante cero
    ],
)
@pytest.mark.parametrize("mode", IVP_MODES)
def test_funciones_dificiles(mode, function):
    outcome = make_outcome(mode, function, 0.0, (0, 2), n=400)
    assert outcome["status"] == "success"
    assert np.all(np.isfinite(outcome["result"]["y"]))


# ═══════════════════════════════════════════════════════════════
# 6. SISTEMAS
# ═══════════════════════════════════════════════════════════════


def test_sistema_maximo_de_ecuaciones():
    """MAX_EQ ecuaciones acopladas en cadena (lo máximo que permite la UI)."""
    funcs = [f"-y{i} + y{(i + 1) % MAX_EQ}" for i in range(MAX_EQ)]
    outcome = make_outcome("rk4_system", funcs, [1.0] * MAX_EQ, (0, 5), n=1_000)
    assert outcome["status"] == "success"
    assert np.array(outcome["result"]["y"]).shape == (1_001, MAX_EQ)
    assert_view("rk4_system", outcome)


def test_sistema_lotka_volterra_largo():
    """Depredador-presa: no lineal, periódico, interval largo."""
    funcs = ["1.1*y0 - 0.4*y0*y1", "0.1*y0*y1 - 0.4*y1"]
    outcome = make_outcome("rk4_system", funcs, [10.0, 10.0], (0, 200), n=20_000)
    assert outcome["status"] == "success"
    Y = np.array(outcome["result"]["y"])
    assert np.all(Y > 0)  # las poblaciones nunca se vuelven negativas


def test_sistema_con_ecuacion_constante():
    outcome = make_outcome("rk4_system", ["1", "y0", "0"], [0, 0, 5], (0, 1), n=10)
    assert outcome["status"] == "success"
    assert np.allclose(outcome["result"]["solution"], [1.0, 0.5, 5.0])


# ═══════════════════════════════════════════════════════════════
# 7. ENTRADAS MALFORMADAS DESDE LA UI → SIEMPRE UNA VISTA
# ═══════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    "exprs, y0s, a, b, step_by, n, h",
    [
        ([None], [None], None, None, "n", None, None),  # formulario vacío
        ([""], [0.5], 0, 2, "n", 10, None),  # función vacía
        (["y +* 2"], [0.5], 0, 2, "n", 10, None),  # sintaxis inválida
        (["import os"], [0.5], 0, 2, "n", 10, None),  # código, no matemáticas
        (["__import__('os')"], [0.5], 0, 2, "n", 10, None),  # intento de inyección
        (["1 if x < 1 else 0"], [0.5], 0, 2, "n", 10, None),  # ternario de Python
        (["y"], [0.5], 2, 2, "n", 10, None),  # intervalo degenerado
        (["y"], [0.5], 0, 2, "n", -5, None),  # n negativo
        (["y"], [0.5], 0, 2, "n", 2.5, None),  # n no entero
        (["y"], [0.5], 0, 2, "h", None, 0),  # h = 0
        (["y"], [0.5], 0, 2, "h", None, 10),  # h > intervalo
        (["y"], [0.5], 0, 2, "h", None, None),  # h vacío
        (["y"], [float("inf")], 0, 2, "n", 10, None),  # y0 infinito
        (["y0", "y1"], [1, 0], 0, 2, "n", 10, None),  # sistema con método escalar
    ],
)
def test_entradas_malformadas_devuelven_vista(exprs, y0s, a, b, step_by, n, h):
    view = run_ivp(build_input_data("rk4", exprs, y0s, a, b, step_by, n, h))
    assert isinstance(view, html.Div)
    assert view.className == "result-error"


# ═══════════════════════════════════════════════════════════════
# 8. DETERMINISMO
# ═══════════════════════════════════════════════════════════════


@pytest.mark.parametrize("mode", IVP_MODES)
def test_determinismo(mode):
    r1 = make_outcome(mode, "sin(x)*y", 1.0, (0, 5), n=500)["result"]
    r2 = make_outcome(mode, "sin(x)*y", 1.0, (0, 5), n=500)["result"]
    assert np.array_equal(r1["y"], r2["y"])


def test_ejecuciones_repetidas_no_acumulan_estado():
    """La misma instancia de UIContract y del registro se reutilizan entre clics."""
    values = [
        make_outcome("rk4", "y - x**2 + 1", 0.5, (0, 2), n=10)["result"]["value"]
        for _ in range(50)
    ]
    assert len(set(values)) == 1
