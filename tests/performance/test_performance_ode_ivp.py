import statistics
import time

import pytest

from app.utils.build_function import build_ode_function
from core.base_method import NumericalMethod

METHODS = [
    "euler",
    "heun",
    "rk2",
    "rk4",
    "adams_bashforth_2",
    "adams_bashforth_3",
    "adams_moulton_2",
]

# Pasos para medir performance (un método O(N) debe resolverlos en < 1 s)
N_PERF = 10_000

# Límites de tiempo por método (segundos). Medidos ~0.15–0.6 s en una VM;
# el margen absorbe ruido del sistema y de CI.
LIMITS = {
    "euler": 1.5,
    "heun": 2.0,
    "rk2": 2.0,
    "rk4": 3.0,
    "adams_bashforth_2": 1.5,
    "adams_bashforth_3": 1.5,
    "adams_moulton_2": 2.5,
    "rk4_system": 3.0,
}


def make_outcome(method, n, function="sin(x) - 0.1*y", y0=1.0, bounds=(0, 10)):
    nm = NumericalMethod(
        method="ode_ivp",
        input_data={
            "mode": "function",
            "calculation_mode": method,
            "function": function,
            "bounds": list(bounds),
            "y0": y0,
            "step": {"mode": "fixed", "n": n},
        },
    )
    nm.validate_input()
    t0 = time.perf_counter()
    outcome = nm.execute()
    elapsed = time.perf_counter() - t0
    return outcome, elapsed


def median_time(method, n, runs=3, **kwargs):
    return statistics.median(make_outcome(method, n, **kwargs)[1] for _ in range(runs))


@pytest.fixture(scope="session", autouse=True)
def warmup_sympy():
    build_ode_function("sin(x) + exp(y)")


# ────────────────────────────────────────────────────────────────
# PERFORMANCE: tiempo de ejecución
# ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", METHODS)
def test_performance(method):
    """Cada método resuelve N_PERF pasos por debajo de su límite."""
    outcome, elapsed = make_outcome(method, N_PERF)
    assert outcome["status"] == "success"
    assert (
        elapsed < LIMITS[method]
    ), f"{method} tardó {elapsed:.4f}s (límite {LIMITS[method]}s)"


def test_performance_sistema():
    """Sistema de 6 ecuaciones (3 osciladores desacoplados) con RK4."""
    funcs = ["y1", "-y0", "y3", "-y2", "y5", "-y4"]
    outcome, elapsed = make_outcome(
        "rk4_system", 5_000, function=funcs, y0=[1, 0, 1, 0, 1, 0]
    )
    assert outcome["status"] == "success"
    assert elapsed < LIMITS["rk4_system"], f"rk4_system tardó {elapsed:.4f}s"


# ────────────────────────────────────────────────────────────────
# PERFORMANCE: complejidad lineal en N
# ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", METHODS)
def test_performance_escala_lineal(method):
    """
    Cuadruplicar N debe costar ~4x (lineal), lejos de ~16x (cuadrático).
    Detecta, por ejemplo, historiales que se copian completos en cada paso.
    """
    t_small = median_time(method, 1_000)
    t_large = median_time(method, 4_000)
    ratio = t_large / t_small
    assert ratio < 8.0, f"{method}: N x4 costó {ratio:.1f}x"


# ────────────────────────────────────────────────────────────────
# PERFORMANCE: determinismo temporal
# ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", METHODS)
def test_performance_determinismo(method):
    """Dos mediciones consecutivas no deben diferir por órdenes de magnitud."""
    t1 = median_time(method, 2_000)
    t2 = median_time(method, 2_000)
    ratio = max(t1, t2) / min(t1, t2)
    assert ratio < 3.0, f"Variación temporal excesiva: {ratio:.2f}x"


# ────────────────────────────────────────────────────────────────
# PERFORMANCE: la función no afecta el orden de magnitud
# ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", ["euler", "rk4", "adams_moulton_2"])
def test_performance_funcion_costosa(method):
    """Una f(x, y) con muchas funciones elementales sigue dentro de 2x el límite."""
    fn = "sin(x)*cos(y) + exp(-x**2)*sqrt(1 + y**2) + log(1 + x**2) - tanh(y)"
    outcome, elapsed = make_outcome(method, N_PERF, function=fn)
    assert outcome["status"] == "success"
    assert elapsed < 2 * LIMITS[method], f"{method} tardó {elapsed:.4f}s con f costosa"
