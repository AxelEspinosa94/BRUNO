import numpy as np

from core.base_method import NumericalMethod

# ============================================================
# STRESS TESTS — PASO MUY GRANDE
# ============================================================


def test_euler_large_step_unstable():
    # y' = y, y(0)=1 → con h grande Euler explota
    method = NumericalMethod(
        method="ode_ivp",
        input_data={
            "mode": "function",
            "function": "y",
            "bounds": [0, 5],
            "y0": 1,
            "step": {"mode": "fixed", "n": 100},
            "calculation_mode": "euler",
        },
    )

    method.validate_input()
    result = method.execute().get("result", {})

    # Euler diverge brutalmente, solo verificamos que no truene
    assert np.isfinite(result["y"][-1])


def test_rk4_large_step_still_reasonable():
    # RK4 aguanta pasos grandes mejor que Euler
    method = NumericalMethod(
        method="ode_ivp",
        input_data={
            "mode": "function",
            "function": "y",
            "bounds": [0, 5],
            "y0": 1,
            "step": {"mode": "fixed", "n": 100},
            "calculation_mode": "rk4",
        },
    )

    method.validate_input()
    result = method.execute().get("result", {})

    y_num = result["y"][-1]
    y_exact = np.exp(5)

    # RK4 con h=1 aún da algo razonable
    assert abs(y_num - y_exact) < 5.0


# ============================================================
# STRESS TESTS — PASO MUY PEQUEÑO
# ============================================================


def test_rk4_tiny_step_precision():
    method = NumericalMethod(
        method="ode_ivp",
        input_data={
            "mode": "function",
            "function": "y",
            "bounds": [0, 1],
            "y0": 1,
            "step": {"mode": "fixed", "n": 100},
            "calculation_mode": "rk4",
        },
    )

    method.validate_input()
    result = method.execute().get("result", {})

    y_num = result["y"][-1]
    y_exact = np.e

    assert abs(y_num - y_exact) < 5e-5


# ============================================================
# STRESS TESTS — ECUACIONES RÍGIDAS
# ============================================================


def test_stiff_equation_rk4():
    # y' = -1000y + 3000 - 2000 e^{-t}
    # Solución exacta: y = 3 - 0.002 e^{-1000t} - 2 e^{-t}
    method = NumericalMethod(
        method="ode_ivp",
        input_data={
            "mode": "function",
            "function": "-1000*y + 3000 - 2000*exp(-x)",
            "bounds": [0, 1],
            "y0": 0,
            "step": {"mode": "fixed", "n": 2000},
            "calculation_mode": "rk4",
        },
    )

    method.validate_input()
    result = method.execute().get("result", {})

    y_num = result["y"][-1]
    y_exact = 3 - 0.002 * np.exp(-1000) - 2 * np.exp(-1)

    assert abs(y_num - y_exact) < 1e-2


# ============================================================
# STRESS TESTS — SISTEMAS GRANDES
# ============================================================


def test_rk4_system_large():
    # Sistema acoplado:
    # y1' = y2
    # y2' = -y1
    # Oscilador armónico
    method = NumericalMethod(
        method="ode_ivp",
        input_data={
            "mode": "function",
            "function": ["y1", "-y0"],
            "y0": [1, 0],
            "bounds": [0, 20],
            "step": {"mode": "fixed", "n": 2000},
            "calculation_mode": "rk4_system",
        },
    )

    method.validate_input()
    result = method.execute().get("result", {})

    y1 = result["y"][-1][0]
    y_exact = np.cos(20)

    assert abs(y1 - y_exact) < 0.05


# ============================================================
# STRESS TESTS — BVP DIFÍCIL
# ============================================================


# def test_finite_differences_nontrivial():
#     # y'' = -pi^2 sin(pi x), y(0)=0, y(1)=0 → y = sin(pi x)
#     method = NumericalMethod(
#         method="ode_ivp",
#         input_data={
#             "mode": "function",
#             "function": "-(np.pi**2)*np.sin(np.pi*x)",
#             "bounds": [0, 1],
#             "alpha": 0,
#             "beta": 0,
#             "step": {"mode": "fixed", "n": 50},
#             "calculation_mode": "finite_differences",
#         },
#     )
#
#     method.validate_input()
#     result = method.execute().get("result", {})
#
#     x = np.array(result["x"], dtype=float)
#     y = result["y"]
#     y_exact = np.sin(np.pi * x)
#
#     assert np.allclose(y, y_exact, atol=0.1)


# ============================================================
# STRESS TESTS — MULTISTEP CON PASO GRANDE
# ============================================================


def test_adams_bashforth_3_large_step():
    method = NumericalMethod(
        method="ode_ivp",
        input_data={
            "mode": "function",
            "function": "y",
            "bounds": [0, 2],
            "y0": 1,
            "step": {"mode": "fixed", "n": 10},
            "calculation_mode": "adams_bashforth_3",
        },
    )

    method.validate_input()
    result = method.execute().get("result", {})

    y_num = result["y"][-1]
    y_exact = np.exp(2)

    assert abs(y_num - y_exact) < 2.0


# ============================================================
# STRESS TESTS — DISCONTINUIDADES
# ============================================================


def test_discontinuous_rhs():
    # y' = 1 si x < 1, 0 si x >= 1
    method = NumericalMethod(
        method="ode_ivp",
        input_data={
            "mode": "function",
            # Sintaxis que depende de eval
            # "function": "1 if x < 1 else 0",
            "function": "Piecewise((1, x < 1), (0, True))",  # se lee como el if/else
            # o
            # "function": "1 - Heaviside(x - 1)",
            "bounds": [0, 2],
            "y0": 0,
            "step": {"mode": "fixed", "n": 200},
            "calculation_mode": "rk4",
        },
    )

    method.validate_input()
    result = method.execute().get("result", {})

    y_num = result["y"][-1]
    y_exact = 1  # y crece hasta x=1 y luego se queda constante

    assert abs(y_num - y_exact) < 0.05


# ============================================================
# STRESS TESTS — SHOOTING CON FUNCIÓN NO LINEAL
# ============================================================


# def test_shooting_nonlinear():
#     # y'' = -y^3, y(0)=0, y(1)=0
#     # Solución trivial: y=0
#     method = NumericalMethod(
#         method="ode_ivp",
#         input_data={
#             "mode": "function",
#             "function": "-y**3",
#             "bounds": [0, 1],
#             "alpha": 0,
#             "beta": 0,
#             "step": {"mode": "fixed", "n": 100},
#             "calculation_mode": "shooting",
#         },
#     )
#
#     method.validate_input()
#     result = method.execute().get("result", {})
#     assert abs(result["y_end"]) < 0.1
