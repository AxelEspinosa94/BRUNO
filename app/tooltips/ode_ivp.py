# ============================================================
# Tooltips — ODE Initial Value Problems (ode_ivp) Module
# ============================================================

TOOLTIPS = {
    # ───────────────────────────────────────────────
    # Selector de método
    # ───────────────────────────────────────────────
    "ivp-method": (
        "Selecciona el método para resolver el problema de valor inicial.\n\n"
        "• Euler: orden 1, el más simple.\n"
        "• Heun / RK2 (punto medio): orden 2.\n"
        "• RK4: orden 4, el más usado en la práctica.\n"
        "• RK4 (Sistema): permite varias ecuaciones acopladas.\n"
        "• Adams-Bashforth 2 / 3: multipaso explícitos, arrancan con RK4.\n"
        "• Adams-Moulton 2: predictor-corrector (AB3 + AM2), orden 3."
    ),
    # ───────────────────────────────────────────────
    # Modo de entrada
    # ───────────────────────────────────────────────
    "ivp-input-mode": (
        "El módulo trabaja con ecuaciones diferenciales escritas como funciones f(x, y). "
        "Un sistema se captura agregando ecuaciones con '+ equation'."
    ),
    # ───────────────────────────────────────────────
    # Ecuaciones y condiciones iniciales
    # ───────────────────────────────────────────────
    "ivp-fn": (
        "Escribe la derivada y' = f(x, y) y su valor inicial y(a).\n\n"
        "• Una ecuación: usa x y y. Ej: y - x**2 + 1\n"
        "• Sistema: usa x, y0, y1, ... Ej: y1  y  -y0\n"
        "• Funciones: sin, cos, exp, log, sqrt; constante e = E\n"
        "• Condicionales: Piecewise((1, x < 1), (0, True))\n\n"
        "Solo RK4 (Sistema) admite más de una ecuación."
    ),
    # ───────────────────────────────────────────────
    # Intervalo
    # ───────────────────────────────────────────────
    "ivp-interval": (
        "Intervalo [a, b] donde se busca la solución. "
        "La condición inicial se da en x = a y debe cumplirse a < b."
    ),
    # ───────────────────────────────────────────────
    # Paso
    # ───────────────────────────────────────────────
    "ivp-step-by": (
        "Define la malla con el número de pasos n o con el tamaño de paso h. "
        "Ambos están ligados por h = (b - a) / n."
    ),
    "ivp-n": (
        "Número de pasos. Adams-Bashforth 2 requiere n ≥ 2; "
        "Adams-Bashforth 3 y Adams-Moulton 2 requieren n ≥ 3."
    ),
    "ivp-h": (
        "Tamaño de paso. Debe dividir exactamente al intervalo; "
        "si no, se sugiere el h válido más cercano. "
        "Un h demasiado grande puede volver inestable al método."
    ),
    # ───────────────────────────────────────────────
    # Botón de ejecución
    # ───────────────────────────────────────────────
    "ivp-run-btn": "Resuelve el problema con el método y los parámetros seleccionados.",
}
