"""
Generic explicit Runge-Kutta solver driven by a Butcher tableau.

Euler, Heun, midpoint (rk2), RK4 and RK4 for systems are the same algorithm with
different coefficients, so they share this code. The pure functions `rk_step`
and `rk_solve` are also used as starters by the Adams methods and, later, by the
shooting method in ode_bvp.
"""

import numpy as np

from .common import check_finite, require_fixed_mesh
from .tableaus import ButcherTableau, get_tableau


def rk_step(F, x, w, h, tableau: ButcherTableau):
    """One step from (x, w). Returns (w_next, K) with K of shape (stages, m)."""
    A, b, c = tableau.A, tableau.b, tableau.c
    K = np.empty((tableau.stages, w.size), dtype=float)
    for i in range(tableau.stages):
        K[i] = F(x + c[i] * h, w + h * (A[i, :i] @ K[:i]))
    return w + h * (b @ K), K


def rk_solve(F, x, y0, tableau: ButcherTableau, method_name: str = "RK"):
    """
    Solves on the given uniform mesh x (length N + 1).
    Returns W of shape (N + 1, m) and K of shape (N, stages, m).
    """
    N = len(x) - 1
    h = x[1] - x[0] if N > 0 else 0.0
    W = np.empty((N + 1, y0.size), dtype=float)
    K_all = np.empty((N, tableau.stages, y0.size), dtype=float)
    W[0] = y0
    for i in range(N):
        W[i + 1], K_all[i] = rk_step(F, x[i], W[i], h, tableau)
        check_finite(W[i + 1], x[i + 1], method_name)
    return W, K_all


def run_explicit_rk(instance, params):
    """Executor entry point. params: {"tableau": "<name>"}."""
    x, _ = require_fixed_mesh(instance)
    tableau = get_tableau(params["tableau"])
    W, K = rk_solve(instance.F, x, instance.y0, tableau, instance.calculation_mode)

    extras = {}
    if instance.m == 1:  # stage slopes only make sense as columns for a scalar IVP
        for s in range(tableau.stages):
            extras[f"k{s + 1}"] = [None] + K[:, s, 0].tolist()

    return {"W": W, "extras": extras}
