"""
Adams-Bashforth (explicit) and Adams-Moulton (PECE) multistep methods.

Starting values come from an explicit RK method (RK4 by default), as described in
docs/theory/ode_ivp/theory_adams-bashford.md and theory_adams-moulton-2.md.
"""

import numpy as np

from core.exceptions import ExecutionError

from .common import check_finite, require_fixed_mesh
from .explicit_rk import rk_solve
from .tableaus import get_tableau

# w_{i+1} = w_i + h * sum_j beta_j f_{i-j}
AB_COEFFS = {
    1: np.array([1.0]),
    2: np.array([3.0, -1.0]) / 2.0,
    3: np.array([23.0, -16.0, 5.0]) / 12.0,
    4: np.array([55.0, -59.0, 37.0, -9.0]) / 24.0,
}

# 2-step Adams-Moulton (order 3): w_{i+1} = w_i + h/12 (5 f_{i+1} + 8 f_i - f_{i-1})
AM2_COEFFS = np.array([5.0, 8.0, -1.0]) / 12.0

# Milne's estimate for the AB3-AM2 pair: y - w_c ≈ -(w_c - w_p) / 10
MILNE_FACTOR = {3: -1.0 / 10.0}


def _start(instance, k, starter):
    """Mesh, step, solution array with w_0..w_{k-1} filled, and slope history."""
    x, h = require_fixed_mesh(instance)
    N = len(x) - 1
    if N < k:
        raise ExecutionError(
            f"{instance.calculation_mode} needs at least {k} steps; got N = {N}."
        )

    W = np.empty((N + 1, instance.m), dtype=float)
    Fv = np.empty_like(W)
    W[:k], _ = rk_solve(instance.F, x[:k], instance.y0, get_tableau(starter), starter)
    for j in range(k):
        Fv[j] = instance.F(x[j], W[j])
    return x, h, N, W, Fv


def _ab_increment(Fv, i, beta):
    """sum_j beta_j f_{i-j}"""
    return sum(beta[j] * Fv[i - j] for j in range(len(beta)))


def run_adams_bashforth(instance, params):
    """params: {"steps": k, "starter": "rk4"}"""
    k = params["steps"]
    beta = AB_COEFFS[k]
    starter = params.get("starter", "rk4")
    x, h, N, W, Fv = _start(instance, k, starter)

    for i in range(k - 1, N):
        W[i + 1] = W[i] + h * _ab_increment(Fv, i, beta)
        check_finite(W[i + 1], x[i + 1], instance.calculation_mode)
        Fv[i + 1] = instance.F(x[i + 1], W[i + 1])

    source = [starter if j < k else f"ab{k}" for j in range(N + 1)]
    source[0] = "y0"
    return {"W": W, "extras": {"source": source}}


def run_adams_moulton_pece(instance, params):
    """
    AM2 corrector with an AB predictor, one correction per step (PECE).
    params: {"predictor_steps": 3, "starter": "rk4"}
    """
    p = params.get("predictor_steps", 3)
    beta_p = AB_COEFFS[p]
    starter = params.get("starter", "rk4")
    k = max(p, 2)  # AM2 itself needs f_i and f_{i-1}
    x, h, N, W, Fv = _start(instance, k, starter)

    predicted = np.full_like(W, np.nan)
    for i in range(k - 1, N):
        wp = W[i] + h * _ab_increment(Fv, i, beta_p)  # P
        fp = instance.F(x[i + 1], wp)  # E
        W[i + 1] = W[i] + h * (
            AM2_COEFFS[0] * fp + AM2_COEFFS[1] * Fv[i] + AM2_COEFFS[2] * Fv[i - 1]  # C
        )
        check_finite(W[i + 1], x[i + 1], instance.calculation_mode)
        Fv[i + 1] = instance.F(x[i + 1], W[i + 1])  # E
        predicted[i + 1] = wp

    extras = {}
    source = [starter if j < k else "am2" for j in range(N + 1)]
    source[0] = "y0"
    extras["source"] = source
    if instance.m == 1:
        extras["predicted"] = [
            None if np.isnan(v) else float(v) for v in predicted[:, 0]
        ]
        if p in MILNE_FACTOR:
            est = MILNE_FACTOR[p] * (W[:, 0] - predicted[:, 0])
            extras["milne_estimate"] = [None if np.isnan(v) else float(v) for v in est]
    return {"W": W, "extras": extras}
