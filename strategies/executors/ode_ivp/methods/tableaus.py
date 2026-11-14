"""
Butcher tableaus for explicit Runge-Kutta methods.

Coefficients are written as fractions so they read exactly like the theory
documents (docs/theory/ode_ivp/theory_rk2.md, theory_rk4.md) and are converted
to float arrays once, at import time.
"""

from dataclasses import dataclass
from fractions import Fraction as Fr
from typing import Dict

import numpy as np


@dataclass(frozen=True)
class ButcherTableau:
    name: str
    order: int
    A: np.ndarray  # (s, s), strictly lower triangular for explicit methods
    b: np.ndarray  # (s,)
    c: np.ndarray  # (s,)

    @property
    def stages(self) -> int:
        return len(self.b)


def _tableau(name, order, A, b):
    A = np.array([[float(v) for v in row] for row in A], dtype=float)
    b = np.array([float(v) for v in b], dtype=float)
    c = A.sum(axis=1)  # row-sum condition c_i = sum_j a_ij
    return ButcherTableau(name=name, order=order, A=A, b=b, c=c)


TABLEAUS: Dict[str, ButcherTableau] = {
    "euler": _tableau("euler", 1, A=[[0]], b=[1]),
    "heun": _tableau(
        "heun",
        2,
        A=[[0, 0], [1, 0]],
        b=[Fr(1, 2), Fr(1, 2)],
    ),
    "midpoint": _tableau(
        "midpoint",
        2,
        A=[[0, 0], [Fr(1, 2), 0]],
        b=[0, 1],
    ),
    "ralston": _tableau(
        "ralston",
        2,
        A=[[0, 0], [Fr(2, 3), 0]],
        b=[Fr(1, 4), Fr(3, 4)],
    ),
    "rk4": _tableau(
        "rk4",
        4,
        A=[[0, 0, 0, 0], [Fr(1, 2), 0, 0, 0], [0, Fr(1, 2), 0, 0], [0, 0, 1, 0]],
        b=[Fr(1, 6), Fr(1, 3), Fr(1, 3), Fr(1, 6)],
    ),
}


def get_tableau(name: str) -> ButcherTableau:
    try:
        return TABLEAUS[name]
    except KeyError as exc:
        raise KeyError(
            f"Unknown Butcher tableau '{name}'. Available: {sorted(TABLEAUS)}"
        ) from exc
