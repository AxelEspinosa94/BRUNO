from typing import List, Union

import numpy as np
import sympy as sp

from core.exceptions import ConstructionError


def build_function(func_str: str, symbols: List[str] = None):
    """
    Build a numpy-vectorized callable from a sympy expression.
    If symbols is None, assume single variable 'x'.
    If symbols provided, create lambdify with those symbols in order.

    Raises ConstructionError if the expression cannot be parsed or uses
    symbols outside `symbols`. Unused symbols are allowed here; modules that
    need every symbol to appear (e.g. integration) check that separately.
    """
    if not isinstance(func_str, str) or not func_str.strip():
        raise ConstructionError("Function must be a non-empty string.")

    if symbols is None:
        symbols = ["x"]

    syms = sp.symbols(list(symbols))  # always a list, even for one symbol

    try:
        f_sym = sp.sympify(func_str)
    except (sp.SympifyError, SyntaxError, TypeError) as exc:
        hint = (
            " For conditional functions use Piecewise((a, cond), (b, True))."
            if " if " in f" {func_str} "
            else ""
        )
        raise ConstructionError(
            f"Could not parse function '{func_str}'.{hint}"
        ) from exc

    if not isinstance(f_sym, sp.Expr):
        raise ConstructionError(f"'{func_str}' is not a valid mathematical expression.")

    unknown = sorted(str(s) for s in f_sym.free_symbols - set(syms))
    if unknown:
        hint = " For Euler's number use E or exp(...)." if "e" in unknown else ""
        raise ConstructionError(
            f"Function uses unknown symbols {unknown}; allowed: {list(symbols)}.{hint}"
        )

    # If constant
    if f_sym.is_Number:
        const = float(f_sym)

        def const_fn(*args):
            # args are arrays; return array of same shape as first arg
            a0 = np.asarray(args[0])
            return np.full_like(a0, const, dtype=float)

        return const_fn

    return sp.lambdify(syms, f_sym, "numpy")


def build_grid(
    f_str: str,
    symbols: List[str] = None,
    bounds: Union[List, List[List]] = None,
    n: Union[int, List[int]] = 1,
):
    f = build_function(f_str, symbols)
    dim = len(symbols) if symbols else 1

    # Normalizar n: mismo número de puntos por eje o uno distinto por eje
    n_list = [n] * dim if isinstance(n, int) else n

    # Construir ejes 1D
    axes = [np.linspace(a, b, ni + 1) for (a, b), ni in zip(bounds, n_list)]

    if dim == 1:
        x = axes[0]
        y = f(x)
        return x, y

    # Malla N-dimensional con indexing='ij' para que los ejes coincidan
    # con el orden de tus variables y con axis=i en la integración posterior
    grids = np.meshgrid(*axes, indexing="ij")
    y = f(*grids)  # broadcasting elemento a elemento gracias a lambdify

    return grids, y


def ode_symbols(m: int) -> List[str]:
    """Symbol names for an m-dimensional IVP: ['x', 'y'] or ['x', 'y0', ..., 'y{m-1}']."""
    if m < 1:
        raise ConstructionError("The system must have at least one equation.")
    return ["x", "y"] if m == 1 else ["x"] + [f"y{i}" for i in range(m)]


def build_ode_function(funcs: Union[str, List[str]]):
    """F(x, y) -> np.ndarray of shape (m,), built on top of build_function."""
    exprs = [funcs] if isinstance(funcs, str) else list(funcs)
    if not exprs:
        raise ConstructionError("At least one equation is required.")

    names = ode_symbols(len(exprs))
    fs = []
    for i, e in enumerate(exprs):
        try:
            fs.append(build_function(e, names))
        except ConstructionError as exc:
            raise ConstructionError(f"Equation {i}: {exc}") from exc

    def F(x, y):
        y = np.atleast_1d(np.asarray(y, dtype=float))
        return np.array([f(x, *y) for f in fs], dtype=float)

    return F, names
