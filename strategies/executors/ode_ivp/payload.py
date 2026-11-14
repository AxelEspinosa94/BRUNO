import pandas as pd


def build_payload(instance, result, expression=None):
    """
    Unifies the output of every IVP method (same role as integration/payload.py).

    result: {"W": ndarray (N + 1, m), "extras": {column: list of length N + 1}}
    """
    W = result["W"]
    x = instance.x
    m = instance.m
    components = ["y"] if m == 1 else [f"y{j}" for j in range(m)]

    data = {"i": list(range(len(x))), "x": x}
    for j, name in enumerate(components):
        data[name] = W[:, j]
    data.update(result.get("extras", {}))
    table = pd.DataFrame(data)

    payload = {
        "calculation_mode": instance.calculation_mode,
        "expression": expression or "",
        "expression_format": "latex",
        "table": table,
        "x": x.tolist(),
        "components": components,
        "a": instance.a,
        "b": instance.b,
        "N": instance.N,
        "h": instance.h,
        "plot_type": "curve" if m == 1 else "multi_curve",
    }

    if m == 1:
        payload["y"] = W[:, 0].tolist()
        payload["value"] = float(W[-1, 0])  # y(b), rendered as scalar
    else:
        payload["y"] = W.tolist()  # one row per mesh point
        payload["solution"] = W[-1].tolist()  # y(b), rendered as vector

    return payload
