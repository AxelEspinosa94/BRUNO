import pytest

from core.exceptions import CatalogLoadError, ValidationError
from strategies.validators.ode_ivp.ode_ivp_validators import IVPValidator


def make(**overrides):
    data = {
        "mode": "function",
        "calculation_mode": "euler",
        "function": "y - x**2 + 1",
        "bounds": [0, 2],
        "y0": 0.5,
        "step": {"mode": "fixed", "n": 10},
    }
    data.update(overrides)
    return data


SYSTEM = {
    "function": ["-4*y0 + 3*y1 + 6", "-2.4*y0 + 1.6*y1 + 3.6"],
    "y0": [0, 0],
}
ADAPTIVE = {"mode": "adaptive", "tol": 1e-6, "h0": 0.1, "hmin": 1e-6, "hmax": 0.5}


@pytest.fixture(scope="module")
def validator():
    return IVPValidator()


# ---------------------------------------------------------------- valid inputs
@pytest.mark.parametrize(
    "method",
    [
        "euler",
        "heun",
        "rk2",
        "rk4",
        "adams_bashforth_2",
        "adams_bashforth_3",
        "adams_moulton_2",
    ],
)
def test_scalar_methods_accept_standard_problem(validator, method):
    assert validator.validate(make(calculation_mode=method)) is True


def test_rk4_system_accepts_system(validator):
    assert validator.validate(make(calculation_mode="rk4_system", **SYSTEM)) is True


def test_min_steps_uses_h_as_well_as_n(validator):
    data = make(
        calculation_mode="adams_bashforth_3", step={"mode": "fixed", "h": 2 / 3}
    )
    assert validator.validate(data) is True  # N = 3


# ---------------------------------------------------------------- method checks
@pytest.mark.parametrize(
    "overrides, match",
    [
        ({"mode": "table"}, "mode='function'"),
        ({"calculation_mode": "shooting"}, "Not supported method"),
        ({"calculation_mode": None}, "Not supported method"),
    ],
)
def test_rejects_mode_and_method(validator, overrides, match):
    with pytest.raises(ValidationError, match=match):
        validator.validate(make(**overrides))


@pytest.mark.parametrize("method", ["euler", "rk4", "adams_bashforth_2"])
def test_adaptive_step_not_supported_yet(validator, method):
    with pytest.raises(ValidationError, match="step mode 'adaptive'"):
        validator.validate(make(calculation_mode=method, step=ADAPTIVE))


@pytest.mark.parametrize("method", ["euler", "heun", "rk4", "adams_moulton_2"])
def test_systems_only_where_released(validator, method):
    with pytest.raises(ValidationError, match="systems of equations"):
        validator.validate(make(calculation_mode=method, **SYSTEM))


@pytest.mark.parametrize(
    "method, n, ok",
    [
        ("adams_bashforth_2", 1, False),
        ("adams_bashforth_2", 2, True),
        ("adams_bashforth_3", 2, False),
        ("adams_bashforth_3", 3, True),
        ("adams_moulton_2", 2, False),
        ("adams_moulton_2", 3, True),
    ],
)
def test_multistep_min_steps(validator, method, n, ok):
    data = make(calculation_mode=method, step={"mode": "fixed", "n": n})
    if ok:
        assert validator.validate(data) is True
    else:
        with pytest.raises(ValidationError, match="at least"):
            validator.validate(data)


# ---------------------------------------------------------------- catalog
def test_catalog_entry_for_every_supported_method(validator):
    for method in validator.catalog["SUPPORTED_METHODS"]:
        assert isinstance(validator.catalog.get(method), dict), method


def test_every_rule_path_resolves(validator):
    from app.utils.rule_loader import load_rule

    paths = list(validator.catalog.get("common_rules", []))
    for method in validator.catalog["SUPPORTED_METHODS"]:
        paths += validator.catalog[method].get("rules", [])
    for path in paths:
        assert callable(load_rule(path)), path


def test_missing_catalog_entry_is_a_catalog_error(validator, monkeypatch):
    monkeypatch.setitem(
        validator.catalog,
        "SUPPORTED_METHODS",
        validator.catalog["SUPPORTED_METHODS"] + ["ralston"],
    )
    with pytest.raises(CatalogLoadError, match="ralston"):
        validator.validate(make(calculation_mode="ralston"))
