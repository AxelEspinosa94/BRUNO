import os

from app.utils.catalog_loader import load_catalog
from app.utils.rule_loader import load_rule
from core.exceptions import CatalogLoadError, ValidationError


class IVPValidator:
    """
    Method-level validation for initial value problems.

    Runs after the IVP constructor, so `input_data` is already structurally valid
    (function parses, bounds, y0 and step are well formed). Each rule receives
    (input_data, params), where params come from the method's catalog entry.

    Catalog layout:
        SUPPORTED_METHODS : methods released in this layer
        common_rules      : rules applied to every method
        <method>.params   : method data used by the rules (step_modes, systems, ...)
        <method>.rules    : extra rules for that method only
    """

    def __init__(self):
        base = os.path.dirname(__file__)
        path = os.path.join(base, "ode_ivp_validation_catalog.json")
        self.catalog = load_catalog(path)

    def validate(self, input_data):
        if input_data.get("mode") != "function":
            raise ValidationError("ODE IVP only supports mode='function'.")

        calculation_mode = input_data.get("calculation_mode")
        supported = self.catalog.get("SUPPORTED_METHODS", [])
        if calculation_mode not in supported:
            raise ValidationError(
                f"Not supported method for IVP: {calculation_mode}. "
                f"Must be one of: {', '.join(supported)}."
            )

        entry = self.catalog.get(calculation_mode)
        if not isinstance(entry, dict):
            raise CatalogLoadError(
                f"Method '{calculation_mode}' is listed as supported "
                "but has no entry in the IVP validation catalog."
            )

        params = entry.get("params", {})
        rules = self._as_list(self.catalog.get("common_rules")) + self._as_list(
            entry.get("rules")
        )

        for rule_path in rules:
            load_rule(rule_path)(input_data, params)

        return True

    @staticmethod
    def _as_list(rules):
        if rules is None:
            return []
        return [rules] if isinstance(rules, str) else list(rules)
