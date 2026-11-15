import os

from app.utils.catalog_loader import load_catalog
from app.utils.rule_loader import load_rule
from core.exceptions import CatalogLoadError, ExecutionError

from .payload import build_payload


class IVPExecutor:
    """
    Dispatches each IVP method to its executor through the catalog
    (same pattern as IntegrationExecutor). Each executor receives
    (instance, params) and returns {"W": ..., "extras": ...}.
    """

    def __init__(self):
        base = os.path.dirname(__file__)
        path = os.path.join(base, "ode_ivp_executors_catalog.json")
        self.catalog = load_catalog(path)

    def run(self, instance):
        calculation_mode = instance.calculation_mode
        supported = self.catalog.get("SUPPORTED_METHODS", [])
        if calculation_mode not in supported:
            raise ExecutionError(
                f"Not supported method for IVP: {calculation_mode}. "
                f"Must be one of: {', '.join(supported)}."
            )

        entry = self.catalog.get(calculation_mode)
        if not isinstance(entry, dict) or "executor" not in entry:
            raise CatalogLoadError(
                f"Method '{calculation_mode}' is listed as supported "
                "but has no executor in the IVP executors catalog."
            )

        execute = load_rule(entry["executor"])
        result = execute(instance, entry.get("params", {}))
        return build_payload(instance, result, entry.get("expression"))
