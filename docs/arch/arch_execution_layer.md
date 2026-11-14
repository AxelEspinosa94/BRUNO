
---

# **Execution Layer Documentation**

# **Table of Contents**

-   [Purpose](#1-purpose)
-   [Interface](#2-interface)
-   [Internal Structure](#3-internal-structure)
    -   [General Pattern Dispatcher by Calculation Mode](#31-dispatcher-by-calculation_mode)
    -   [Helper Functions](#32-helper-functions)
-   [Return Structure](#4-return-structure)
-   [Error Handling](#5-error-handling)
-   [Existing Executors](#6-existing-executors)
-   [Executor Philosophy](#7-executor-philosophy)
-   [Adding a New Method](#8-adding-a-new-method-new-method)
-   [Adding a New Module](#9-adding-a-new-module-new-module)
-   [What the Executor Does not do](#10-what-the-executor-does-not-do)
-   [Summary](#11-summary)

## **1. Purpose**

The Execution Layer is where **numerical computation actually happens**.  
It receives a fully constructed and validated method instance and executes the algorithm associated with its `calculation_mode`.

It is the **only** layer that performs mathematical work.  
All other layers (constructor, validator, renderer, UI contract) exist to support it.

---

## **2. Interface**

All executors implement `ExecutorProtocol`:

```python
class ExecutorProtocol(Protocol):
    def run(self, instance: Any) -> Any: ...
```

`NumericalMethod` invokes it as:

```python
result = self.executor.run(self.method_instance)
```

Where `self.method_instance` is the constructor instance containing all parsed and validated attributes (`df`, `xk`, `interval`, `n`, `calculation_mode`, etc.).

---

## **3. Internal Structure**

### **3.1 Dispatcher by `calculation_mode`**

Each executor uses an internal dispatcher to route execution to the correct private method:

```python
class <Module>Executor:
    def run(self, instance) -> Dict[str, Any]:
        dispatch = {
            "method_a": self._run_method_a,
            "method_b": self._run_method_b,
        }
        fn = dispatch.get(instance.calculation_mode)
        if fn is None:
            raise ExecutionError(f"Unknown calculation_mode: '{instance.calculation_mode}'")
        return fn(instance)
```

Each `_run_<mode>()` is **self‑contained**:

- reads from the instance  
- performs the numerical computation  
- returns a structured result dict  

Methods are not overloaded because each numerical technique has distinct algorithms and helper functions.

---

### **3.2 Helper Functions**

Executors may define private helper functions for reusable sub‑computations:

```python
def _lagrange_multiplier(self, df, i, xk) -> tuple: ...
def _eval_lagrange(self, df, x) -> float: ...
def _newton_expression(self, coef, x) -> str: ...
def _eval_newton(self, coef, x_nodes, xk) -> float: ...
def _spline_expression(self, M, x, y, k, hk) -> str: ...
def _eval_spline_curve(self, M, x, y, h) -> tuple: ...
def _hermite_expression(self, Q, z, m) -> str: ...
def _eval_hermite(self, Q, z, m, xk) -> float: ...
```

Helpers are **private** and belong exclusively to the executor that uses them.  
They are not shared across modules.

---

## **4. Return Structure**

Each `_run_<mode>()` returns a dict.  
The keys vary by method, but the general structure for interpolation is:

```python
{
    "value":      float,
    "expression": str,
    "table":      pd.DataFrame,
    "x":          list[float],
    "y":          list[float],
    "x_nodes":    list[float],
    "y_nodes":    list[float],
}
```

Other domains (integration, ODE, linear algebra, etc.) define their own key sets.

The contract is:

> **Whatever keys the executor returns, `UIContract._build_blocks()` will render them.**

Executors do not need to coordinate with the UI layer.

---

## **5. Error Handling**

Executors raise `ExecutionError` when a computation fails:

```python
raise ExecutionError("Descriptive message about what failed")
```

`NumericalMethod.execute()` wraps the executor call:

```python
def execute(self):
    try:
        result = self.executor.run(self.method_instance)
        return {"status": "success", "result": result}
    except Exception as e:
        return ErrorNormalizer.normalize(
            exception=e,
            method_name=self.method,
            input_data=self.input
        )
```

Executors **never** return error dicts — they raise exceptions.  
Normalization is handled by `NumericalMethod`.

---

## **6. Existing Executors**

| Module            | Class                       | File Path                                                   |
|-------------------|-----------------------------|-------------------------------------------------------------|
| Interpolation     | `InterpolationExecutor`     | `strategies/executors/interpolation_executors.py`          |
| Integration       | `IntegrationExecutor`       | `strategies/executors/integration/integration_executors.py` |
| ODE (legacy)      | `ODEExecutor`               | `strategies/executors/ode_executors.py`                    |
| ODE IVP           | `IVPExecutor`               | `strategies/executors/ode_ivp/ode_ivp_executors.py`        |
| Linear Algebra    | `LinearAlgebraExecutor`     | `strategies/executors/linear_algebra_executors.py`         |
| Non‑Linear        | `NonLinearExecutor`         | `strategies/executors/non_linear_executors.py`             |
| Numerical Deriv.  | `NumericalDerivativeExecutor` | `strategies/executors/numerical_derivative_executors.py` |

### **6.1 Catalog-based executors (Integration, ODE IVP)**

Catalog-based modules replace the internal `dispatch` dict with an executor catalog and
move each algorithm to its own file:

```
strategies/executors/<module>/
    <module>_executors.py           # <Module>Executor: catalog lookup + payload
    <module>_executors_catalog.json # method → executor path (+ params)
    payload.py                      # unifies the output of every method
    methods/
        <family>.py                 # algorithms, plain functions
```

```
run() → catalog lookup → load_rule(executor) → execute(instance[, params]) → build_payload()
```

ODE IVP extends this with per-method `params`, so several methods share one implementation:

```json
"heun": {
  "executor": "strategies.executors.ode_ivp.methods.explicit_rk.run_explicit_rk",
  "params": { "tableau": "heun" },
  "expression": "$w_{i+1} = w_i + \\frac{h}{2}(k_1 + k_2)$"
}
```

- Euler, Heun, RK2 (midpoint), RK4 and RK4 for systems are one generic explicit RK solver
  driven by Butcher tableaus (`methods/tableaus.py`). Adding Ralston is a catalog entry.
- Executors work on vectors: a scalar IVP is a system with `m = 1`.
- `rk_step` / `rk_solve` are pure functions reused as starters by the Adams methods
  (and, later, by the shooting method in `ode_bvp`).
- Method functions return `{"W": ndarray (N + 1, m), "extras": {column: values}}`;
  `payload.py` turns that into `x`, `y`, `value` (scalar) or `solution` (system),
  `expression` and `table`.
- Executors raise `ExecutionError` when the solution stops being finite (divergence).

For catalog-based modules, principle 2 below ("one `_run_<mode>()` per method") becomes
"one executor function per algorithm family, selected through the catalog".

---

# **7. Executor Philosophy**

The Execution Layer follows three core principles:

**1. One executor per module**
Each numerical domain (interpolation, integration, ODE, etc.) has its own executor class.

This keeps algorithms isolated and prevents cross‑module coupling.

**2. One `_run_<mode>()` per method**
Each numerical method is implemented as a private function:

- `_run_simpson_1_3`
- `_run_gauss_legendre`
- `_run_clenshaw_curtis`
- `_run_<new-method>`

This ensures clarity, maintainability, and testability.

**3. Dispatcher + helpers**
The executor is structured as:

```
run() → dispatcher → _run_<mode>() → helpers → result dict
```

This pattern is consistent across all modules.

---

# **8. Adding a New Method (`<new-method>`)**

When adding a new calculation mode to an existing module:

**1. Add a private `_run_<new-method>()`**
Example:

```python
def _run_<new-method>(self, instance):
    ...
    return {"value": result, ...}
```

**2. Add helper functions if needed**
Example:

```python
def _eval_<new-method>(self, ...): ...
```

**3. Register the method in the dispatcher**

```python
dispatch = {
    "simpson_1_3": self._run_simpson_1_3,
    "clenshaw_curtis": self._run_clenshaw_curtis,
    "<new-method>": self._run_<new-method>,
}
```

**4. Ensure the Validator supports the new method**
Add validation rules in:

```
strategies/validators/<module>/<module>_validation_catalog.json
```

**5. Add the method to `SUPPORTED_MODES`**
Always append at the end:

```json
"supported_modes": [
  "simpson_1_3",
  "clenshaw_curtis",
  "<new-method>"
]
```

**6. Register the executor in `method_catalog.json`**

```json
{
  "<new-method>": {
    "classExecutor": "strategies.executors.<module>_executors.<Module>Executor",
    ...
  }
}
```

---

# **9. Adding a New Module (`<new-module>`)**

When adding a new numerical domain:

**1. Create the folder structure**

```
strategies/
  executors/
    <new-module>/
      __init__.py
      <new-module>_executors.py
```

**2. Create the executor class**

```python
class <NewModule>Executor:
    def run(self, instance):
        dispatch = {
            "<new-method>": self._run_<new-method>,
        }
        ...
```

**3. Add validation rules**
Create:

```
strategies/validators/<new-module>/
    <new-module>_validators.py
    <new-module>_validation_catalog.json
```

**4. Register the module in `method_catalog.json`**

```json
{
  "<new-method>": {
    "classConstructor": "strategies.constructors.<new-module>.<NewModule>Constructor",
    "classInputValidator": "strategies.validators.<new-module>.<NewModule>Validator",
    "classExecutor": "strategies.executors.<new-module>.<NewModule>Executor"
  }
}
```

**5. Add the new method to `SUPPORTED_MODES`**

**6. Write unit tests for the new executor**

---

# **10. What the Executor Does NOT Do**

- Does **not** validate input  
- Does **not** parse or transform `input_data`  
- Does **not** render UI output  
- Does **not** return error dicts  
- Does **not** modify the constructor instance  

Executors **only** compute and return results.

---

# **11. Summary**

The Execution Layer is:

- modular  
- scalable  
- catalog‑driven  
- easy to extend with `<new-module>` and `<new-method>`  
- consistent across all numerical domains  
- cleanly separated from validation and UI layers  

This architecture ensures NA‑Engine remains maintainable, predictable, and ready for future numerical methods.

---
