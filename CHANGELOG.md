# Changelog

---

## BRUNO v0.4.0 (Unreleased) — Stage 5 Minor 4: ODE IVP Module Refactor into Catalog-based approach
**Release date:** TBD

### Added
- New `ode_ivp` module (constructor, validator and executor) following the catalog-based architecture, built in parallel to the legacy `ode` module so the latter keeps working during the migration.
  - Methods: Euler, Heun, RK2 (midpoint), RK4, RK4 for systems, Adams–Bashforth 2 and 3, Adams–Moulton 2 (PECE with AB3 predictor).
  - Scalar IVPs (`function` as a string, symbols `x`, `y`) and systems (`function` as a list, symbols `x`, `y0`, …, `y{m-1}`).
  - Canonical step definition: `step = {"mode": "fixed", "n" | "h"}`; `"adaptive"` mode is accepted by the constructor and reserved for future methods (e.g. RKF45).
- Generic explicit Runge–Kutta executor driven by Butcher tableaus (`methods/tableaus.py`): Euler, Heun, midpoint and RK4 share one implementation.
- Adams methods use RK4 starting values; Adams–Moulton 2 reports the predicted value and Milne's local error estimate.
- Divergence detection in IVP executors (`ExecutionError` with the `x` where the solution overflowed).
- `build_ode_function` and `ode_symbols` in `app/utils/build_function.py`.
- Theory documents for IVP fundamentals and every IVP method (`docs/theory/ode_ivp/`).
- Unit tests: `test_build_function.py`, `test_ivp_constructor.py`, `test_ivp_validator.py`, `test_ivp_executor.py` (theory examples and observed order of convergence).
- Integration, performance and stress tests for `ode_ivp`:
  - `tests/integration/test_ode_ivp_STAGE-5.py`: catalog consistency across layers, accuracy against analytic solutions, observed order, error type per layer, renderer/contract output, UI form ↔ API equivalence and parity with the legacy module (Euler, Heun, RK4).
  - `tests/performance/test_performance_ode_ivp.py`: time limits for N = 10,000 steps, linear scaling in N, timing stability and costly right-hand sides.
  - `tests/stress/test_ode_ivp_UI_stress.py`: large intervals, tiny h, 50,000 steps, stiff/unstable problems, finite-time blow-up, NaN-producing functions, discontinuous and oscillatory functions, systems up to the UI limit and malformed UI inputs.
- Catalog-based UI for `ode_ivp` (`app/layout/ode_ivp/`, `app/callbacks/ode_ivp/`, `app/tooltips/ode_ivp.py`) following the Integration layout pattern:
  - Dynamic equation rows (`y' = f(x, y)`, `y(a)`); extra equations are enabled only for methods listed in `system_methods`.
  - Step defined either by number of steps `n` or step size `h`.
  - Mounted as a temporary "EDO (IVP)" tab next to the legacy ODE tab, so both can be compared until `ode_bvp` is released.
- `build_function` suggests `Piecewise((a, cond), (b, True))` when a Python ternary (`a if cond else b`) is used.

### Improved
- `build_function` now rejects unknown symbols (via sympy `free_symbols`) and unparseable expressions with a `ConstructionError` at construction time instead of failing during evaluation; it hints at `E`/`exp(...)` when `e` is used.
- Integration and ODE IVP constructors check `calculation_mode` against the `SUPPORTED_METHODS` list of their catalog before building anything. Each layer keeps its own list, so a method can be held back in any layer.
- ODE IVP validator rules receive per-method `params` from the catalog (`step_modes`, `systems`, `min_steps`) plus module-wide `common_rules`.
- IVP mesh is built as `x_i = a + i·h` from a normalized `(N, h)`, removing the floating-point drift of the legacy `while x < x_end` loop.
- Expressions can declare `"expression_format": "latex"`: the renderer forwards the format and `UIContract` renders them with MathJax (display mode); expressions without it keep the code-block style. ODE IVP and Integration now send real LaTeX (Integration previously sent pseudo-LaTeX text such as `∫_0^1 f(x) dx ≈ …`).

### Fixed (relative to the legacy `ode` module)
- Legacy `adams_moulton_2` actually implemented Heun's method (order 2); `ode_ivp` implements the 2-step Adams–Moulton formula (order 3).
- Legacy Adams–Bashforth 3 used a first-order starter, losing its third order; `ode_ivp` uses RK4 starting values.
- Legacy `rk2` used `k2 = f(x + h, y + h·k1)` with `y + h·k2`, a first-order scheme; `ode_ivp` implements the midpoint method (order 2).

### Notes
The legacy ODE tab is still available. Once `ode_ivp` reaches parity and `ode_bvp` is released, the temporary "EDO (IVP)" tab will replace it, and BVP methods (shooting, finite differences) will move to `ode_bvp`, reusing the `ode_ivp` RK solver.
Upcoming work will target:

- Renderer / UIContract generalization for ODE IVP results (shared with other modules):
  - Value block label is hardcoded as `f(x_k)`; IVPs should show `y(b)` (use `a`, `b` from the payload).
  - Systems: `y` is a list of rows and is plotted as a single trace; plot one trace per entry of `components` when `plot_type == "multi_curve"`.
  - Plot labels are hardcoded as `P(x)` (interpolation).
  - The `table` DataFrame is not rendered by any module (only `x_nodes`/`y_nodes` tables are); IVP extra columns contain `None` for rows where they do not apply.
- `ode_bvp` module (shooting, finite differences)
- Adaptive methods (RKF45) on top of the `"adaptive"` step mode

---

## BRUNO v0.3.0 — Stage 5 Minor 3: Integral N-Dimensional Extension + Monte Carlo addition in Integration Module
**Release date:** 2026-10-10

### Added
- Extended Trapezoid, Simpson, Romberg, Gauss-Legendre and Clenshaw-Curtis to the calculation of multivariable functions
  - $x_{i}$ use of variables added when integrating
  - Extended interval attribute to bounds (set of intervals) for calculation
  - $\mathcal{R}^{2}$ and $\mathcal{R}^{3}$ plots renderization
- New Monte Carlo integration method added to the Integration module that follows current architecture.

### Improved
- Extended Calculation to multiple dimensions
- Renderer compatibility improved for $\mathcal{R}^{2}$ and $\mathcal{R}^{3}$ functions.

### Notes
This minor release completes the third phase of Stage 5, focusing on UI/UX modernization through catalog-driven architecture.
Upcoming work will target:

- Adaptative Simpson & Gauss Quadratures
- Advanced Romberg Integration module
- Interpolation refactor to the catalog-based layouts
- Numerical Derivative refactor to the catalog-based layouts
- Full migration of ODE/BVP modules to catalog-based layouts
- Export tools (PDF, LaTeX, CSV)
- Full i18n migration (English/Spanish)
- Stage 5: Advanced numerical methods (RKF45, Adams–Bashforth, Clenshaw–Curtis extensions)
- Stage 5: Non linear refactor + New Methods
- Stage 6: Linear Algebra refactor to the catalog-base aproach


---

## BRUNO v0.2.0 — Stage 5 Minor 2: Integral, Documentation and About Refactoring into Catalog-based approach + Clenshaw-Curtis Quadrature addition in Integration Module
**Release date:** 2026-09-02

### Added
- Catalog‑based architecture for:
  - Integration module UI
  - Documentation module (theory viewer)
  - About module (project metadata, badges, changelog, credits)
- Full refactor of Integration module layouts into declarative catalog format.
- New Clenshaw–Curtis Quadrature method added to the Integration module.
- New catalog-driven content loader for About section (header, badges, version, changelog, credits, notes).
- New catalog-driven documentation renderer with multi-section support.
- Unified catalog schema for future modules (ODE, BVP, Nonlinear, Linear Algebra).

### Improved
- Integration module UI now uses dynamic catalog-driven layouts, reducing hardcoded components.
- Documentation module now supports structured multi-block theory pages with improved readability.
- About module fully migrated to catalog-based rendering with component mapping and parameter injection.
- Renderer compatibility improved for catalog-based modules (plot, table, expression, markdown).
- Better separation of content and layout across all refactored modules.
- Enhanced maintainability: modules now follow a consistent declarative pattern.

### Notes
This minor release completes the second phase of Stage 5, focusing on UI/UX modernization through catalog-driven architecture.
Upcoming work will target:

- Full migration of ODE/BVP modules to catalog-based layouts
- Multi-curve plotting improvements for ODE systems
- Export tools (PDF, LaTeX, CSV)
- Full i18n migration (English/Spanish)
- Stage 6: Advanced numerical methods (RKF45, Adams–Bashforth, Clenshaw–Curtis extensions)

---

## BRUNO v0.1.5 — Stage 5 Patch 2: Unified Output Renderer + Matrix LaTeX + Header Upgrade
**Release date:** 2026-08-04

### Added
- Unified output rendering system with multi-block support.
- Tooltip integration for all output blocks (scalar, vector, matrix, plots, tables, markdown).
- LaTeX rendering for:
  - Matrices (L, U, P)
  - Vectors (solution vectors)
  - Matrix groups (LU, QR)
- New `matrix_expression` block type for mathematical outputs. (Will substitute the original `matrix` block type)
- Dynamic version display in the header using Git tags.
- New header design with gradient typography and theme-aware styling.

### Improved
- Output blocks now follow consistent styling per module (Interpolation, Integration, Linear Algebra, etc.).
- Enhanced KEY_DISPATCH to support multi-output detection (value, expression, table, plot, matrix groups).
- Renderer now produces structured multi-block payloads for the UIContract.
- Matrix group renderer upgraded to include raw values + LaTeX expressions.
- Cleaned and standardized result-area layout for readability and coherence.
- Improved dark/light mode compatibility across output components.

### Notes
This patch completes the foundational work for Stage 5’s rendering architecture.  
Upcoming patches will focus on:
- Advanced ODE/BVP visualization
- Export tools (PDF, LaTeX, CSV)
- Responsive layout improvements
- Full i18n migration (English/Spanish)

---

## BRUNO v0.1.4 — Stage 5 Patch 1: UI/UX Foundation + Tooltip System
**Release date:** 2026-07-28

### Added
- Global UI/UX foundation for all modules.
- Unified styling for:
  - Numeric inputs
  - Text inputs
  - Dropdowns
  - RadioItems
  - Buttons
  - Error states
  - Hover/disabled states
- Complete tooltip system with:
  - Tooltip icon component
  - Tooltip input wrapper
  - Light/dark mode support
  - Hover/click interactions
- Full tooltip coverage for:
  - Nonlinear Equations
  - Numerical Derivatives
  - Numerical Integration
  - Interpolation
  - Linear Algebra

### Improved
- Standardized layout containers across modules.
- Cleaned legacy CSS and removed inconsistent styles.
- Improved readability and consistency of UI components.
- Prepared structure for future i18n (English migration not yet active).

### Notes
This patch establishes the visual and structural foundation for Stage 5.  
Upcoming patches will focus on new numerical methods, advanced ODE/BVP features, responsiveness, export tools, and i18n.

---

## BRUNO v0.1.3 — About Section

### Added
- About section

### Improved
- Consistent styling across modules
- Better separation of UI concerns
- Cleaner callback structure

### Notes
This version is a patch before improvements to take place in Stage 5.

---

## BRUNO v0.1.0 — MVP

### Added
- New dynamic documentation module
- Markdown rendering with MathJax support
- Dark/light theme support for theoretical content
- New folder structure under docs/theory/<module>/<method>.md
- UI integration with module/method selectors
- Architecture documentation for the documentation module

### Improved
- Consistent styling across modules
- Better separation of UI concerns
- Cleaner callback structure

### Notes
This version completes Stage 4 of BRUNO.

---