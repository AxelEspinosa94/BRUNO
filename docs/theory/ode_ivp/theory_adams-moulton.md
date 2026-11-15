---

# Adams–Bashforth Methods — Theory

---

# **Table of Contents**

- [Overview](#overview)
- [Mathematical Foundations](#mathematical-foundations)
  - [Linear Multistep Methods](#linear-multistep-methods)
  - [Derivation by Polynomial Interpolation](#derivation-by-polynomial-interpolation)
  - [Starting Values](#starting-values)
- [Adams–Bashforth 2-Step (AB2)](#adamsbashforth-2-step-ab2)
- [Adams–Bashforth 3-Step (AB3)](#adamsbashforth-3-step-ab3)
- [Algorithm](#algorithm)
- [Error Analysis](#error-analysis)
  - [Local Truncation Error](#local-truncation-error)
  - [Zero-Stability and Convergence](#zero-stability-and-convergence)
- [Stability](#stability)
- [Properties](#properties)
- [Example](#example)
- [Why Adams–Bashforth Methods Matter](#why-adamsbashforth-methods-matter)
- [Sources](#sources)

---

# Overview

Let $f(x,y)$ be continuous and Lipschitz in $y$ on $[a,b] \times \mathbb{R}$. The goal is to approximate the solution of

$$
y' = f(x,y), \qquad a \le x \le b, \qquad y(a) = \alpha
$$

at equally spaced mesh points $x_i = a + ih$.

Runge–Kutta methods discard everything except $w_i$ after each step and evaluate $f$ several times inside the new step. **Multistep methods** do the opposite: they **reuse** slopes already computed at previous mesh points. The **Adams–Bashforth** (AB) methods are the explicit members of the Adams family. They need only **one new evaluation of $f$ per step**, regardless of their order.

This document covers the two- and three-step versions implemented in BRUNO (`adams_bashforth_2`, `adams_bashforth_3`). See `theory_ivp.md` for the general theory.

---

# Mathematical Foundations

## Linear Multistep Methods

A linear $k$-step method has the form

$$
w_{i+1} = \sum_{j=0}^{k-1} a_j\, w_{i-j} + h \sum_{j=-1}^{k-1} b_j\, f_{i-j}, \qquad f_j = f(x_j, w_j)
$$

- If $b_{-1} = 0$, $w_{i+1}$ appears only on the left: the method is **explicit**.
- If $b_{-1} \neq 0$, $f_{i+1} = f(x_{i+1}, w_{i+1})$ appears on the right: the method is **implicit**.

Adams methods are the subfamily with $a_0 = 1$ and $a_1 = \dots = a_{k-1} = 0$:

$$
w_{i+1} = w_i + h \sum_j b_j\, f_{i-j}
$$

## Derivation by Polynomial Interpolation

Integrating the ODE over $[x_i, x_{i+1}]$:

$$
y(x_{i+1}) = y(x_i) + \int_{x_i}^{x_{i+1}} f\big(x, y(x)\big)\,dx
$$

The integrand is unknown, but its values $f_i, f_{i-1}, \dots, f_{i-k+1}$ at the previous $k$ mesh points are known. AB methods replace the integrand by the **interpolating polynomial** $P_{k-1}(x)$ through those $k$ points and integrate it over the next interval. Since the interval $[x_i, x_{i+1}]$ lies **outside** the interpolation nodes, this is an **extrapolation**.

With $x = x_i + sh$, the Newton backward-difference form gives

$$
\int_{x_i}^{x_{i+1}} P_{k-1}(x)\,dx = h \sum_{j=0}^{k-1} \gamma_j\, \nabla^j f_i, \qquad \gamma_j = (-1)^j \int_0^1 \binom{-s}{j}\,ds
$$

with $\gamma_0 = 1$, $\gamma_1 = \tfrac12$, $\gamma_2 = \tfrac{5}{12}$, $\gamma_3 = \tfrac38$. Expanding the backward differences $\nabla f_i = f_i - f_{i-1}$, $\nabla^2 f_i = f_i - 2f_{i-1} + f_{i-2}$ produces the coefficients below.

The connection with the interpolation module (Newton's divided differences, `theory_newton.md`) is direct.

## Starting Values

A $k$-step method needs $w_0, w_1, \dots, w_{k-1}$ before it can be applied. Only $w_0 = \alpha$ is given; the rest must come from a **one-step method**.

For the multistep method of order $p$ to keep its order, the starting values must be accurate to $O(h^p)$. A one-step method of order $p - 1$ is formally enough, since it is used for a fixed, finite number of steps. In practice **RK4** is the usual choice because its starting values are much more accurate than needed and it is already available.

---

# Adams–Bashforth 2-Step (AB2)

Linear interpolation through $(x_{i-1}, f_{i-1})$ and $(x_i, f_i)$, integrated over $[x_i, x_{i+1}]$:

$$
\boxed{
w_{i+1} = w_i + \frac{h}{2}\big(3 f_i - f_{i-1}\big), \qquad i = 1, 2, \dots, N-1
}
$$

Starting values: $w_0 = \alpha$, $w_1$ from a one-step method.

Geometric reading: the slope over the next interval is extrapolated linearly from the two most recent slopes, $f_i + \tfrac12(f_i - f_{i-1})$.

---

# Adams–Bashforth 3-Step (AB3)

Quadratic interpolation through $(x_{i-2}, f_{i-2})$, $(x_{i-1}, f_{i-1})$, $(x_i, f_i)$:

$$
\boxed{
w_{i+1} = w_i + \frac{h}{12}\big(23 f_i - 16 f_{i-1} + 5 f_{i-2}\big), \qquad i = 2, 3, \dots, N-1
}
$$

Starting values: $w_0 = \alpha$, $w_1$ and $w_2$ from a one-step method.

In each AB formula the coefficients add up to 1 ($\tfrac{3-1}{2} = 1$, $\tfrac{23-16+5}{12} = 1$). This is the consistency condition: for constant $f$ the method must reproduce $w_{i+1} = w_i + hf$.

---

# Algorithm

**Input:** $f$, endpoints $a, b$, initial value $\alpha$, number of steps $N$, number of steps $k \in \{2, 3\}$.
**Output:** approximations $w_i \approx y(x_i)$ for $i = 0, \dots, N$.

```
h ← (b − a) / N
w₀ ← α
compute w₁, ..., w_{k−1} with RK4
F[j] ← f(x_j, w_j)   for j = 0, ..., k−1        # slope history
for i = k−1, ..., N−1:
    if k = 2: w_{i+1} ← w_i + (h/2)  · (3F[i] − F[i−1])
    if k = 3: w_{i+1} ← w_i + (h/12) · (23F[i] − 16F[i−1] + 5F[i−2])
    F[i+1] ← f(x_{i+1}, w_{i+1})                   # only new evaluation
return (x, w)
```

**Cost:** after the starter, **one** evaluation of $f$ per step. Only the last $k$ slopes need to be stored.

**Validation requirements:** $N \ge k$ (otherwise the method never runs past the starter), and the mesh must be **uniform**: the coefficients are derived for equal spacing.

---

# Error Analysis

## Local Truncation Error

For $y$ sufficiently smooth, there exists $\mu_i \in (x_{i-k+1}, x_{i+1})$ such that

| Method | Local truncation error $\tau_{i+1}(h)$ | Order |
|---|---|---|
| AB2 | $\dfrac{5}{12}\, h^2\, y'''(\mu_i)$ | 2 |
| AB3 | $\dfrac{3}{8}\, h^3\, y^{(4)}(\mu_i)$ | 3 |

In general, the $k$-step Adams–Bashforth method has **order $k$**, with error constant $\gamma_k$ (the next coefficient in the series above).

## Zero-Stability and Convergence

A linear multistep method has the **characteristic polynomial**

$$
\rho(\zeta) = \zeta^k - \sum_{j=0}^{k-1} a_j\, \zeta^{k-1-j}
$$

It is **zero-stable** if all roots of $\rho$ satisfy $|\zeta| \le 1$ and those with $|\zeta| = 1$ are simple (the **root condition**). For Adams methods

$$
\rho(\zeta) = \zeta^k - \zeta^{k-1} = \zeta^{k-1}(\zeta - 1)
$$

with roots $1$ (simple) and $0$: all Adams methods are zero-stable. By the **Dahlquist equivalence theorem**, consistency plus zero-stability gives convergence, so

$$
|y(x_i) - w_i| = O(h^k)
$$

provided the starting values are accurate enough.

---

# Stability

On $y' = \lambda y$, the method becomes a linear recurrence with characteristic equation

$$
\rho(\zeta) - z\,\sigma(\zeta) = 0, \qquad z = h\lambda
$$

where $\sigma(\zeta) = \sum_j b_j \zeta^{k-1-j}$. The method is absolutely stable at $z$ if all roots satisfy $|\zeta| \le 1$.

| Method | $\sigma(\zeta)$ | Real stability interval |
|---|---|---|
| AB2 | $\tfrac12(3\zeta - 1)$ | $-1 \le h\lambda \le 0$ |
| AB3 | $\tfrac{1}{12}(23\zeta^2 - 16\zeta + 5)$ | $-\tfrac{6}{11} \approx -0.545 \le h\lambda \le 0$ |

The stability regions of AB methods are **small** and **shrink as the order grows**, considerably smaller than those of RK methods of the same order (Euler: $-2$, RK4: $-2.785$). This is the price of using just one function evaluation per step.

---

# Properties

- Explicit, linear multistep, **not** self-starting.
- One new function evaluation per step, independent of the order.
- $k$-step AB has order $k$: AB2 order 2, AB3 order 3.
- Zero-stable (root condition holds for all Adams methods).
- Small stability regions that shrink with the order.
- Require a uniform mesh; changing $h$ requires restarting or interpolating the history.
- Require $N \ge k$ and $k - 1$ starting values.

---

# Example

Approximate the solution of

$$
y' = y - x^2 + 1, \qquad 0 \le x \le 2, \qquad y(0) = 0.5
$$

with $h = 0.2$ ($N = 10$), using RK4 for the starting values. Exact solution: $y(x) = (x+1)^2 - 0.5\,e^{x}$.

**Starting values (RK4).** $w_1 = 0.8292933$, $w_2 = 1.2140762$.

**AB2, first step** ($i = 1$):

$$
f_0 = f(0, 0.5) = 1.5, \qquad f_1 = f(0.2,\, 0.8292933) = 1.7892933
$$

$$
w_2 = 0.8292933 + 0.1\,\big(3 \cdot 1.7892933 - 1.5\big) = 1.2160813
$$

(AB2 only takes $w_1$ from RK4; from $w_2$ on, it uses its own values.)

**Results.**

| $x_i$ | AB2 $w_i$ | AB2 $\lvert e \rvert$ | AB3 $w_i$ | AB3 $\lvert e \rvert$ |
|---|---|---|---|---|
| 0.0 | 0.5000000 | 0 | 0.5000000 | 0 |
| 0.2 | 0.8292933* | $5.29 \times 10^{-6}$ | 0.8292933* | $5.29 \times 10^{-6}$ |
| 0.4 | 1.2160813 | $1.99 \times 10^{-3}$ | 1.2140762* | $1.14 \times 10^{-5}$ |
| 0.6 | 1.6539764 | $5.04 \times 10^{-3}$ | 1.6493272 | $3.87 \times 10^{-4}$ |
| 0.8 | 2.1365612 | $9.33 \times 10^{-3}$ | 2.1282568 | $1.03 \times 10^{-3}$ |
| 1.0 | 2.6561319 | $1.53 \times 10^{-2}$ | 2.6427743 | $1.92 \times 10^{-3}$ |
| 1.2 | 3.2033154 | $2.34 \times 10^{-2}$ | 3.1830799 | $3.14 \times 10^{-3}$ |
| 1.4 | 3.7666968 | $3.43 \times 10^{-2}$ | 3.7372087 | $4.81 \times 10^{-3}$ |
| 1.6 | 4.3323743 | $4.89 \times 10^{-2}$ | 4.2905487 | $7.07 \times 10^{-3}$ |
| 1.8 | 4.8834169 | $6.82 \times 10^{-2}$ | 4.8252600 | $1.01 \times 10^{-2}$ |
| 2.0 | 5.3992045 | $9.37 \times 10^{-2}$ | 5.3195640 | $1.41 \times 10^{-2}$ |

\* starting value computed with RK4.

**Convergence.** Global error at $x = 2$:

| $h$ | AB2 $\lvert e \rvert$ | ratio | AB3 $\lvert e \rvert$ | ratio |
|---|---|---|---|---|
| 0.2 | $9.373255 \times 10^{-2}$ | — | $1.409209 \times 10^{-2}$ | — |
| 0.1 | $2.718050 \times 10^{-2}$ | 3.45 | $2.235020 \times 10^{-3}$ | 6.31 |
| 0.05 | $7.257185 \times 10^{-3}$ | 3.75 | $3.117684 \times 10^{-4}$ | 7.17 |
| 0.025 | $1.870205 \times 10^{-3}$ | 3.88 | $4.109527 \times 10^{-5}$ | 7.59 |

The ratios approach $2^2 = 4$ and $2^3 = 8$ more slowly than for one-step methods (the error constants are larger and the asymptotic regime starts at smaller $h$), but the trend confirms orders 2 and 3.

**Cost comparison.** At $h = 0.1$, AB2 uses about 20 evaluations of $f$ (plus 4 for the RK4 starter) and reaches $|e(2)| = 2.72 \times 10^{-2}$; Heun with the same budget ($h = 0.2$) gives $7.24 \times 10^{-2}$.

---

# Why Adams–Bashforth Methods Matter

- They introduce the **multistep** idea: reuse past information instead of recomputing it. Order grows without increasing the cost per step.
- They are built from **polynomial interpolation + quadrature**, tying the ODE module to BRUNO's interpolation and integration modules.
- They are the **predictors** in Adams predictor–corrector schemes (`theory_adams-moulton-2.md`).
- In BRUNO they impose specific Validator rules (uniform mesh, $N \ge k$) and a dependency on a one-step starter, which is why RK4 should be an importable pure function in `ode_ivp`.
- In the catalog, AB2 and AB3 can share one executor parametrized by the coefficient vector $(b_0, b_1, \dots)$ and the number of steps $k$.

---

# Sources

- Burden, R. L., Faires, J. D., & Burden, A. M. (2016). *Numerical Analysis* (10th ed.). Cengage Learning. Section 5.6: Multistep Methods.
- Chapra, S. C., & Canale, R. P. (2015). *Numerical Methods for Engineers* (7th ed.). McGraw-Hill. Section 26.2: Multistep Methods.
- Hairer, E., Nørsett, S. P., & Wanner, G. (1993). *Solving Ordinary Differential Equations I: Nonstiff Problems* (2nd ed.). Springer. Chapter III.
- Süli, E., & Mayers, D. (2003). *An Introduction to Numerical Analysis*. Cambridge University Press. Section 12.6.
- Bashforth, F., & Adams, J. C. (1883). *An Attempt to Test the Theories of Capillary Action*. Cambridge University Press.