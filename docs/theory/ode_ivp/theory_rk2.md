
---

# Runge–Kutta Order 2 (Midpoint Method) — Theory

---

# **Table of Contents**

- [Overview](#overview)
- [Mathematical Foundations](#mathematical-foundations)
  - [Runge–Kutta Methods](#rungekutta-methods)
  - [Derivation of the Order Conditions](#derivation-of-the-order-conditions)
  - [The RK2 Family](#the-rk2-family)
  - [Derivation from the Midpoint Rule](#derivation-from-the-midpoint-rule)
- [Midpoint Method](#midpoint-method)
  - [Algorithm](#algorithm)
  - [Systems of Equations](#systems-of-equations)
- [Error Analysis](#error-analysis)
- [Stability](#stability)
- [Properties](#properties)
- [Example](#example)
- [Why RK2 Matters](#why-rk2-matters)
- [Sources](#sources)

---

# Overview

Let $f(x,y)$ be continuous and Lipschitz in $y$ on $[a,b] \times \mathbb{R}$. The goal is to approximate the solution of

$$
y' = f(x,y), \qquad a \le x \le b, \qquad y(a) = \alpha
$$

at the mesh points $x_i = a + ih$, $i = 0, \dots, N$.

**Runge–Kutta (RK) methods** raise the order of Euler's method by evaluating $f$ at several intermediate points inside each step and combining those slopes. They achieve the accuracy of a high-order Taylor method **without computing derivatives of $f$**.

"RK2" denotes any explicit two-stage RK method of order 2. There is a whole one-parameter family of them. In BRUNO, `rk2` is the **midpoint method**; the **Heun** member of the family has its own document (`theory_heun.md`).

---

# Mathematical Foundations

## Runge–Kutta Methods

An explicit $s$-stage Runge–Kutta method has the form

$$
\begin{aligned}
k_1 &= f(x_i,\; w_i) \\
k_j &= f\!\left(x_i + c_j h,\; w_i + h \sum_{l=1}^{j-1} a_{jl}\, k_l\right), \qquad j = 2, \dots, s \\
w_{i+1} &= w_i + h \sum_{j=1}^{s} b_j\, k_j
\end{aligned}
$$

The coefficients are collected in a **Butcher tableau**:

$$
\begin{array}{c|cccc}
0 & & & & \\
c_2 & a_{21} & & & \\
\vdots & \vdots & \ddots & & \\
c_s & a_{s1} & \cdots & a_{s,s-1} & \\
\hline
 & b_1 & \cdots & b_{s-1} & b_s
\end{array}
$$

with the (customary) row-sum condition $c_j = \sum_l a_{jl}$. Euler is the 1-stage tableau with $b_1 = 1$.

The tableau completely defines the method. This is important for implementation: a **single generic RK routine** driven by tableaus can run Euler, Heun, midpoint, Ralston and RK4.

## Derivation of the Order Conditions

For two stages, write

$$
k_1 = f(x_i, w_i), \qquad k_2 = f(x_i + c\,h,\; w_i + c\,h\,k_1), \qquad w_{i+1} = w_i + h\,(b_1 k_1 + b_2 k_2)
$$

Expanding $k_2$ about $(x_i, y_i)$, with $f, f_x, f_y$ evaluated there:

$$
k_2 = f + c\,h\,(f_x + f\,f_y) + O(h^2)
$$

so the numerical increment is

$$
\phi = (b_1 + b_2)\, f + b_2\, c\, h\,(f_x + f\,f_y) + O(h^2)
$$

The exact solution satisfies

$$
\frac{y(x_{i+1}) - y(x_i)}{h} = f + \frac{h}{2}\,(f_x + f\,f_y) + O(h^2)
$$

Matching terms gives the **order conditions**:

$$
b_1 + b_2 = 1, \qquad b_2\, c = \frac{1}{2}
$$

## The RK2 Family

Choosing $c \neq 0$ as free parameter, $b_2 = \frac{1}{2c}$ and $b_1 = 1 - \frac{1}{2c}$:

| Method | $c$ | $b_1$ | $b_2$ | In BRUNO |
|---|---|---|---|---|
| Midpoint | $\tfrac12$ | $0$ | $1$ | `rk2` |
| Heun | $1$ | $\tfrac12$ | $\tfrac12$ | `heun` |
| Ralston | $\tfrac23$ | $\tfrac14$ | $\tfrac34$ | — |

Midpoint tableau:

$$
\begin{array}{c|cc}
0 & & \\
\tfrac12 & \tfrac12 & \\
\hline
 & 0 & 1
\end{array}
$$

No choice of $c$ can also match the $O(h^2)$ terms for every $f$, so **two stages can never give order 3**.

## Derivation from the Midpoint Rule

Integrating the ODE over one step:

$$
y(x_{i+1}) = y(x_i) + \int_{x_i}^{x_{i+1}} f\big(x, y(x)\big)\,dx
$$

The **midpoint quadrature rule** approximates the integral by $h\, f\big(x_{i+1/2}, y(x_{i+1/2})\big)$, which is exact for linear integrands. The value $y(x_{i+1/2})$ is unknown, so it is estimated with half an Euler step:

$$
y(x_{i+1/2}) \approx w_i + \frac{h}{2}\, f(x_i, w_i)
$$

Geometrically: take a half step along the initial tangent, measure the slope **there**, go back to $(x_i, w_i)$ and take the full step using that midpoint slope.

---

# Midpoint Method

$$
\boxed{
\begin{aligned}
w_0 &= \alpha \\
k_1 &= f(x_i, w_i) \\
k_2 &= f\!\left(x_i + \frac{h}{2},\; w_i + \frac{h}{2}\,k_1\right) \\
w_{i+1} &= w_i + h\,k_2, \qquad i = 0, 1, \dots, N-1
\end{aligned}
}
$$

Note that $k_1$ is **not** used directly in the update; it only serves to locate the midpoint.

## Algorithm

**Input:** $f$, endpoints $a, b$, initial value $\alpha$, number of steps $N$ (or step size $h$).
**Output:** approximations $w_i \approx y(x_i)$ for $i = 0, \dots, N$.

```
h  ← (b − a) / N
w₀ ← α
for i = 0, 1, ..., N−1:
    x_i ← a + i · h
    k1  ← f(x_i, w_i)
    k2  ← f(x_i + h/2, w_i + (h/2) · k1)
    w_{i+1} ← w_i + h · k2
return (x, w)
```

**Cost:** two evaluations of $f$ per step.

## Systems of Equations

Replace scalars by vectors: $\mathbf{k}_1 = \mathbf{f}(x_i, \mathbf{w}_i)$, $\mathbf{k}_2 = \mathbf{f}(x_i + h/2,\ \mathbf{w}_i + \tfrac{h}{2}\mathbf{k}_1)$, $\mathbf{w}_{i+1} = \mathbf{w}_i + h\,\mathbf{k}_2$. The full vector $\mathbf{k}_1$ must be computed before any component of $\mathbf{k}_2$.

---

# Error Analysis

By construction the $O(1)$ and $O(h)$ terms of the local truncation error vanish, so

$$
\tau_{i+1}(h) = O(h^2), \qquad |y(x_i) - w_i| = O(h^2)
$$

The midpoint method has **order 2**. The leading error term is

$$
\tau = h^2\left[\frac{1}{6}\,y''' - \frac{1}{8}\left(f_{xx} + 2f\,f_{xy} + f^2 f_{yy}\right)\right] + O(h^3)
$$

Heun's leading term has the same structure with $\tfrac14$ in place of $\tfrac18$. Which member of the family is more accurate depends on $f$; Ralston's $c = \tfrac23$ minimizes a bound on this term over all $f$.

---

# Stability

On $y' = \lambda y$, $z = h\lambda$:

$$
k_1 = \lambda w_i, \quad k_2 = \lambda\left(1 + \frac{z}{2}\right) w_i, \quad w_{i+1} = \left(1 + z + \frac{z^2}{2}\right) w_i
$$

$$
R(z) = 1 + z + \frac{z^2}{2}
$$

This is the **same** stability function as Heun's method. It is a general fact: every explicit $s$-stage method of order $p = s$ (possible for $s \le 4$) has $R(z) = \sum_{j=0}^{s} z^j / j!$. Real stability interval: $-2 \le h\lambda \le 0$.

---

# Properties

- Explicit, one-step, self-starting.
- Two function evaluations per step.
- Local truncation error $O(h^2)$; global error $O(h^2)$ — **order 2**.
- Uses only the midpoint slope in the update; $k_1$ is an auxiliary predictor.
- Same stability function and region as Heun.
- Requires $y \in C^3[a,b]$ for the error estimate.

---

# Example

Approximate the solution of

$$
y' = y - x^2 + 1, \qquad 0 \le x \le 2, \qquad y(0) = 0.5
$$

with $N = 10$, $h = 0.2$. Exact solution: $y(x) = (x+1)^2 - 0.5\,e^{x}$.

**First step.**

$$
k_1 = f(0,\, 0.5) = 1.5
$$

$$
k_2 = f(0.1,\; 0.5 + 0.1 \cdot 1.5) = f(0.1,\, 0.65) = 0.65 - 0.01 + 1 = 1.64
$$

$$
w_1 = 0.5 + 0.2 \cdot 1.64 = 0.828
$$

**Full table.**

| $i$ | $x_i$ | $k_1$ | $k_2$ | $w_i$ | $y(x_i)$ | $\lvert y(x_i) - w_i \rvert$ |
|---|---|---|---|---|---|---|
| 0 | 0.0 | — | — | 0.5000000 | 0.5000000 | 0.0000000 |
| 1 | 0.2 | 1.5000000 | 1.6400000 | 0.8280000 | 0.8292986 | 0.0012986 |
| 2 | 0.4 | 1.7880000 | 1.9168000 | 1.2113600 | 1.2140877 | 0.0027277 |
| 3 | 0.6 | 2.0513600 | 2.1664960 | 1.6446592 | 1.6489406 | 0.0042814 |
| 4 | 0.8 | 2.2846592 | 2.3831251 | 2.1212842 | 2.1272295 | 0.0059453 |
| 5 | 1.0 | 2.4812842 | 2.5594126 | 2.6331668 | 2.6408591 | 0.0076923 |
| 6 | 1.2 | 2.6331668 | 2.6864834 | 3.1704634 | 3.1799415 | 0.0094781 |
| 7 | 1.4 | 2.7304634 | 2.7535098 | 3.7211654 | 3.7324000 | 0.0112346 |
| 8 | 1.6 | 2.7611654 | 2.7472819 | 4.2706218 | 4.2834838 | 0.0128620 |
| 9 | 1.8 | 2.7106218 | 2.6516840 | 4.8009586 | 4.8151763 | 0.0142177 |
| 10 | 2.0 | 2.5609586 | 2.4470544 | 5.2903695 | 5.3054720 | 0.0151025 |

**Comparison within the RK2 family** at $x = 2$, $h = 0.2$:

| Method | $\lvert e(2) \rvert$ |
|---|---|
| Heun | 0.0724173 |
| Midpoint | 0.0151025 |

Same order and cost, but for this $f$ the midpoint method is almost 5 times more accurate. The error *constant* depends on the method and the problem; the *rate* does not.

**Convergence.** Global error at $x = 2$:

| $h$ | $\lvert e \rvert$ | ratio |
|---|---|---|
| 0.2 | $1.510249 \times 10^{-2}$ | — |
| 0.1 | $3.747074 \times 10^{-3}$ | 4.03 |
| 0.05 | $9.277142 \times 10^{-4}$ | 4.04 |
| 0.025 | $2.304037 \times 10^{-4}$ | 4.03 |

Ratios ≈ 4 = $2^2$: order 2.

---

# Why RK2 Matters

- It is the smallest example where the **Runge–Kutta order conditions** can be derived by hand, and it shows why they leave free parameters.
- It introduces the **Butcher tableau**, the data structure that fully describes an RK method.
- In BRUNO, the tableau view suggests a catalog-driven design for `ode_ivp`: one generic explicit-RK executor, and Euler, Heun, midpoint, Ralston and RK4 registered as tableaus in the JSON catalog.
- It shows that two methods with the same order, cost and stability can still differ several-fold in accuracy.

---

# Sources

- Burden, R. L., Faires, J. D., & Burden, A. M. (2016). *Numerical Analysis* (10th ed.). Cengage Learning. Section 5.4: Runge–Kutta Methods.
- Chapra, S. C., & Canale, R. P. (2015). *Numerical Methods for Engineers* (7th ed.). McGraw-Hill. Sections 25.2–25.3.
- Hairer, E., Nørsett, S. P., & Wanner, G. (1993). *Solving Ordinary Differential Equations I: Nonstiff Problems* (2nd ed.). Springer. Section II.1.
- Butcher, J. C. (2016). *Numerical Methods for Ordinary Differential Equations* (3rd ed.). Wiley.
- Ralston, A. (1962). Runge–Kutta methods with minimum error bounds. *Mathematics of Computation*, 16(80), 431–437.