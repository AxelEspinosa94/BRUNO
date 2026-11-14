
---

# Heun's Method — Theory

---

# **Table of Contents**

- [Overview](#overview)
- [Mathematical Foundations](#mathematical-foundations)
  - [Derivation from the Trapezoidal Rule](#derivation-from-the-trapezoidal-rule)
  - [Predictor–Corrector Interpretation](#predictorcorrector-interpretation)
  - [The Second-Order Runge–Kutta Family](#the-second-order-rungekutta-family)
- [Heun's Method](#heuns-method)
  - [Algorithm](#algorithm)
  - [Systems of Equations](#systems-of-equations)
  - [Iterated Heun](#iterated-heun)
- [Error Analysis](#error-analysis)
  - [Local Truncation Error](#local-truncation-error)
  - [Global Error](#global-error)
- [Stability](#stability)
- [Properties](#properties)
- [Example](#example)
- [Why Heun's Method Matters](#why-heuns-method-matters)
- [Sources](#sources)

---

# Overview

Let $f(x,y)$ be continuous and Lipschitz in $y$ on $[a,b] \times \mathbb{R}$. The goal is to approximate the solution of the initial value problem

$$
y' = f(x,y), \qquad a \le x \le b, \qquad y(a) = \alpha
$$

at the mesh points $x_i = a + ih$, $i = 0, 1, \dots, N$.

Euler's method uses only the slope at the **start** of each step. **Heun's method** (also called the *improved Euler* or *explicit trapezoidal* method) fixes this by averaging the slope at the start with an estimate of the slope at the **end** of the step. The extra function evaluation raises the order from 1 to 2.

> **Naming warning.** The literature is inconsistent:
>
> | Source | What it calls "Heun" |
> |---|---|
> | Chapra & Canale; most engineering texts | this method (order 2) |
> | Burden & Faires | a **third-order** RK method; they call this method *Modified Euler* |
> | Hairer, Nørsett & Wanner | this method (explicit trapezoidal rule) |
>
> In BRUNO, `heun` is the order-2 method described here, and `rk2` is the **midpoint** method. Both are members of the same RK2 family (see below).

See `theory_ivp.md` for the definitions of local/global error, order and stability, and `theory_euler.md` for the base method.

---

# Mathematical Foundations

## Derivation from the Trapezoidal Rule

Integrating the ODE over one step gives the exact relation

$$
y(x_{i+1}) = y(x_i) + \int_{x_i}^{x_{i+1}} f\big(x, y(x)\big)\,dx
$$

Euler approximates the integral with the left-endpoint rectangle. Using the **trapezoidal rule** instead,

$$
\int_{x_i}^{x_{i+1}} f\big(x, y(x)\big)\,dx \approx \frac{h}{2}\Big[f\big(x_i, y(x_i)\big) + f\big(x_{i+1}, y(x_{i+1})\big)\Big]
$$

yields the **implicit trapezoidal method**

$$
w_{i+1} = w_i + \frac{h}{2}\Big[f(x_i, w_i) + f(x_{i+1}, w_{i+1})\Big]
$$

The unknown $w_{i+1}$ appears on both sides, so for nonlinear $f$ each step requires solving an equation.

## Predictor–Corrector Interpretation

Heun avoids the implicit equation by replacing $w_{i+1}$ on the right-hand side with an Euler estimate:

1. **Predictor** (Euler):
$$
\tilde{w}_{i+1} = w_i + h\,f(x_i, w_i)
$$
2. **Corrector** (trapezoid with the predicted value):
$$
w_{i+1} = w_i + \frac{h}{2}\Big[f(x_i, w_i) + f(x_{i+1}, \tilde{w}_{i+1})\Big]
$$

Geometrically: take the slope at the start, take the slope at the point Euler would land on, and step with their **average**. For a convex solution the start slope underestimates and the end slope overestimates, so the errors largely cancel.

## The Second-Order Runge–Kutta Family

Consider the general two-stage explicit method

$$
\begin{aligned}
k_1 &= f(x_i, w_i) \\
k_2 &= f(x_i + c\,h,\; w_i + c\,h\,k_1) \\
w_{i+1} &= w_i + h\,(b_1 k_1 + b_2 k_2)
\end{aligned}
$$

Expanding $k_2$ in a Taylor series about $(x_i, w_i)$ and matching with the expansion of the exact solution,

$$
y(x_{i+1}) = y + h f + \frac{h^2}{2}\left(f_x + f\,f_y\right) + O(h^3)
$$

gives the **order conditions**

$$
b_1 + b_2 = 1, \qquad b_2\,c = \frac{1}{2}
$$

Two equations, three unknowns: there is a **one-parameter family** of order-2 methods. In Butcher tableau form:

| Method | $c$ | $b_1$ | $b_2$ | In BRUNO |
|---|---|---|---|---|
| Heun (explicit trapezoid) | $1$ | $\tfrac12$ | $\tfrac12$ | `heun` |
| Midpoint | $\tfrac12$ | $0$ | $1$ | `rk2` |
| Ralston | $\tfrac23$ | $\tfrac14$ | $\tfrac34$ | — |

All three have the same order. Ralston's choice minimizes the leading error constant; Heun's is the one with the clearest geometric meaning.

---

# Heun's Method

$$
\boxed{
\begin{aligned}
w_0 &= \alpha \\
k_1 &= f(x_i, w_i) \\
k_2 &= f(x_i + h,\; w_i + h\,k_1) \\
w_{i+1} &= w_i + \frac{h}{2}\,(k_1 + k_2), \qquad i = 0, 1, \dots, N-1
\end{aligned}
}
$$

Increment function:

$$
\phi(x, w; h) = \frac{1}{2}\Big[f(x, w) + f\big(x + h,\; w + h\,f(x,w)\big)\Big]
$$

Note that $\phi(x, w; 0) = f(x,w)$, so the method is consistent.

## Algorithm

**Input:** $f$, endpoints $a, b$, initial value $\alpha$, number of steps $N$ (or step size $h$).
**Output:** approximations $w_i \approx y(x_i)$ for $i = 0, \dots, N$.

```
h  ← (b − a) / N
x₀ ← a
w₀ ← α
for i = 0, 1, ..., N−1:
    k1 ← f(x_i, w_i)
    k2 ← f(x_i + h, w_i + h · k1)
    w_{i+1} ← w_i + (h / 2) · (k1 + k2)
    x_{i+1} ← a + (i+1) · h
return (x, w)
```

**Cost:** two evaluations of $f$ per step, $2N$ in total. Since $k_2$ is evaluated at $x_{i+1}$ but at the *predicted* value, it **cannot** be reused as $k_1$ of the next step (unlike in some embedded RK pairs).

Exposing $k_1$ and $k_2$ in the output table makes intermediate results auditable, which is useful for teaching and for debugging.

## Systems of Equations

For $\mathbf{y}' = \mathbf{f}(x, \mathbf{y})$ with $\mathbf{y} \in \mathbb{R}^m$, each stage becomes a vector:

$$
\begin{aligned}
\mathbf{k}_1 &= \mathbf{f}(x_i, \mathbf{w}_i) \\
\mathbf{k}_2 &= \mathbf{f}(x_i + h,\; \mathbf{w}_i + h\,\mathbf{k}_1) \\
\mathbf{w}_{i+1} &= \mathbf{w}_i + \frac{h}{2}\,(\mathbf{k}_1 + \mathbf{k}_2)
\end{aligned}
$$

$\mathbf{k}_1$ must be computed **completely** (all $m$ components) before any component of $\mathbf{k}_2$, since every component of $\mathbf{k}_2$ depends on the full predicted vector.

## Iterated Heun

Some texts (e.g. Chapra & Canale) apply the corrector repeatedly:

$$
w_{i+1}^{(j+1)} = w_i + \frac{h}{2}\Big[f(x_i, w_i) + f\big(x_{i+1}, w_{i+1}^{(j)}\big)\Big], \qquad w_{i+1}^{(0)} = \tilde{w}_{i+1}
$$

until $|w_{i+1}^{(j+1)} - w_{i+1}^{(j)}|$ falls below a tolerance. The iteration converges, for $h L / 2 < 1$, to the **implicit trapezoidal** solution, **not** to $y(x_{i+1})$. It does not raise the order beyond 2, and it is not necessarily more accurate than a single correction. The standard (non-iterated) form is what BRUNO implements.

---

# Error Analysis

## Local Truncation Error

Write $f$, $f_x$, $f_y$ for values at $(x_i, y(x_i))$. Expanding $k_2$:

$$
k_2 = f + h\left(f_x + f\,f_y\right) + \frac{h^2}{2}\left(f_{xx} + 2 f\,f_{xy} + f^2 f_{yy}\right) + O(h^3)
$$

so

$$
\phi = f + \frac{h}{2}\left(f_x + f\,f_y\right) + \frac{h^2}{4}\left(f_{xx} + 2 f\,f_{xy} + f^2 f_{yy}\right) + O(h^3)
$$

The exact increment is

$$
\frac{y(x_{i+1}) - y(x_i)}{h} = y' + \frac{h}{2}\,y'' + \frac{h^2}{6}\,y''' + O(h^3)
$$

with $y' = f$ and $y'' = f_x + f f_y$. The $O(1)$ and $O(h)$ terms match exactly, so

$$
\tau_{i+1}(h) = O(h^2)
$$

The leading term involves $y'''$ and second partial derivatives of $f$, so the method needs $y \in C^3[a,b]$ for this estimate.

## Global Error

By the general theory for one-step methods with a Lipschitz increment function (`theory_ivp.md`):

$$
|y(x_i) - w_i| = O(h^2)
$$

Heun's method has **order 2**: halving $h$ divides the error by approximately 4.

Since a step costs two $f$-evaluations, the fair comparison with Euler is at **equal work**: Heun with step $h$ against Euler with step $h/2$. Even then Heun wins clearly whenever the solution is smooth, because $O(h^2)$ beats $O(h)$.

---

# Stability

On the test equation $y' = \lambda y$, with $z = h\lambda$:

$$
k_1 = \lambda w_i, \qquad k_2 = \lambda\,(w_i + h\lambda w_i) = \lambda(1 + z)\,w_i
$$

$$
w_{i+1} = w_i + \frac{h}{2}\big(\lambda w_i + \lambda(1+z) w_i\big) = \left(1 + z + \frac{z^2}{2}\right) w_i
$$

The stability function

$$
R(z) = 1 + z + \frac{z^2}{2}
$$

is the degree-2 Taylor polynomial of $e^{z}$. **Every** explicit 2-stage, order-2 Runge–Kutta method (Heun, midpoint, Ralston) has this same $R(z)$, and therefore the same stability region.

On the negative real axis, $|R(z)| \le 1$ for

$$
-2 \le z \le 0 \quad\Longleftrightarrow\quad h \le \frac{2}{|\lambda|}
$$

the same interval as Euler. The region is slightly larger than Euler's disk in the complex plane, but it still contains no part of the imaginary axis except the origin: oscillatory problems are amplified, though far more slowly than with Euler since $|R(i\omega h)|^2 = 1 + (\omega h)^4/4$.

---

# Properties

- Explicit, one-step, self-starting.
- Two function evaluations per step.
- Local truncation error $O(h^2)$; global error $O(h^2)$ — **order 2**.
- Equivalent to the trapezoidal rule with an Euler predictor.
- Member of the RK2 family ($c = 1$, $b_1 = b_2 = \tfrac12$).
- Stability function $R(z) = 1 + z + z^2/2$; real interval $-2 \le h\lambda \le 0$.
- Requires $y \in C^3[a,b]$ for the error estimate.

---

# Example

Approximate the solution of

$$
y' = y - x^2 + 1, \qquad 0 \le x \le 2, \qquad y(0) = 0.5
$$

with $N = 10$, $h = 0.2$. The exact solution is $y(x) = (x+1)^2 - 0.5\,e^{x}$.

**First step.**

$$
k_1 = f(0,\, 0.5) = 0.5 - 0 + 1 = 1.5
$$

$$
k_2 = f(0.2,\; 0.5 + 0.2 \cdot 1.5) = f(0.2,\, 0.8) = 0.8 - 0.04 + 1 = 1.76
$$

$$
w_1 = 0.5 + \frac{0.2}{2}\,(1.5 + 1.76) = 0.5 + 0.1 \cdot 3.26 = 0.826
$$

Compare with Euler's $w_1 = 0.8$ and the exact $y(0.2) = 0.8292986$: the error drops from $0.0293$ to $0.0033$.

**Full table.**

| $i$ | $x_i$ | $k_1$ | $k_2$ | $w_i$ | $y(x_i)$ | $\lvert y(x_i) - w_i \rvert$ |
|---|---|---|---|---|---|---|
| 0 | 0.0 | — | — | 0.5000000 | 0.5000000 | 0.0000000 |
| 1 | 0.2 | 1.5000000 | 1.7600000 | 0.8260000 | 0.8292986 | 0.0032986 |
| 2 | 0.4 | 1.7860000 | 2.0232000 | 1.2069200 | 1.2140877 | 0.0071677 |
| 3 | 0.6 | 2.0469200 | 2.2563040 | 1.6372424 | 1.6489406 | 0.0116982 |
| 4 | 0.8 | 2.2772424 | 2.4526909 | 2.1102357 | 2.1272295 | 0.0169938 |
| 5 | 1.0 | 2.4702357 | 2.6042829 | 2.6176876 | 2.6408591 | 0.0231715 |
| 6 | 1.2 | 2.6176876 | 2.7012251 | 3.1495789 | 3.1799415 | 0.0303627 |
| 7 | 1.4 | 2.7095789 | 2.7314946 | 3.6936862 | 3.7324000 | 0.0387138 |
| 8 | 1.6 | 2.7336862 | 2.6804234 | 4.2350972 | 4.2834838 | 0.0483866 |
| 9 | 1.8 | 2.6750972 | 2.5301166 | 4.7556185 | 4.8151763 | 0.0595577 |
| 10 | 2.0 | 2.5156185 | 2.2587423 | 5.2330546 | 5.3054720 | 0.0724173 |

Here $k_1, k_2$ in row $i$ are the slopes used to go from $w_{i-1}$ to $w_i$.

**Comparison with Euler at $x = 2$.**

| Method | $h$ | $f$-evaluations | $\lvert e(2) \rvert$ |
|---|---|---|---|
| Euler | 0.2 | 10 | 0.4396874 |
| Euler | 0.1 | 20 | 0.2419719 |
| Heun | 0.2 | 20 | 0.0724173 |

At equal cost (20 evaluations), Heun is about 3.3 times more accurate than Euler.

**Convergence.** At $x = 2$: $h = 0.1$ gives $|e| = 1.890478 \times 10^{-2}$ and $h = 0.05$ gives $|e| = 4.819865 \times 10^{-3}$. Ratios 3.83 and 3.92 approach $2^2 = 4$, confirming order 2.

---

# Why Heun's Method Matters

- It is the simplest example of how **sampling more slopes per step** raises the order: the core idea of all Runge–Kutta methods.
- Its derivation connects ODE solvers to quadrature (trapezoidal rule) and to **predictor–corrector** schemes, which reappear in the Adams–Bashforth/Adams–Moulton pairs.
- It shows that order alone does not determine stability: Heun doubles Euler's order but keeps the same real stability interval.
- In BRUNO, `heun` and `rk2` are two instances of the **same** two-stage formula with different coefficients $(c, b_1, b_2)$. The `ode_ivp` executor can implement one generic RK2 routine parametrized by its Butcher tableau, and register Heun, midpoint and Ralston as catalog entries rather than separate code paths.
- The order-2 convergence ratio (≈ 4) gives a regression test that distinguishes a correct Heun from an accidental Euler (≈ 2).

---

# Sources

- Burden, R. L., Faires, J. D., & Burden, A. M. (2016). *Numerical Analysis* (10th ed.). Cengage Learning. Section 5.4: Runge–Kutta Methods (Modified Euler method).
- Chapra, S. C., & Canale, R. P. (2015). *Numerical Methods for Engineers* (7th ed.). McGraw-Hill. Section 25.2: Improvements of Euler's Method (Heun's method).
- Hairer, E., Nørsett, S. P., & Wanner, G. (1993). *Solving Ordinary Differential Equations I: Nonstiff Problems* (2nd ed.). Springer. Section II.1.
- Süli, E., & Mayers, D. (2003). *An Introduction to Numerical Analysis*. Cambridge University Press. Section 12.5.
- Butcher, J. C. (2016). *Numerical Methods for Ordinary Differential Equations* (3rd ed.). Wiley.