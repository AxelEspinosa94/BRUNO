
---

# Euler's Method — Theory

---

# **Table of Contents**

- [Overview](#overview)
- [Mathematical Foundations](#mathematical-foundations)
  - [Derivation from Taylor's Theorem](#derivation-from-taylors-theorem)
  - [Geometric Interpretation](#geometric-interpretation)
- [Euler's Method](#eulers-method)
  - [Algorithm](#algorithm)
  - [Systems of Equations](#systems-of-equations)
- [Error Analysis](#error-analysis)
  - [Local Truncation Error](#local-truncation-error)
  - [Global Error Bound](#global-error-bound)
- [Stability](#stability)
- [Properties](#properties)
- [Example](#example)
- [Why Euler's Method Matters](#why-eulers-method-matters)
- [Sources](#sources)

---

# Overview

Let $f(x,y)$ be continuous and Lipschitz in $y$ on $[a,b] \times \mathbb{R}$. The goal is to approximate the solution of the initial value problem

$$
y' = f(x,y), \qquad a \le x \le b, \qquad y(a) = \alpha
$$

at the mesh points $x_i = a + ih$, $i = 0, 1, \dots, N$, with $h = (b-a)/N$.

**Euler's method** is the simplest numerical method for IVPs: at each step it follows the tangent line of the solution for a distance $h$. It is rarely the method of choice in practice because of its low accuracy, but every other one-step method can be understood as an improvement on it.

See `theory_ivp.md` for the definitions of local/global error, order and stability used below.

---

# Mathematical Foundations

## Derivation from Taylor's Theorem

Assume the exact solution $y(x)$ has two continuous derivatives on $[a,b]$. Expanding $y$ about $x_i$:

$$
y(x_{i+1}) = y(x_i) + h\,y'(x_i) + \frac{h^2}{2}\,y''(\xi_i), \qquad \xi_i \in (x_i, x_{i+1})
$$

Since $y$ solves the differential equation, $y'(x_i) = f\big(x_i, y(x_i)\big)$:

$$
y(x_{i+1}) = y(x_i) + h\,f\big(x_i, y(x_i)\big) + \frac{h^2}{2}\,y''(\xi_i)
$$

Dropping the remainder term and replacing $y(x_i)$ by the approximation $w_i$ gives Euler's method.

The same formula can be obtained by integrating the ODE over one step,

$$
y(x_{i+1}) = y(x_i) + \int_{x_i}^{x_{i+1}} f\big(x, y(x)\big)\,dx
$$

and approximating the integral with the **left-endpoint rectangle rule**. This link to quadrature is the bridge to Heun's method (trapezoidal rule) and the Adams family (interpolating polynomials).

## Geometric Interpretation

At $(x_i, w_i)$ the direction field prescribes slope $f(x_i, w_i)$. Euler's method moves along that straight line until $x_{i+1}$:

$$
w_{i+1} = w_i + \underbrace{h}_{\text{run}} \cdot \underbrace{f(x_i, w_i)}_{\text{slope}}
$$

Because the slope is sampled only at the **start** of the interval, the method systematically lags behind curving solutions: for a convex solution ($y'' > 0$) it underestimates, for a concave one it overestimates.

After the first step, $w_i$ no longer lies on the true solution curve, so subsequent slopes are taken from a *neighboring* solution of the ODE. This is how local errors propagate into global error.

---

# Euler's Method

$$
\boxed{
\begin{aligned}
w_0 &= \alpha \\
w_{i+1} &= w_i + h\,f(x_i, w_i), \qquad i = 0, 1, \dots, N-1
\end{aligned}
}
$$

In the one-step notation of `theory_ivp.md`, the increment function is simply

$$
\phi(x, w; h) = f(x, w)
$$

## Algorithm

**Input:** $f$, endpoints $a, b$, initial value $\alpha$, number of steps $N$ (or step size $h$).
**Output:** approximations $w_i \approx y(x_i)$ for $i = 0, \dots, N$.

```
h  ← (b − a) / N
x₀ ← a
w₀ ← α
for i = 0, 1, ..., N−1:
    w_{i+1} ← w_i + h · f(x_i, w_i)
    x_{i+1} ← a + (i+1) · h
return (x, w)
```

Computing $x_{i+1} = a + (i+1)h$ instead of accumulating $x_{i+1} = x_i + h$ avoids drift from repeated floating-point addition, so the final mesh point lands on $b$ as exactly as possible.

**Cost:** one evaluation of $f$ per step, $N$ evaluations in total.

## Systems of Equations

For $\mathbf{y}' = \mathbf{f}(x, \mathbf{y})$, $\mathbf{y}(a) = \boldsymbol{\alpha} \in \mathbb{R}^m$, the method is identical in vector form:

$$
\mathbf{w}_{i+1} = \mathbf{w}_i + h\,\mathbf{f}(x_i, \mathbf{w}_i)
$$

All components must be updated using the slopes evaluated at the **same** old state $\mathbf{w}_i$. Updating $w^{(1)}$ first and then using the new value to compute the slope of $w^{(2)}$ is a different (and generally inconsistent) scheme.

---

# Error Analysis

## Local Truncation Error

From the Taylor expansion above:

$$
\tau_{i+1}(h) = \frac{y(x_{i+1}) - y(x_i)}{h} - f\big(x_i, y(x_i)\big) = \frac{h}{2}\,y''(\xi_i)
$$

If $|y''(x)| \le M$ on $[a,b]$, then

$$
|\tau_{i+1}(h)| \le \frac{hM}{2} = O(h)
$$

The error introduced in a single step is $h\,\tau_{i+1} = O(h^2)$.

## Global Error Bound

**Theorem.** Suppose $f$ is continuous and satisfies a Lipschitz condition with constant $L$ on $D = \{(x,y): a \le x \le b,\ -\infty < y < \infty\}$, and that $|y''(x)| \le M$ on $[a,b]$. Then for each $i = 0, 1, \dots, N$:

$$
|y(x_i) - w_i| \le \frac{hM}{2L}\left[e^{L(x_i - a)} - 1\right]
$$

Consequences:

- The global error is $O(h)$: Euler's method has **order 1**.
- Halving $h$ roughly halves the error.
- The bound grows **exponentially** with the length of the interval through $e^{L(x_i - a)}$. It is usually pessimistic, but it shows why long integrations with Euler accumulate large errors.

The bound needs $M$, which in turn needs the unknown solution. In practice $M$ is estimated from $y'' = \frac{\partial f}{\partial x} + \frac{\partial f}{\partial y}\,f$, or the error is estimated empirically by repeating the computation with $h/2$.

**With rounding.** If each step carries a rounding error bounded by $\delta$, the bound becomes

$$
|y(x_i) - w_i| \le \frac{1}{L}\left(\frac{hM}{2} + \frac{\delta}{h}\right)\left[e^{L(x_i - a)} - 1\right] + |\delta_0|\,e^{L(x_i - a)}
$$

which is minimized at $h = \sqrt{2\delta/M}$. Taking $h$ smaller than that makes the result worse, not better.

---

# Stability

On the test equation $y' = \lambda y$:

$$
w_{i+1} = w_i + h\lambda\,w_i = (1 + z)\,w_i, \qquad z = h\lambda
$$

The stability function is $R(z) = 1 + z$, and the region of absolute stability is the disk

$$
\mathcal{S} = \{ z \in \mathbb{C} : |1 + z| \le 1 \}
$$

of radius 1 centered at $-1$. For real $\lambda < 0$ this requires

$$
h \le \frac{2}{|\lambda|}
$$

and for monotone (non-oscillating) decay, $h \le 1/|\lambda|$. Purely imaginary $\lambda$ (undamped oscillators) lie outside $\mathcal{S}$ for every $h > 0$: Euler always makes the amplitude of an oscillator grow.

---

# Properties

- Explicit, one-step, self-starting.
- One function evaluation per step.
- Local truncation error $O(h)$; global error $O(h)$ — **order 1**.
- Equivalent to the left-endpoint rectangle rule applied to $\int f\,dx$.
- Stability interval on the real axis: $-2 \le h\lambda \le 0$.
- Unconditionally unstable for purely oscillatory problems.
- Requires $y \in C^2[a,b]$ for the error bound to hold.

---

# Example

Approximate the solution of

$$
y' = y - x^2 + 1, \qquad 0 \le x \le 2, \qquad y(0) = 0.5
$$

with $N = 10$, so $h = 0.2$. The exact solution is $y(x) = (x+1)^2 - 0.5\,e^{x}$.

**First step.**

$$
w_1 = w_0 + h\,f(x_0, w_0) = 0.5 + 0.2\,(0.5 - 0^2 + 1) = 0.5 + 0.2\,(1.5) = 0.8
$$

**Second step.**

$$
w_2 = w_1 + h\,f(x_1, w_1) = 0.8 + 0.2\,(0.8 - 0.04 + 1) = 0.8 + 0.2\,(1.76) = 1.152
$$

**Full table.**

| $i$ | $x_i$ | $w_i$ | $y(x_i)$ | $\lvert y(x_i) - w_i \rvert$ |
|---|---|---|---|---|
| 0 | 0.0 | 0.5000000 | 0.5000000 | 0.0000000 |
| 1 | 0.2 | 0.8000000 | 0.8292986 | 0.0292986 |
| 2 | 0.4 | 1.1520000 | 1.2140877 | 0.0620877 |
| 3 | 0.6 | 1.5504000 | 1.6489406 | 0.0985406 |
| 4 | 0.8 | 1.9884800 | 2.1272295 | 0.1387495 |
| 5 | 1.0 | 2.4581760 | 2.6408591 | 0.1826831 |
| 6 | 1.2 | 2.9498112 | 3.1799415 | 0.2301303 |
| 7 | 1.4 | 3.4517734 | 3.7324000 | 0.2806266 |
| 8 | 1.6 | 3.9501281 | 4.2834838 | 0.3333557 |
| 9 | 1.8 | 4.4281538 | 4.8151763 | 0.3870225 |
| 10 | 2.0 | 4.8657845 | 5.3054720 | 0.4396874 |

**Error bound.** Here $L = 1$ and $y''(x) = 2 - 0.5\,e^{x}$, so

$$
M = \max_{[0,2]} |y''(x)| = 0.5\,e^{2} - 2 = 1.6945280
$$

$$
|y(x_i) - w_i| \le \frac{0.2 \cdot 1.6945280}{2}\left[e^{x_i} - 1\right] = 0.1694528\left[e^{x_i} - 1\right]
$$

| $x_i$ | bound | actual error |
|---|---|---|
| 1.0 | 0.2911677 | 0.1826831 |
| 2.0 | 1.0826435 | 0.4396874 |

The bound holds and, as expected, is pessimistic. The error grows steadily because the solution is convex on most of the interval and Euler always lags behind it.

**Convergence.** At $x = 2$: $h = 0.1$ gives $|e| = 0.2419719$ and $h = 0.05$ gives $|e| = 0.1274657$. Each halving of $h$ roughly halves the error (ratios 1.82 and 1.90 approaching 2), confirming order 1.

---

# Why Euler's Method Matters

- It is the **baseline** against which every other IVP method is measured.
- Its derivation (Taylor truncation / rectangle quadrature) generalizes directly to Runge–Kutta (more slope samples per step) and Adams methods (more past points).
- Its error and stability analysis is simple enough to do by hand, making it the reference case for understanding local vs global error and stability regions.
- In BRUNO it is the natural **first method to migrate** in `ode_ivp`: it exercises the whole Constructor → Validator → Executor → Renderer pipeline with the minimum numerical complexity.
- Its first-order convergence ratio (≈ 2) is an easy regression test that the pipeline is wired correctly.

---

# Sources

- Burden, R. L., Faires, J. D., & Burden, A. M. (2016). *Numerical Analysis* (10th ed.). Cengage Learning. Section 5.2: Euler's Method.
- Chapra, S. C., & Canale, R. P. (2015). *Numerical Methods for Engineers* (7th ed.). McGraw-Hill. Section 25.1: Euler's Method.
- Hairer, E., Nørsett, S. P., & Wanner, G. (1993). *Solving Ordinary Differential Equations I: Nonstiff Problems* (2nd ed.). Springer. Section I.7.
- Süli, E., & Mayers, D. (2003). *An Introduction to Numerical Analysis*. Cambridge University Press. Section 12.2.