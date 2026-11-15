---

# Classical Runge–Kutta Method (RK4) — Theory

---

# **Table of Contents**

- [Overview](#overview)
- [Mathematical Foundations](#mathematical-foundations)
  - [Order Conditions for Four Stages](#order-conditions-for-four-stages)
  - [Butcher Tableau](#butcher-tableau)
  - [Connection with Simpson's Rule](#connection-with-simpsons-rule)
  - [Butcher Barriers](#butcher-barriers)
- [RK4 Method](#rk4-method)
  - [Algorithm](#algorithm)
- [Error Analysis](#error-analysis)
  - [Local and Global Error](#local-and-global-error)
  - [Error Estimation by Step Doubling](#error-estimation-by-step-doubling)
- [Stability](#stability)
- [Properties](#properties)
- [Example](#example)
- [Why RK4 Matters](#why-rk4-matters)
- [Sources](#sources)

---

# Overview

Let $f(x,y)$ be continuous and Lipschitz in $y$ on $[a,b] \times \mathbb{R}$. The goal is to approximate the solution of

$$
y' = f(x,y), \qquad a \le x \le b, \qquad y(a) = \alpha
$$

at the mesh points $x_i = a + ih$, $i = 0, \dots, N$.

The **classical fourth-order Runge–Kutta method** (RK4) samples the slope four times per step (once at the start, twice at the midpoint and once at the end) and combines them with weights $1, 2, 2, 1$. It reaches **order 4** at a cost of four evaluations of $f$ per step.

RK4 offers one of the best balances between accuracy, simplicity and cost, and is the standard general-purpose fixed-step solver. In BRUNO it is also the natural **starter** for the multistep methods and the inner solver of the **shooting** method.

See `theory_ivp.md` for the general theory and `theory_rk2.md` for the Runge–Kutta framework and Butcher tableaus.

---

# Mathematical Foundations

## Order Conditions for Four Stages

For an explicit RK method to have order $p$, its Taylor expansion must match that of the exact solution up to $h^p$. Each term corresponds to a **rooted tree** (Butcher's theory); the number of conditions grows quickly:

| Order $p$ | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| Conditions | 1 | 2 | 4 | 8 | 17 | 37 |

For order 4 the 8 conditions are (with $c_j = \sum_l a_{jl}$):

$$
\begin{aligned}
&\textstyle\sum b_j = 1, &\quad &\textstyle\sum b_j c_j = \frac12, \\
&\textstyle\sum b_j c_j^2 = \frac13, & &\textstyle\sum_{j,l} b_j a_{jl} c_l = \frac16, \\
&\textstyle\sum b_j c_j^3 = \frac14, & &\textstyle\sum_{j,l} b_j c_j a_{jl} c_l = \frac18, \\
&\textstyle\sum_{j,l} b_j a_{jl} c_l^2 = \frac{1}{12}, & &\textstyle\sum_{j,l,m} b_j a_{jl} a_{lm} c_m = \frac{1}{24}
\end{aligned}
$$

A four-stage explicit method has 10 free coefficients, so there is a two-parameter family of solutions. The classical RK4 is the best-known one; Kutta's **3/8 rule** is another.

## Butcher Tableau

$$
\begin{array}{c|cccc}
0 & & & & \\
\tfrac12 & \tfrac12 & & & \\
\tfrac12 & 0 & \tfrac12 & & \\
1 & 0 & 0 & 1 & \\
\hline
 & \tfrac16 & \tfrac13 & \tfrac13 & \tfrac16
\end{array}
$$

Each stage only uses the immediately previous one, so the method is simple to implement and needs little storage.

## Connection with Simpson's Rule

If $f$ does not depend on $y$, the ODE reduces to the quadrature $y(x_{i+1}) = y(x_i) + \int_{x_i}^{x_{i+1}} f(x)\,dx$. Then $k_2 = k_3 = f(x_{i+1/2})$ and

$$
w_{i+1} = w_i + \frac{h}{6}\Big[f(x_i) + 4 f(x_{i+1/2}) + f(x_{i+1})\Big]
$$

which is exactly **Simpson's 1/3 rule** on $[x_i, x_{i+1}]$ (with half-width $h/2$). RK4 is the natural ODE generalization of Simpson's rule, which explains the weights $1, 2, 2, 1$ (the midpoint weight 4 is split between the two midpoint stages). See `theory_simpson.md` in the integration theory.

## Butcher Barriers

For explicit RK methods, the minimum number of stages $s$ needed to reach order $p$ is:

| Order $p$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| Min. stages $s$ | 1 | 2 | 3 | 4 | 6 | 7 | 9 | 11 |

Up to order 4, each extra stage buys one order. Beyond that, extra stages are needed per order. This is why RK4 is the sweet spot of fixed-step explicit RK methods.

---

# RK4 Method

$$
\boxed{
\begin{aligned}
w_0 &= \alpha \\
k_1 &= f(x_i,\; w_i) \\
k_2 &= f\!\left(x_i + \tfrac{h}{2},\; w_i + \tfrac{h}{2}\,k_1\right) \\
k_3 &= f\!\left(x_i + \tfrac{h}{2},\; w_i + \tfrac{h}{2}\,k_2\right) \\
k_4 &= f\!\left(x_i + h,\; w_i + h\,k_3\right) \\
w_{i+1} &= w_i + \frac{h}{6}\,\big(k_1 + 2k_2 + 2k_3 + k_4\big)
\end{aligned}
}
$$

Interpretation:

- $k_1$: slope at the start (Euler).
- $k_2$: slope at the midpoint, reached with $k_1$.
- $k_3$: slope at the midpoint again, reached with the better estimate $k_2$.
- $k_4$: slope at the end, reached with $k_3$.

The update is a weighted average of these four slopes.

## Algorithm

**Input:** $f$, endpoints $a, b$, initial value $\alpha$, number of steps $N$ (or step size $h$).
**Output:** approximations $w_i \approx y(x_i)$ for $i = 0, \dots, N$.

```
h  ← (b − a) / N
w₀ ← α
for i = 0, 1, ..., N−1:
    x_i ← a + i · h
    k1 ← f(x_i,       w_i)
    k2 ← f(x_i + h/2, w_i + (h/2) · k1)
    k3 ← f(x_i + h/2, w_i + (h/2) · k2)
    k4 ← f(x_i + h,   w_i + h · k3)
    w_{i+1} ← w_i + (h/6) · (k1 + 2·k2 + 2·k3 + k4)
return (x, w)
```

**Cost:** four evaluations of $f$ per step. For systems the same algorithm applies with vectors (see `theory_rk4-system.md`).

---

# Error Analysis

## Local and Global Error

When $y \in C^5[a,b]$, the order conditions guarantee

$$
\tau_{i+1}(h) = O(h^4), \qquad |y(x_i) - w_i| = O(h^4)
$$

RK4 has **order 4**: halving $h$ divides the error by about 16.

**Equal-work comparison.** With a fixed budget of $4N$ function evaluations, RK4 with step $h$ competes against Heun/midpoint with step $h/2$ and Euler with step $h/4$. For smooth problems RK4 wins by orders of magnitude (see the example).

The order-4 estimate needs smoothness. If $f$ has discontinuities or kinks inside the interval, the observed order drops; the mesh should then be aligned with the discontinuities.

## Error Estimation by Step Doubling

Since the exact solution is generally unknown, the error can be estimated with **Richardson extrapolation**. Compute $w^{(h)}$ with one step of size $h$ and $w^{(h/2)}$ with two steps of size $h/2$. For an order-$p$ method:

$$
y(x_{i+1}) - w^{(h/2)} \approx \frac{w^{(h/2)} - w^{(h)}}{2^p - 1} = \frac{w^{(h/2)} - w^{(h)}}{15}
$$

This is the simplest basis for **adaptive step size control**. Embedded pairs such as RKF45 do the same job more cheaply by sharing stages between an order-4 and an order-5 formula.

---

# Stability

On $y' = \lambda y$, $z = h\lambda$:

$$
R(z) = 1 + z + \frac{z^2}{2} + \frac{z^3}{6} + \frac{z^4}{24}
$$

the degree-4 Taylor polynomial of $e^{z}$.

- Real stability interval: $-2.785 \lesssim h\lambda \le 0$, slightly larger than Euler/RK2.
- Unlike Euler and RK2, the region **includes a segment of the imaginary axis**, $|h\omega| \le 2\sqrt{2} \approx 2.83$. RK4 can integrate undamped oscillations with only a slight amplitude loss.
- The region is still bounded: RK4 is not suitable for stiff problems.

---

# Properties

- Explicit, one-step, self-starting.
- Four function evaluations per step.
- Local truncation error $O(h^4)$; global error $O(h^4)$ — **order 4**.
- Reduces to Simpson's 1/3 rule when $f = f(x)$.
- Reaches the maximum order possible for 4 stages.
- Real stability interval $\approx [-2.785, 0]$; includes part of the imaginary axis.
- Requires $y \in C^5[a,b]$ for the error estimate.

---

# Example

Approximate the solution of

$$
y' = y - x^2 + 1, \qquad 0 \le x \le 2, \qquad y(0) = 0.5
$$

with $N = 10$, $h = 0.2$. Exact solution: $y(x) = (x+1)^2 - 0.5\,e^{x}$.

**First step.**

$$
\begin{aligned}
k_1 &= f(0,\; 0.5) = 1.5 \\
k_2 &= f(0.1,\; 0.5 + 0.1 \cdot 1.5) = f(0.1,\; 0.65) = 1.64 \\
k_3 &= f(0.1,\; 0.5 + 0.1 \cdot 1.64) = f(0.1,\; 0.664) = 1.654 \\
k_4 &= f(0.2,\; 0.5 + 0.2 \cdot 1.654) = f(0.2,\; 0.8308) = 1.7908
\end{aligned}
$$

$$
w_1 = 0.5 + \frac{0.2}{6}\,\big(1.5 + 2(1.64) + 2(1.654) + 1.7908\big) = 0.8292933
$$

**Full table.**

| $i$ | $x_i$ | $k_1$ | $k_2$ | $k_3$ | $k_4$ | $w_i$ | $\lvert y(x_i) - w_i \rvert$ |
|---|---|---|---|---|---|---|---|
| 0 | 0.0 | — | — | — | — | 0.5000000 | 0 |
| 1 | 0.2 | 1.5000000 | 1.6400000 | 1.6540000 | 1.7908000 | 0.8292933 | $5.288 \times 10^{-6}$ |
| 2 | 0.4 | 1.7892933 | 1.9182227 | 1.9311156 | 2.0555165 | 1.2140762 | $1.144 \times 10^{-5}$ |
| 3 | 0.6 | 2.0540762 | 2.1694838 | 2.1810246 | 2.2902811 | 1.6489220 | $1.858 \times 10^{-5}$ |
| 4 | 0.8 | 2.2889220 | 2.3878142 | 2.3977034 | 2.4884627 | 2.1272027 | $2.685 \times 10^{-5}$ |
| 5 | 1.0 | 2.4872027 | 2.5659230 | 2.5737950 | 2.6419617 | 2.6408227 | $3.639 \times 10^{-5}$ |
| 6 | 1.2 | 2.6408227 | 2.6949050 | 2.7003132 | 2.7408853 | 3.1798942 | $4.737 \times 10^{-5}$ |
| 7 | 1.4 | 2.7398942 | 2.7638836 | 2.7662825 | 2.7731507 | 3.7323401 | $5.994 \times 10^{-5}$ |
| 8 | 1.6 | 2.7723401 | 2.7595741 | 2.7582975 | 2.7239996 | 4.2834095 | $7.429 \times 10^{-5}$ |
| 9 | 1.8 | 2.7234095 | 2.6657504 | 2.6599845 | 2.5754064 | 4.8150857 | $9.057 \times 10^{-5}$ |
| 10 | 2.0 | 2.5750857 | 2.4625943 | 2.4513451 | 2.3053547 | 5.3053630 | $1.089 \times 10^{-4}$ |

**Equal-work comparison** at $x = 2$ (40 evaluations of $f$ each):

| Method | $h$ | $\lvert e(2) \rvert$ |
|---|---|---|
| Euler | 0.05 | $1.274657 \times 10^{-1}$ |
| Heun | 0.1 | $1.890478 \times 10^{-2}$ |
| Midpoint | 0.1 | $3.747074 \times 10^{-3}$ |
| RK4 | 0.2 | $1.089498 \times 10^{-4}$ |

**Convergence.** Global error at $x = 2$:

| $h$ | $\lvert e \rvert$ | ratio |
|---|---|---|
| 0.2 | $1.089498 \times 10^{-4}$ | — |
| 0.1 | $6.990307 \times 10^{-6}$ | 15.59 |
| 0.05 | $4.421339 \times 10^{-7}$ | 15.81 |
| 0.025 | $2.778991 \times 10^{-8}$ | 15.91 |

Ratios ≈ 16 = $2^4$: order 4.

---

# Why RK4 Matters

- It is the **workhorse** of fixed-step IVP solvers: high accuracy, no derivatives of $f$, simple code.
- It is the conceptual bridge between quadrature (Simpson) and ODE solvers.
- In BRUNO it is reused in three places, so it should be a **pure, importable function** in `ode_ivp`:
  - as a method on its own (`rk4`),
  - as the **starter** for Adams–Bashforth and Adams–Moulton,
  - as the inner IVP solver of the **shooting** method in the future `ode_bvp`.
- Its step-doubling estimate is the first step towards adaptive methods (RKF45).

---

# Sources

- Burden, R. L., Faires, J. D., & Burden, A. M. (2016). *Numerical Analysis* (10th ed.). Cengage Learning. Section 5.4: Runge–Kutta Methods.
- Chapra, S. C., & Canale, R. P. (2015). *Numerical Methods for Engineers* (7th ed.). McGraw-Hill. Section 25.3.
- Hairer, E., Nørsett, S. P., & Wanner, G. (1993). *Solving Ordinary Differential Equations I: Nonstiff Problems* (2nd ed.). Springer. Sections II.1–II.3.
- Butcher, J. C. (2016). *Numerical Methods for Ordinary Differential Equations* (3rd ed.). Wiley.
- Kutta, W. (1901). Beitrag zur näherungsweisen Integration totaler Differentialgleichungen. *Zeitschrift für Mathematik und Physik*, 46, 435–453.