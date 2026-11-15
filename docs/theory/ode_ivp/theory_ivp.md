
---

# Initial Value Problems — Theory

---

# **Table of Contents**

- [Overview](#overview)
- [Mathematical Foundations](#mathematical-foundations)
  - [Problem Statement](#problem-statement)
  - [Systems and Higher-Order Equations](#systems-and-higher-order-equations)
  - [Existence and Uniqueness](#existence-and-uniqueness)
  - [Well-Posedness](#well-posedness)
- [Discretization](#discretization)
  - [The Mesh](#the-mesh)
  - [One-Step vs Multistep Methods](#one-step-vs-multistep-methods)
- [Error Analysis](#error-analysis)
  - [Local Truncation Error](#local-truncation-error)
  - [Global Error](#global-error)
  - [Consistency, Order and Convergence](#consistency-order-and-convergence)
- [Stability](#stability)
  - [The Test Equation](#the-test-equation)
  - [Region of Absolute Stability](#region-of-absolute-stability)
  - [Stiffness](#stiffness)
- [Properties](#properties)
- [Example](#example)
- [Why IVP Theory Matters](#why-ivp-theory-matters)
- [Sources](#sources)

---

# Overview

An **initial value problem** (IVP) asks for a function $y(x)$ that satisfies a differential equation **and** takes a prescribed value at a single starting point $x_0$.

All the information needed to determine the solution is known at **one point**, so the solution can be built by *marching forward* from $x_0$, step by step. This is what distinguishes IVPs from **boundary value problems** (BVPs), where conditions are split between the two ends of the interval (see `theory_shooting.md` and `theory_finite-differences-bvp.md`).

Numerical IVP methods do not produce a formula for $y(x)$; they produce a table of approximations

$$
w_i \approx y(x_i), \qquad i = 0, 1, \dots, N
$$

on a discrete set of points called the **mesh**.

---

# Mathematical Foundations

## Problem Statement

Let $f : [a,b] \times \mathbb{R} \to \mathbb{R}$. The scalar first-order IVP is

$$
y'(x) = f\big(x,\, y(x)\big), \qquad a \le x \le b, \qquad y(a) = \alpha
$$

where:

- $x$ is the independent variable,
- $y$ is the unknown function,
- $f(x,y)$ gives the **slope** of the solution at each point $(x,y)$,
- $\alpha$ is the **initial condition**.

Geometrically, $f$ defines a **direction field**: at every point of the plane it prescribes the slope a solution passing through that point must have. The initial condition selects *one* curve from the infinite family of solutions.

## Systems and Higher-Order Equations

The same framework covers **systems** of $m$ equations by letting $\mathbf{y}$ be a vector:

$$
\mathbf{y}'(x) = \mathbf{f}\big(x, \mathbf{y}(x)\big), \qquad \mathbf{y}(a) = \boldsymbol{\alpha}, \qquad \mathbf{y}, \mathbf{f}, \boldsymbol{\alpha} \in \mathbb{R}^m
$$

Every one-step method written for the scalar case applies to systems **component-wise, with the same formula**, replacing scalars by vectors.

A higher-order equation can always be reduced to a first-order system. For

$$
y^{(m)} = f\big(x, y, y', \dots, y^{(m-1)}\big), \qquad y(a)=\alpha_1,\; y'(a)=\alpha_2,\; \dots,\; y^{(m-1)}(a)=\alpha_m
$$

define $u_1 = y,\; u_2 = y',\; \dots,\; u_m = y^{(m-1)}$. Then

$$
\begin{aligned}
u_1' &= u_2 \\
u_2' &= u_3 \\
&\;\;\vdots \\
u_m' &= f(x, u_1, u_2, \dots, u_m)
\end{aligned}
$$

with $u_j(a) = \alpha_j$. This is exactly the form consumed by the `rk4_system` method.

> **Note.** A second-order equation with **both** conditions at $x = a$ (value and slope) is an IVP. With one condition at $x=a$ and another at $x=b$ it is a BVP, even though the differential equation is identical.

## Existence and Uniqueness

**Definition (Lipschitz condition).** $f(x,y)$ satisfies a Lipschitz condition in $y$ on

$$
D = \{(x,y) : a \le x \le b,\; -\infty < y < \infty\}
$$

if there exists a constant $L > 0$ such that

$$
|f(x, y_1) - f(x, y_2)| \le L\,|y_1 - y_2| \qquad \forall\, (x,y_1),(x,y_2) \in D
$$

$L$ is called the **Lipschitz constant**.

A practical sufficient condition: if $\partial f / \partial y$ exists and is bounded on $D$, then

$$
L = \sup_{(x,y)\in D} \left| \frac{\partial f}{\partial y}(x,y) \right|
$$

works as a Lipschitz constant (by the Mean Value Theorem).

**Theorem (Picard–Lindelöf).** If $f$ is continuous on $D$ and Lipschitz in $y$, then the IVP

$$
y' = f(x,y), \qquad y(a) = \alpha
$$

has a **unique** solution $y(x)$ for $a \le x \le b$.

Without the Lipschitz condition uniqueness may fail. Classic example:

$$
y' = \sqrt{|y|}, \qquad y(0) = 0
$$

has both $y \equiv 0$ and $y = x^2/4$ (for $x \ge 0$) as solutions; $\partial f/\partial y$ is unbounded near $y=0$.

## Well-Posedness

An IVP is **well-posed** if:

1. a unique solution exists, and
2. small perturbations in the data produce small changes in the solution.

Formally, for the perturbed problem

$$
z' = f(x,z) + \delta(x), \qquad z(a) = \alpha + \varepsilon_0
$$

with $|\varepsilon_0| < \varepsilon$ and $|\delta(x)| < \varepsilon$, there exists $k(\varepsilon) > 0$ such that

$$
|z(x) - y(x)| < k(\varepsilon)\,\varepsilon \qquad \forall\, x \in [a,b]
$$

**Theorem.** If $f$ is continuous and Lipschitz in $y$ on $D$, the IVP is well-posed.

This matters for numerical work because every method introduces perturbations (rounding, truncation). Well-posedness guarantees those perturbations do not, by themselves, destroy the solution.

---

# Discretization

## The Mesh

The interval $[a,b]$ is divided into $N$ subintervals of equal length

$$
h = \frac{b - a}{N}, \qquad x_i = a + i\,h, \qquad i = 0, 1, \dots, N
$$

where $h$ is the **step size**. The notation used across BRUNO's ODE theory documents is:

| Symbol | Meaning |
|---|---|
| $y(x_i)$ | exact solution at the mesh point |
| $w_i$ | numerical approximation of $y(x_i)$ |
| $f_i$ | $f(x_i, w_i)$ |
| $h$ | step size |
| $N$ | number of steps |

The user may specify either $h$ or $N$; they are related by $N = (b-a)/h$, which must be an integer for the mesh to land exactly on $b$.

## One-Step vs Multistep Methods

**One-step methods** compute $w_{i+1}$ using only $(x_i, w_i)$:

$$
w_{i+1} = w_i + h\,\phi(x_i, w_i; h)
$$

where $\phi$ is the **increment function**, a weighted average of slopes. They are *self-starting*: $w_0 = \alpha$ is enough.

**Multistep methods** use several previous values $w_i, w_{i-1}, \dots, w_{i-k+1}$:

$$
w_{i+1} = \sum_{j=0}^{k-1} a_j\, w_{i-j} + h \sum_{j=-1}^{k-1} b_j\, f_{i-j}
$$

If $b_{-1} = 0$ the method is **explicit** (Adams–Bashforth); otherwise it is **implicit** (Adams–Moulton). Multistep methods are **not self-starting**: the first $k-1$ values must come from a one-step method.

Classification of the IVP methods in BRUNO:

| Method | Type | Order | Self-starting |
|---|---|---|---|
| Euler | one-step, explicit | 1 | yes |
| Heun | one-step, explicit (RK2 family) | 2 | yes |
| RK2 (midpoint) | one-step, explicit (RK2 family) | 2 | yes |
| RK4 | one-step, explicit | 4 | yes |
| RK4 system | one-step, explicit, vector | 4 | yes |
| Adams–Bashforth 2 | 2-step, explicit | 2 | no (needs $w_1$) |
| Adams–Bashforth 3 | 3-step, explicit | 3 | no (needs $w_1, w_2$) |
| Adams–Moulton 2 | 2-step, implicit (PECE in practice) | 3* | no (needs $w_1$, plus $w_2$ with AB3 predictor) |

\* Naming varies between sources: "AM2" may denote the 1-step trapezoidal rule (order 2) or the 2-step Adams–Moulton formula (order 3). BRUNO's theory documents use the 2-step, order-3 convention (see `theory_adams-moulton-2.md`); the executor must match it.

---

# Error Analysis

## Local Truncation Error

The **local truncation error** measures how badly the *exact* solution fails to satisfy the method's formula for one step. For a one-step method:

$$
\tau_{i+1}(h) = \frac{y(x_{i+1}) - y(x_i)}{h} - \phi\big(x_i, y(x_i); h\big)
$$

It answers: *if I start a step at the exact value, how much error does that single step introduce (per unit $h$)?*

The error committed in one step is $h\,\tau_{i+1}(h)$. Many textbooks call this quantity the local error directly; be careful which convention a source uses.

## Global Error

The **global error** is what the user actually sees:

$$
e_i = y(x_i) - w_i
$$

It accumulates the local errors of all previous steps, each one propagated (and possibly amplified) by the differential equation.

Since reaching $b$ takes $N = (b-a)/h$ steps, a method with local error $O(h^{p+1})$ per step produces, under the Lipschitz hypothesis,

$$
|e_i| = O(h^p)
$$

i.e. **one power of $h$ is lost** when local errors accumulate.

## Consistency, Order and Convergence

- **Consistency.** A method is consistent if $\tau_i(h) \to 0$ as $h \to 0$. For one-step methods this is equivalent to $\phi(x, y; 0) = f(x,y)$.
- **Order.** A method has order $p$ if $\tau_i(h) = O(h^p)$. Then the global error is also $O(h^p)$.
- **Convergence.** A method is convergent if

$$
\lim_{h \to 0} \max_{0 \le i \le N} |y(x_i) - w_i| = 0
$$

**Theorem.** For a one-step method whose increment function $\phi$ is Lipschitz in $w$, consistency implies convergence.

**Dahlquist Equivalence Theorem.** For linear multistep methods:

$$
\text{consistency} + \text{zero-stability} \iff \text{convergence}
$$

**Practical test of order.** If a method has order $p$, halving $h$ divides the global error by approximately $2^p$:

$$
\frac{|e(h)|}{|e(h/2)|} \approx 2^p
$$

This is the basis for regression tests that verify an implementation reaches its theoretical order.

**Rounding error.** In floating-point arithmetic the total error behaves like

$$
|e| \approx C_1\, h^p + \frac{C_2\, \varepsilon_{\text{mach}}}{h}
$$

so decreasing $h$ indefinitely eventually *increases* the error. There is an optimal $h$ below which refinement is counter-productive.

---

# Stability

Convergence is an asymptotic statement ($h \to 0$). In practice $h$ is finite, and a convergent method can still produce garbage if $h$ is too large. Stability analysis quantifies that.

## The Test Equation

Stability is studied on the **Dahlquist test equation**

$$
y' = \lambda y, \qquad y(0) = 1, \qquad \lambda \in \mathbb{C},\; \operatorname{Re}(\lambda) < 0
$$

whose exact solution $y = e^{\lambda x}$ decays to zero. A reasonable method should reproduce that decay.

Applying any explicit one-step method to the test equation yields

$$
w_{i+1} = R(z)\, w_i, \qquad z = h\lambda
$$

where $R(z)$ is the **stability function** of the method. Therefore $w_i = R(z)^i$.

## Region of Absolute Stability

The **region of absolute stability** is

$$
\mathcal{S} = \{ z \in \mathbb{C} : |R(z)| \le 1 \}
$$

If $z = h\lambda \notin \mathcal{S}$, the numerical solution grows without bound even though the exact solution decays.

| Method | $R(z)$ | Real stability interval |
|---|---|---|
| Euler | $1 + z$ | $-2 \le z \le 0$ |
| Heun / RK2 | $1 + z + \tfrac{z^2}{2}$ | $-2 \le z \le 0$ |
| RK4 | $1 + z + \tfrac{z^2}{2} + \tfrac{z^3}{6} + \tfrac{z^4}{24}$ | $\approx -2.785 \le z \le 0$ |
| AB2 | — (multistep) | $-1 \le z \le 0$ |
| AB3 | — (multistep) | $\approx -0.545 \le z \le 0$ |
| AB3–AM2 PECE | — (multistep) | $\approx -1.73 \le z \le 0$ |

Multistep methods have no single stability function; their regions come from the roots of $\rho(\zeta) - z\,\sigma(\zeta)$ (see `theory_adams-bashford.md`).

For a nonlinear problem, $\lambda$ plays the role of $\partial f/\partial y$ along the solution (for systems, the eigenvalues of the Jacobian).

## Stiffness

A problem is **stiff** when it contains components that decay very fast ($|\lambda|$ large, $\operatorname{Re}(\lambda) \ll 0$) alongside components that evolve slowly.

The fast components become negligible almost immediately, yet they still force $|h\lambda|$ to stay inside $\mathcal{S}$. Every explicit method has a **bounded** stability region, so a stiff problem forces a tiny $h$ for the *entire* interval, even where the solution is smooth.

Stiff problems call for implicit methods (backward Euler, BDF, implicit Runge–Kutta), whose stability regions contain the whole left half-plane. None of BRUNO's current IVP methods are designed for stiff problems.

---

# Properties

- An IVP is determined by $f$, the interval $[a,b]$ and **one** initial condition per unknown.
- Lipschitz continuity of $f$ in $y$ guarantees existence, uniqueness and well-posedness.
- Higher-order equations are always reduced to first-order systems.
- Global error is one order lower in $h$ than the per-step local error.
- Order $p$ means halving $h$ divides the error by about $2^p$.
- Explicit methods have bounded stability regions: $h$ must satisfy $h|\lambda| \lesssim 2$ (about $2.785$ for RK4, and much less for Adams–Bashforth).
- Multistep methods need a one-step starter for the first values.

---

# Example

Consider the IVP used throughout BRUNO's ODE theory documents:

$$
y' = y - x^2 + 1, \qquad 0 \le x \le 2, \qquad y(0) = 0.5
$$

**Exact solution.**

$$
y(x) = (x+1)^2 - 0.5\,e^{x}, \qquad y(2) = 5.3054720
$$

**Lipschitz check.**

$$
\frac{\partial f}{\partial y} = 1 \quad\Rightarrow\quad L = 1
$$

$f$ is continuous and Lipschitz on $D$, so the IVP has a unique solution and is well-posed.

**Stability check.** Since $\lambda \approx \partial f/\partial y = 1 > 0$, the exact solution is *not* decaying in the sense of the test equation; there is no stability restriction on $h$ here, only accuracy.

**Empirical order.** Global error at $x = 2$ for Euler (order 1) and Heun (order 2):

| $h$ | $N$ | Euler $\lvert e\rvert$ | ratio | Heun $\lvert e\rvert$ | ratio |
|---|---|---|---|---|---|
| 0.2 | 10 | $4.396874 \times 10^{-1}$ | — | $7.241732 \times 10^{-2}$ | — |
| 0.1 | 20 | $2.419719 \times 10^{-1}$ | 1.82 | $1.890478 \times 10^{-2}$ | 3.83 |
| 0.05 | 40 | $1.274657 \times 10^{-1}$ | 1.90 | $4.819865 \times 10^{-3}$ | 3.92 |
| 0.025 | 80 | $6.549505 \times 10^{-2}$ | 1.95 | $1.216136 \times 10^{-3}$ | 3.96 |

The ratios approach $2^1 = 2$ for Euler and $2^2 = 4$ for Heun, confirming their theoretical orders.

**Instability illustration.** For $y' = -20y$, $y(0)=1$ with Euler, $R(z) = 1 - 20h$:

| $h$ | $z = -20h$ | $R(z)$ | $w_1,\dots,w_5$ |
|---|---|---|---|
| 0.15 | $-3$ | $-2$ | $-2,\; 4,\; -8,\; 16,\; -32$ (diverges) |
| 0.10 | $-2$ | $-1$ | $-1,\; 1,\; -1,\; 1,\; -1$ (oscillates) |
| 0.05 | $-1$ | $0$ | $0,\; 0,\; 0,\; 0,\; 0$ (stable) |

The exact solution $e^{-20x}$ decays monotonically; only the last step size is inside the stability region.

---

# Why IVP Theory Matters

These results translate directly into the layers of BRUNO's `ode_ivp` module:

- **Constructor.** Requires $a < b$, a parseable $f(x,y)$ (or $f(x, y_0, \dots, y_{m-1})$ for systems), and one initial value per equation.
- **Validator.**
  - $h > 0$ and $N \ge 1$; if both are given, they must satisfy $N = (b-a)/h$.
  - For systems, the number of equations must equal the number of initial values.
  - Multistep methods require $N \ge k$ (enough steps after the starter values).
  - Optional warning when $h\,|\partial f/\partial y|$ falls outside the method's real stability interval.
- **Executor.** Methods are pure functions of $(f, a, \alpha, h, N)$, so the same RK4 can be reused by multistep starters and, later, by the shooting method in `ode_bvp`.
- **Tests.** The convergence ratio $|e(h)|/|e(h/2)| \approx 2^p$ is a robust, implementation-independent regression test.

---

# Sources

- Burden, R. L., Faires, J. D., & Burden, A. M. (2016). *Numerical Analysis* (10th ed.). Cengage Learning. Chapter 5: Initial-Value Problems for Ordinary Differential Equations.
- Chapra, S. C., & Canale, R. P. (2015). *Numerical Methods for Engineers* (7th ed.). McGraw-Hill. Part Seven: Ordinary Differential Equations.
- Hairer, E., Nørsett, S. P., & Wanner, G. (1993). *Solving Ordinary Differential Equations I: Nonstiff Problems* (2nd ed.). Springer.
- Süli, E., & Mayers, D. (2003). *An Introduction to Numerical Analysis*. Cambridge University Press. Chapter 12.
- Atkinson, K. E. (1989). *An Introduction to Numerical Analysis* (2nd ed.). Wiley. Chapter 6.