
---

# Adaptive Simpson Quadrature — Theory

---

# **Table of Contents**

- [Overview](#overview)
- [Mathematical Foundations](#mathematical-foundations)
  - [Composite Simpson Recap](#simpson-rule-recap)
  - [Error Estimation by Interval Halving](#error-estimation-by-interval-halving)
  - [Richardson Extrapolation](#richardson-extrapolation)
- [Adaptive Algorithm](#adaptive-algorithm)
  - [Stopping Criterion](#stopping-criterion)
  - [Tolerance Distribution](#tolerance-distribution)
  - [Recursive Formulation](#recursive-formulation)
  - [Priority-Queue Formulation](#priority-queue-formulation)
- [Cost and Function Evaluations](#cost-and-function-evaluations)
- [Properties](#properties)
- [Example](#worked-example)
- [Why Adaptive Simpson Is Powerful](#why-adaptive-simpson-is-powerful)
- [Limitations and Pitfalls](#limitations-and-pitfalls)
- [Sources](#sources)

---

# Overview

Adaptive Simpson quadrature is a **refinement strategy**, not a new quadrature rule. It uses the ordinary Simpson rule but decides **where** to place nodes automatically: it subdivides only the regions of $[a,b]$ where the integrand is hard to approximate, and leaves smooth regions coarse.

A fixed-step rule (Simpson, trapezoid) uses the same step $h$ everywhere. If the integrand has a sharp peak, an oscillation, or a steep derivative in a small part of the domain, $h$ must be tiny **everywhere** to resolve it. Adaptive Simpson instead:

1. estimates the error on each subinterval,
2. accepts the subinterval if the estimate is below its tolerance,
3. bisects it otherwise and repeats.

The goal is to reach a user-specified tolerance $\varepsilon$ with as few evaluations of $f$ as possible:

$$
\left| I - I_{\text{approx}} \right| \le \varepsilon,
\qquad
I = \int_a^b f(x)\,dx.
$$

---

# Mathematical Foundations

## Simpson Rule Recap

On a single interval $[a,b]$, with midpoint $c = \tfrac{a+b}{2}$ and $h = \tfrac{b-a}{2}$, the Simpson rule is:

$$
S(a,b) = \frac{b-a}{6}\Big[f(a) + 4f(c) + f(b)\Big].
$$

Its error, for $f \in C^4[a,b]$, is:

$$
I - S(a,b) = -\frac{h^5}{90}\, f^{(4)}(\xi), \qquad \xi \in (a,b).
$$

The error is $O(h^5)$ per interval, and $S$ is **exact for polynomials up to degree 3**.

---

## Error Estimation by Interval Halving

The true error depends on $f^{(4)}$, which we do not know. The adaptive method avoids computing it by using **two** Simpson approximations:

- $S_1 = S(a,b)$, one Simpson panel on the whole interval (step $h$),
- $S_2 = S(a,c) + S(c,b)$, two Simpson panels on the halves (step $h/2$).

Assuming $f^{(4)}$ is approximately constant on $[a,b]$, halving the step divides the error by $2^5 = 32$ per panel, and since two panels are added, the total error shrinks by $2^4 = 16$:

$$
I - S_1 \approx 16\,(I - S_2).
$$

Writing $I - S_1 = (I - S_2) + (S_2 - S_1)$ and substituting:

$$
15\,(I - S_2) \approx S_2 - S_1,
$$

which gives the error estimate:

$$
\boxed{
I - S_2 \approx \frac{S_2 - S_1}{15}
}
$$

The difference between the two approximations, divided by 15, measures how wrong $S_2$ is. No derivative is ever evaluated.

---

## Richardson Extrapolation

Since the error of $S_2$ is known to leading order, it can be **subtracted**:

$$
\boxed{
I \approx S_2 + \frac{S_2 - S_1}{15}
}
$$

This corrected value is one order higher: it is equivalent to **Boole's rule** (Newton–Cotes with 5 points), exact for polynomials up to degree 5. It costs nothing extra, so adaptive implementations return it instead of $S_2$.

---

# Adaptive Algorithm

## Stopping Criterion

A subinterval $[a,b]$ is **accepted** when the estimated error is below its tolerance $\varepsilon_{[a,b]}$:

$$
\frac{|S_2 - S_1|}{15} \le \varepsilon_{[a,b]}
\quad\Longleftrightarrow\quad
|S_2 - S_1| \le 15\,\varepsilon_{[a,b]}.
$$

Otherwise it is bisected into $[a,c]$ and $[c,b]$, and each half is processed with the same test.

---

## Tolerance Distribution

The errors of the accepted subintervals add up, so the local tolerances must add up to the global one. The classical rule splits the tolerance in half at every bisection:

$$
\varepsilon_{[a,c]} = \varepsilon_{[c,b]} = \frac{\varepsilon_{[a,b]}}{2}.
$$

An interval at depth $k$ therefore has tolerance $\varepsilon / 2^k$, and the sum over all accepted intervals never exceeds $\varepsilon$.

A common safeguard: the error estimate is a heuristic, so many implementations use a slightly stricter test (for example $|S_2-S_1| \le 10\,\varepsilon_{[a,b]}$) to compensate for the case where $f^{(4)}$ is **not** constant.

---

## Recursive Formulation

```text
adaptive_simpson(f, a, b, eps, S1, fa, fc, fb, depth):
    c   = (a + b) / 2
    d   = (a + c) / 2
    e   = (c + b) / 2
    fd  = f(d)
    fe  = f(e)

    S_left  = (c - a) / 6 * (fa + 4*fd + fc)
    S_right = (b - c) / 6 * (fc + 4*fe + fb)
    S2      = S_left + S_right

    if depth <= 0 or |S2 - S1| <= 15 * eps:
        return S2 + (S2 - S1) / 15

    return adaptive_simpson(f, a, c, eps/2, S_left,  fa, fd, fc, depth-1)
         + adaptive_simpson(f, c, b, eps/2, S_right, fc, fe, fb, depth-1)
```

Note that the values $f(a)$, $f(c)$, $f(b)$ are **passed down** to the children, so each bisection only needs **two new evaluations** ($f(d)$ and $f(e)$).

---

## Priority-Queue Formulation

Recursion fixes the tolerance split in advance and can hit the recursion limit. An equivalent iterative alternative maintains a **max-heap** keyed by the error estimate:

```text
1. Compute S1, S2 and err = |S2 - S1| / 15 on [a, b]; push (err, a, b, S2).
2. While  sum(err over heap) > max(abs_tol, rel_tol * |I_current|)
          and iterations < max_iter:
       pop the interval with the largest err
       bisect it, compute S1, S2, err for both halves
       push both halves
3. Return  sum(S2 + (S2 - S1)/15)  over all intervals in the heap.
```

Advantages:

- always refines the **worst** interval first (globally optimal use of evaluations),
- no recursion limit, and a hard cap on iterations or evaluations is trivial to add,
- the final heap is the list of subintervals, which is useful to **visualize where the method refined**.

---

# Cost and Function Evaluations

- The first panel costs 5 evaluations ($a$, $d$, $c$, $e$, $b$).
- Every bisection adds **2** new evaluations per refined interval.
- If $m$ intervals end up accepted, the total is:

$$
N_{\text{eval}} = 5 + 2\,(m - 1) = 2m + 3.
$$

For a smooth integrand, $m$ is small. For an integrand with a localized difficulty, $m$ grows only near it, while a fixed-step rule would need $m$ proportional to the finest required resolution **over the whole domain**.

---

# Properties

### Advantages
- Automatic: the user specifies a tolerance, not a step size  
- Concentrates evaluations where the integrand is difficult  
- Reuses evaluations between levels (only 2 new per bisection)  
- Returns the subintervals, useful for diagnostics and plots  
- Richardson correction gives Boole-level accuracy for free  

### Disadvantages
- The error estimate is a heuristic (assumes $f^{(4)}$ nearly constant)  
- Can be **fooled** by functions that look smooth at the 5 sampled points  
- Evaluations are sequential and data-dependent, so it is harder to vectorize than a fixed grid  
- Natural for 1D; multidimensional versions need a different formulation  

---

# Example

Consider:

$$
f(x) = x^4, \quad [0,1].
$$

The exact integral is:

$$
I = \frac{1}{5} = 0.2.
$$

**Step 1: one Simpson panel on $[0,1]$** ($c = 0.5$):

$$
S_1 = \frac{1}{6}\left[f(0) + 4f(0.5) + f(1)\right]
    = \frac{1}{6}\left[0 + 0.25 + 1\right]
    = 0.208333.
$$

**Step 2: two panels on $[0,0.5]$ and $[0.5,1]$:**

| Panel | Nodes | Value |
|---|---|---|
| $[0, 0.5]$ | $f(0)=0,\; f(0.25)=0.00390625,\; f(0.5)=0.0625$ | $0.006510$ |
| $[0.5, 1]$ | $f(0.5)=0.0625,\; f(0.75)=0.31640625,\; f(1)=1$ | $0.194010$ |

$$
S_2 = 0.006510 + 0.194010 = 0.200521.
$$

**Step 3: error estimate:**

$$
\frac{S_2 - S_1}{15} = \frac{0.200521 - 0.208333}{15} = -0.000521.
$$

The actual error of $S_2$ is $0.200521 - 0.2 = 0.000521$, so the estimate matches in magnitude (the integrand's fourth derivative is constant, so the assumption holds exactly).

**Step 4: decision.**

| Tolerance $\varepsilon$ | Estimated error | Decision |
|---|---|---|
| $10^{-3}$ | $5.2\times10^{-4} \le 10^{-3}$ | Accept the interval |
| $10^{-4}$ | $5.2\times10^{-4} > 10^{-4}$ | Bisect, continue with $\varepsilon/2$ on each half |

**Step 5: Richardson correction.** Accepting the interval and applying the correction:

$$
I \approx S_2 + \frac{S_2 - S_1}{15} = 0.200521 - 0.000521 = 0.200000.
$$

The result is exact: $x^4$ has degree $\le 5$, and the corrected value is Boole's rule. The whole computation used only **5 evaluations**.

---

# Why Adaptive Simpson Is Powerful

### ✔ Localized difficulty  
Integrands such as $\sqrt{x}$ near $0$, $\sin(1/x)$ near $0$, or a narrow Gaussian peak need small steps in one region only. Adaptive Simpson places nodes there and nowhere else.

### ✔ Tolerance-driven  
The user controls the accuracy directly. A fixed-step rule requires guessing $n$ by trial and error (as in a batch sweep over $n$ to find the optimum).

### ✔ Built-in diagnostics  
The final subdivision shows where the integrand was hard to integrate, and the estimated total error is a by-product, with no reference value needed.

### ✔ Foundation for production integrators  
The same refine-the-worst-interval structure is used by general-purpose integrators. The difference is the local rule: **Gauss–Kronrod** pairs (such as G7–K15) are preferred in libraries because their nodes nest, giving a high-order error estimate with fewer evaluations.

---

# Limitations and Pitfalls

- **Maximum depth is mandatory.** A singularity or discontinuity can make the estimate never fall below the tolerance, so a recursion depth (or iteration) limit must stop the process and report that the tolerance was not reached.
- **Premature acceptance.** If $f$ has a feature narrower than the spacing of the 5 initial nodes (a thin spike), $S_1$ and $S_2$ can agree by coincidence and the spike is never detected. Starting from a few initial panels (for example splitting $[a,b]$ into 4 or 8 pieces) mitigates this.
- **Roundoff floor.** When $\varepsilon$ approaches machine precision ($\sim10^{-15}$ relative), $S_2 - S_1$ is dominated by roundoff and the test may never pass. Use $\max(\varepsilon_{\text{abs}}, \varepsilon_{\text{rel}} |I|)$ and a sensible lower bound.
- **Non-integrable singularities.** For $\int_0^1 \frac{1}{x}\,dx$ the integral diverges and the method will exhaust its depth limit; this should be reported as an error, not returned as a number.
- **Dimension.** The scheme is naturally 1D. In $d$ dimensions, the analogue is adaptive cubature (error-ranked subboxes, each bisected along one axis at a time), or nested 1D adaptive integration axis by axis.

---

# Sources

- Burden, R. L. and Faires, J. D. (2011) "Adaptive Quadrature Methods." §4.6 in Numerical Analysis, 9th ed. Boston, MA: Brooks/Cole.
- McKeeman, W. M. (1962) "Algorithm 145: Adaptive Numerical Integration by Simpson's Rule." Communications of the ACM 5(12), 604.
- Lyness, J. N. (1969) "Notes on the Adaptive Simpson Quadrature Routine." Journal of the ACM 16(3), 483–495.
- Press, W. H. et al. (1992) "Romberg Integration" and "Gaussian Quadratures and Orthogonal Polynomials." §4.3 and §4.5 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: Cambridge University Press.

---