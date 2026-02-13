# Mathematical Background

This document explains the mathematical theory behind Gelfgren's algorithms.

## Table of Contents

- [Bernstein Polynomials](#bernstein-polynomials)
- [Rational Functions](#rational-functions)
- [Padé Approximants](#padé-approximants)
- [Lagrange-Hermite Interpolation](#lagrange-hermite-interpolation)
- [Piecewise Rational Interpolation](#piecewise-rational-interpolation)
- [Boundary Value Problems](#boundary-value-problems)
- [Numerical Properties](#numerical-properties)

## Bernstein Polynomials

### Definition

A Bernstein polynomial of degree $n$ on the interval $[a,b]$ is defined as:

$$
B_n(x) = \sum_{i=0}^{n} b_i \cdot B_{i,n}(t)
$$

where $t = (x-a)/(b-a)$ normalizes $x$ to $[0,1]$, and

$$
B_{i,n}(t) = \binom{n}{i} t^i (1-t)^{n-i}
$$

are the Bernstein basis polynomials.

### Properties

1. **Partition of Unity**: $\sum_{i=0}^{n} B_{i,n}(t) = 1$ for all $t \in [0,1]$

2. **Endpoint Interpolation**:
   - $B_n(a) = b_0$
   - $B_n(b) = b_n$

3. **Convex Hull Property**: The polynomial lies within the convex hull of its control points $(i/(n), b_i)$

4. **Variation Diminishing**: The polynomial oscillates less than its control polygon

5. **Degree Elevation**: A degree $n$ polynomial can be exactly represented at degree $n+1$:

$$
b_i' = \frac{i}{n+1} b_{i-1} + \frac{n+1-i}{n+1} b_i
$$

### Why Bernstein Basis?

The Bernstein basis provides superior numerical stability compared to the power basis $(1, x, x^2, \ldots, x^n)$:

**Condition Number**: For degree $n$ polynomials on $[0,1]$:
- Power basis: $\kappa \approx 2^n$ (exponential growth)
- Bernstein basis: $\kappa \approx O(\sqrt{n})$ (mild growth)

This is proven in Farouki & Rajan (1987).

### Operations

**Derivative**:

$$
\frac{d}{dx} B_n(x) = \frac{n}{b-a} \sum_{i=0}^{n-1} (b_{i+1} - b_i) B_{i,n-1}(t)
$$

**Integral**:

$$
\int B_n(x) dx = (b-a) \sum_{i=0}^{n+1} \left(\frac{1}{n+1} \sum_{j=0}^{i-1} b_j \right) B_{i,n+1}(t) + C
$$

**Addition**: If $P(x) = \sum b_i B_{i,n}(t)$ and $Q(x) = \sum c_i B_{i,n}(t)$:

$$
P(x) + Q(x) = \sum (b_i + c_i) B_{i,n}(t)
$$

**Multiplication**: Requires degree elevation to common degree and convolution of coefficients.

## Rational Functions

### Definition

A rational function is the ratio of two polynomials:

$$
R(x) = \frac{P(x)}{Q(x)} = \frac{\sum_{i=0}^{n} p_i B_{i,n}(t)}{\sum_{j=0}^{m} q_j B_{j,m}(t)}
$$

where $P(x)$ and $Q(x)$ are Bernstein polynomials.

### Properties

1. **Poles**: $R(x)$ has poles where $Q(x) = 0$

2. **Asymptotic Behavior**: As $x \to \infty$:
   - If $n < m$: $R(x) \to 0$
   - If $n = m$: $R(x) \to p_n/q_m$
   - If $n > m$: $|R(x)| \to \infty$

3. **Flexibility**: Rational functions can approximate a wider class of functions than polynomials alone, including those with singularities

### Operations

**Evaluation**: Direct substitution with pole checking

**Derivative**: Using the quotient rule:

$$
R'(x) = \frac{P'(x)Q(x) - P(x)Q'(x)}{Q(x)^2}
$$

**Arithmetic**:

$$
\frac{P_1}{Q_1} + \frac{P_2}{Q_2} = \frac{P_1 Q_2 + P_2 Q_1}{Q_1 Q_2}
$$

$$
\frac{P_1}{Q_1} \cdot \frac{P_2}{Q_2} = \frac{P_1 P_2}{Q_1 Q_2}
$$

## Padé Approximants

### Definition

A Padé approximant $[n/m]$ to a function $f(x)$ is a rational function $R(x) = P_n(x)/Q_m(x)$ such that:

$$
f(x) - R(x) = O(x^{n+m+1})
$$

near $x = 0$ (or another expansion point).

### Construction from Power Series

Given power series coefficients $f(x) = \sum_{k=0}^{\infty} c_k x^k$, find $P_n(x) = \sum_{i=0}^{n} p_i x^i$ and $Q_m(x) = \sum_{j=0}^{m} q_j x^j$ with $q_0 = 1$ such that:

$$
f(x) Q_m(x) - P_n(x) = O(x^{n+m+1})
$$

This gives the linear system:

$$
\begin{bmatrix}
c_{n+1} & c_{n} & \cdots & c_{n-m+2} \\
c_{n+2} & c_{n+1} & \cdots & c_{n-m+3} \\
\vdots & \vdots & \ddots & \vdots \\
c_{n+m} & c_{n+m-1} & \cdots & c_{n+1}
\end{bmatrix}
\begin{bmatrix}
q_m \\
q_{m-1} \\
\vdots \\
q_1
\end{bmatrix}
=
-\begin{bmatrix}
c_{n+2} \\
c_{n+3} \\
\vdots \\
c_{n+m+1}
\end{bmatrix}
$$

Then $p_i = \sum_{j=0}^{i} c_{i-j} q_j$ for $i = 0, \ldots, n$.

### Construction from Derivative Data

Given derivative values $f^{(k)}(x_0)$ for $k = 0, \ldots, n+m$, first construct the Taylor series:

$$
f(x) \approx \sum_{k=0}^{n+m} \frac{f^{(k)}(x_0)}{k!} (x - x_0)^k
$$

Then apply the power series construction above.

### Symmetric Padé Approximants (Gelfgren's Method)

For a function on $[a, b]$, construct a symmetric Padé approximant matching derivative data at both endpoints $a$ and $b$. This provides better approximation across the entire interval.

## Lagrange-Hermite Interpolation

### Problem Statement

Given data points $(x_i, y_i)$ and derivative values $y_i^{(k)}$ for $i = 1, \ldots, N$ and $k = 0, \ldots, m_i - 1$, find a polynomial $P(x)$ such that:

$$
P^{(k)}(x_i) = y_i^{(k)} \quad \text{for all } i, k
$$

The polynomial has degree $n = \sum_{i=1}^{N} m_i - 1$.

### Traub's Formula

Traub (1964) derived an explicit formula using Bell polynomials. For the case of two points with derivatives up to order $m-1$:

$$
P(x) = \sum_{k=0}^{m-1} \left[ y_a^{(k)} S_k(x; a, b) + y_b^{(k)} S_k(x; b, a) \right]
$$

where $S_k(x; a, b)$ are the fundamental functions satisfying:

$$
S_k^{(j)}(a) = \delta_{jk}, \quad S_k^{(j)}(b) = 0 \quad \text{for } j = 0, \ldots, m-1
$$

### Bell Polynomials

The Bell polynomial $B_n(x_1, x_2, \ldots, x_n)$ appears in Faà di Bruno's formula for derivatives of composed functions. It's used in Traub's formula to construct the fundamental functions.

Recursive definition:

$$
B_0 = 1
$$

$$
B_{n+1}(x_1, \ldots, x_{n+1}) = \sum_{k=0}^{n} \binom{n}{k} x_{k+1} B_{n-k}(x_1, \ldots, x_{n-k})
$$

## Piecewise Rational Interpolation

### Gelfgren's Method

Given a mesh $a = x_0 < x_1 < \cdots < x_N = b$ and function values at mesh points, Gelfgren's method constructs a piecewise rational function $R(x)$ such that:

1. On each subinterval $[x_i, x_{i+1}]$, $R(x)$ is a rational function of type $[n/m]$

2. $R(x)$ matches function (and possibly derivative) values at mesh points

3. $R(x)$ has specified smoothness (continuity of derivatives) at interior mesh points

### Construction Algorithm

For each subinterval $[x_i, x_{i+1}]$:

1. **Select degrees** $n$ and $m$ (often $n = m$ for symmetric approximants)

2. **Construct Padé approximant** using derivative data at endpoints $x_i$ and $x_{i+1}$:
   - Use symmetric two-point Padé construction
   - Match function and derivative values

3. **Check smoothness**: Verify that left and right derivatives match at interior points

4. **Refine mesh** if needed: If error estimates are too large, subdivide intervals

### Error Estimation

On each subinterval, the error can be estimated using Hermite's formula:

$$
f(x) - R(x) = \frac{f^{(n+m+1)}(\xi)}{(n+m+1)!} \omega(x)
$$

where $\omega(x)$ is a product of linear factors involving the interpolation points, and $\xi \in [x_i, x_{i+1}]$.

### Mesh Strategies

**Uniform Mesh**: $x_i = a + i h$ where $h = (b-a)/N$
- Simple
- May require many points for functions with varying behavior

**Chebyshev Mesh**: $x_i = \frac{a+b}{2} + \frac{b-a}{2} \cos\left(\frac{(2i+1)\pi}{2(N+1)}\right)$
- Clustered near endpoints
- Near-optimal for polynomial approximation
- Minimizes Runge phenomenon

**Adaptive Mesh**: Start with coarse mesh, refine where error estimates are large
- Efficient
- Concentrates effort where needed

## Boundary Value Problems

### Problem Statement

Find $y(x)$ satisfying:

$$
\mathcal{L}[y] = f(x) \quad \text{on } [a, b]
$$

with boundary conditions:

$$
\begin{aligned}
\alpha_0 y(a) + \alpha_1 y'(a) &= \beta_a \\
\gamma_0 y(b) + \gamma_1 y'(b) &= \beta_b
\end{aligned}
$$

where $\mathcal{L}$ is a differential operator (e.g., $\mathcal{L}[y] = y'' + p(x)y' + q(x)y$).

### Collocation Method

Approximate $y(x)$ with a piecewise rational function $R(x)$ on a mesh. At collocation points $\{t_j\}$:

$$
\mathcal{L}[R](t_j) = f(t_j)
$$

Combined with boundary conditions, this gives a system of equations for the coefficients of $R(x)$.

### Advantages of Rational Approximation

- Can handle singularities in the solution
- Higher accuracy per degree of freedom than polynomials
- Natural for problems with boundary layers or rapid variations

## Numerical Properties

### Conditioning and Stability

**Bernstein Basis**: The condition number of the Bernstein basis matrix grows only as $O(\sqrt{n})$, making it ideal for high-degree polynomials.

**Power Basis**: Condition number grows exponentially, leading to severe numerical instability.

**Comparison**: For degree 20 on $[0,1]$:
- Bernstein: $\kappa \approx 10$
- Power: $\kappa \approx 10^6$

### Error Analysis

**Polynomial Approximation**: For smooth functions on $[a, b]$:

$$
\|f - P_n\| \leq C \frac{M_{n+1}}{(n+1)!} h^{n+1}
$$

where $M_{n+1} = \max |f^{(n+1)}(x)|$ and $h = b - a$.

**Rational Approximation**: Can achieve exponential convergence for analytic functions:

$$
\|f - R_{n,n}\| \leq C e^{-\alpha n}
$$

### Roundoff Error

Bernstein polynomials minimize the effect of coefficient perturbations:

$$
\|\delta B_n\| \leq \max_i |\delta b_i|
$$

This is the **optimal** property among all polynomial bases.

## Computational Complexity

| Operation | Bernstein | Power Basis |
|-----------|-----------|-------------|
| Evaluation | $O(n)$ | $O(n)$ |
| Degree elevation | $O(n)$ | $O(n^2)$ |
| Derivative | $O(n)$ | $O(n)$ |
| Multiplication | $O(n^2)$ | $O(n^2)$ |
| Root finding | $O(n^2)$ | $O(n^2)$ |

The Bernstein basis has the same or better complexity for all operations, with vastly superior numerical stability.

## References

1. **Gelfgren, J. (1975)**. "Piecewise Rational Interpolation". *BIT Numerical Mathematics*, 15, 382-393.
   - Original development of piecewise rational methods
   - Symmetric Padé approximants on subintervals

2. **Farouki, R. T., & Rajan, V. T. (1987)**. "Algorithms for Polynomials in Bernstein Form". *Computer Aided Geometric Design*, 5(1), 1-26.
   - Complete set of Bernstein polynomial algorithms
   - Proof of superior conditioning

3. **Traub, J. F. (1964)**. "On Lagrange-Hermite Interpolation". *SIAM Journal on Numerical Analysis*, 1(1), 1-15.
   - Bell polynomial formulation
   - Explicit formulas for Hermite interpolation

4. **Baker, G. A., & Graves-Morris, P. (1996)**. *Padé Approximants* (2nd ed.). Cambridge University Press.
   - Comprehensive treatment of Padé theory
   - Applications to physical problems

5. **Farin, G. (2002)**. *Curves and Surfaces for CAGD: A Practical Guide* (5th ed.). Morgan Kaufmann.
   - Geometric interpretation of Bernstein polynomials
   - Computer-aided design applications

6. **Powell, M. J. D. (1981)**. *Approximation Theory and Methods*. Cambridge University Press.
   - General approximation theory
   - Error estimates and convergence

## Notation Summary

| Symbol | Meaning |
|--------|---------|
| $B_{i,n}(t)$ | Bernstein basis polynomial of degree $n$, index $i$ |
| $b_i$ | Bernstein coefficient (control point) |
| $[n/m]$ | Padé approximant with numerator degree $n$, denominator degree $m$ |
| $\mathcal{L}$ | Differential operator |
| $\omega(x)$ | Node polynomial for error estimation |
| $\kappa$ | Condition number |
| $O(\cdot)$ | Big-O notation (asymptotic behavior) |
