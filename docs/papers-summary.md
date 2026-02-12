# Research Papers Summary

This document summarizes the three key papers that form the theoretical foundation of the Gelfgren library.

## 1. Gelfgren 1975: "Piecewise Rational Interpolation"

**Core Concept**: Constructing piecewise rational approximations by using Padé approximants on subintervals.

**Key Ideas**:
- **Padé Approximant [n,m]**: Rational function P_n(z)/Q_m(z) matching a power series to order m+n+1
  - P_n, Q_m are polynomials of degree n, m
  - Coefficients determined by linear system of equations
- **Symmetric Padé Approximant**: At two points z₀ and z₁ using Newton series
  - n+m+1 = 2p for symmetric case
  - Provides contact of highest order at both endpoints
- **Piecewise Rational Function S_{n,m}(τ)(z)**:
  - Partition interval into equal subintervals Δⱼ = [x_{j-1}, x_j]
  - Construct symmetric [n,m] Padé approximant on each subinterval
  - If P_{n,j} and Q_{m,j} have no common factors → p-1 continuous derivatives at nodes
- **Error Estimation**: Uses Hermite's interpolation formula and capacity theory
  - Requires n ≥ k·m for convergence
  - Better approximation near boundaries than in middle

**Mathematical Components**:
1. Newton series: f(z) ~ Σ a_k w_k(z) where w_k(z) = (z-z₀)^{k-⌊k/2⌋}(z-z₁)^{⌊k/2⌋}
2. Linear equations for coefficients: Σ a_{j-k}β_k = {a_j if 0≤j≤n; 0 if n+1≤j≤n+m}
3. Piecewise construction on partition τ: a = x₀ < x₁ < ... < x_N = b

## 2. Traub 1964: "On Lagrange-Hermite Interpolation"

**Core Concept**: Formula for interpolating polynomial matching function and derivatives at multiple points.

**Key Ideas**:
- **Lagrange-Hermite Problem**: Given p(n+1) values y_i^(m), find polynomial P_{n,p}(t) such that:
  - P_{n,p}^(m)(x_i) = y_i^(m) for i=0,...,n and m=0,...,p-1
  - Polynomial has degree p(n+1)-1
- **Bell Polynomials B_n(ω; g_1,...,g_n)**:
  - Defined by e^{-ωg}D_t^k e^{ωg} = B_k(ω)
  - Used to express derivatives in terms of function values
- **Interpolation Formula** (equation 3.6):
  ```
  P_{n,p}(t) = Σᵢ Lᵢᵖ(t) Σₘ [(t-xᵢ)ᵐ/m!] yᵢ^(m) · Σᵣ [(t-xᵢ)ʳ/r!] Bᵣ(p; S₁,...,Sᵣ)
  ```
  where S_r(x_i) = (-1)^r(r-1)! Σ_{v≠i} 1/(x_i - x_v)^r

**Applications**:
- Confluent divided differences
- Generalized Cauchy relations
- Foundation for computing Padé approximants with derivative information

**Connection to Gelfgren**: Provides the theoretical machinery for computing rational approximants at interval endpoints when derivative values are given.

## 3. Farouki & Rajan 1987: "Algorithms for Polynomials in Bernstein Form"

**Core Concept**: Algorithms for polynomial operations in Bernstein basis, emphasizing numerical stability.

**Key Ideas**:
- **Bernstein Basis on [0,1]**: b_k^n(x) = (n choose k)(1-x)^{n-k}x^k
- **Superior Numerical Stability**:
  - Bernstein basis has better root condition numbers than power basis
  - **Critical**: Do NOT convert power→Bernstein; formulate directly in Bernstein
  - Explicit conversion amplifies errors and cancels stability benefits

**Essential Algorithms**:

1. **Degree Elevation** (Section 3.2):
   - Single stage: C_k^{n+1} = (k/(n+1))C_{k-1}^n + (1-k/(n+1))C_k^n
   - General r-fold: equation (27)
   - **Degree Reduction**: equations (29)-(30)

2. **Arithmetic Operations** (Section 4):
   - **Addition**: C_k^m = A_k^m ± B_k^m (after degree matching)
   - **Multiplication**: C_k^{m+n} = Σ (m choose j)(n choose k-j)/(m+n choose k) A_j^m B_{k-j}^n
   - **Division**: System of m+1 linear equations (47)

3. **Differentiation & Integration** (Section 5.1):
   - **Derivative**: D_k^{n-1} = n[C_{k+1}^n - C_k^n]
   - **Integral**: I_k^{n+1} = 1/(n+1) Σ_{j=0}^{k-1} C_j^n

4. **Evaluation** (Section 3.4):
   - de Casteljau algorithm (stable, O(n²))
   - Horner-like method (faster but less stable)

5. **Interpolation** (Section 5.3):
   - Vandermonde matrix approach
   - Lagrange polynomial construction

6. **GCD & Resultants** (Section 5.4-5.5):
   - Euclid's algorithm in Bernstein form
   - Sylvester resultant (77)-(78)
   - Bezout form using scaled coefficients

**Scaled Coefficients**: Ĉ_k^n = (n choose k)C_k^n
- Simplifies many operations
- Closer to power form structure
- Recommended for implementation

**Connection to Gelfgren**: The numerator P_n and denominator Q_m polynomials in Gelfgren's rational approximants should be represented in Bernstein form for numerical stability. All polynomial operations (evaluation, differentiation, arithmetic) can be performed directly in this basis.

---

## Implementation Strategy

Based on these papers, the library should:

1. **Core Polynomial Type**: `BernsteinPolynomial<T>` with:
   - Scaled Bernstein coefficients
   - Degree elevation/reduction
   - Arithmetic operations (add, multiply, divide)
   - Differentiation and integration
   - Evaluation via de Casteljau

2. **Rational Function Type**: `RationalFunction<T>`:
   - Numerator and denominator as BernsteinPolynomials
   - Evaluation
   - Differentiation
   - Common factor removal (GCD)

3. **Padé Approximant Construction**:
   - From power series coefficients (solve linear system)
   - Symmetric form at two points (Gelfgren's method)
   - With derivative constraints (using Traub's formulas)

4. **Piecewise Rational Interpolation**:
   - Mesh management (interval partition)
   - Per-subinterval rational approximant construction
   - Continuity enforcement at nodes
   - Error estimation

5. **API Entry Points**:
   - `construct_constraint_system(mesh, degrees)` → Linear system
   - `construct_approximant(mesh, function_data)` → Piecewise rational function
   - `evaluate(piecewise_rational, point)` → Value
   - `solve_bvp(...)` → Boundary value problem solution
