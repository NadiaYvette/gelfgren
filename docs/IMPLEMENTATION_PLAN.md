# Gelfgren Library Implementation Plan

## Project Overview

Implementation of Jan Gelfgren's piecewise rational interpolation algorithm (1975) with support from Traub's Lagrange-Hermite interpolation (1964) and Farouki & Rajan's Bernstein polynomial algorithms (1987).

**Core Purpose**: Construct piecewise rational approximants for univariate functions, particularly for solving boundary value problems, by:
1. Building rational approximants (P_n/Q_m) on mesh subintervals
2. Matching function values and derivatives at mesh points
3. Using Bernstein polynomial representation for numerical stability

---

## Architecture Overview

### Core Type Hierarchy

```
BernsteinPolynomial<T>          (Fundamental building block)
    ↓
RationalFunction<T>              (Numerator/Denominator pair)
    ↓
PadeApproximant<T>              (Specific rational function from series/derivatives)
    ↓
PiecewiseRational<T>            (Collection over mesh intervals)
    ↓
BVPSolution<T>                  (Solution to boundary value problem)
```

### Supporting Components

```
Mesh                             (Interval partition management)
BellPolynomial                   (For Traub's formulas)
LinearConstraintSystem           (For BVP formulation)
NewtonSeries                     (For symmetric Padé approximants)
```

---

## Phase 1: Bernstein Polynomial Foundation

**Goal**: Robust Bernstein polynomial implementation with all necessary operations

### 1.1 Core Type: `BernsteinPolynomial<T>`

```rust
pub struct BernsteinPolynomial<T> {
    /// Scaled Bernstein coefficients: Ĉₖⁿ = (n choose k) Cₖⁿ
    coefficients: Vec<T>,
    /// Degree n
    degree: usize,
    /// Interval [a, b] (defaults to [0, 1])
    interval: (T, T),
}
```

**Key Features**:
- Generic over `T: Float` (f32, f64, or eventually complex)
- Scaled coefficients for numerical stability
- Support for arbitrary interval [a, b]

### 1.2 Essential Operations

**Construction**:
```rust
impl<T: Float> BernsteinPolynomial<T> {
    /// Create from scaled Bernstein coefficients on [0,1]
    pub fn new(coefficients: Vec<T>) -> Self;

    /// Create from scaled coefficients on [a, b]
    pub fn with_interval(coefficients: Vec<T>, interval: (T, T)) -> Self;

    /// Create from power form coefficients (avoid when possible)
    pub fn from_power_form(power_coeffs: &[T]) -> Self;

    /// Zero polynomial of degree n
    pub fn zeros(degree: usize) -> Self;

    /// Constant polynomial
    pub fn constant(value: T, degree: usize) -> Self;
}
```

**Degree Operations** (Farouki-Rajan Section 3.2):
```rust
/// Single degree elevation (equation 24)
pub fn degree_elevate(&self) -> Self;

/// r-fold degree elevation (equation 27)
pub fn degree_elevate_by(&self, r: usize) -> Self;

/// Degree reduction (equations 29-30)
/// Returns None if reduction not possible
pub fn degree_reduce(&self) -> Option<Self>;

/// r-fold degree reduction
pub fn degree_reduce_by(&self, r: usize) -> Option<Self>;

/// Get actual degree (considering trailing zeros)
pub fn actual_degree(&self) -> usize;
```

**Arithmetic** (Farouki-Rajan Section 4):
```rust
/// Addition (equation 41)
pub fn add(&self, other: &Self) -> Result<Self>;

/// Subtraction
pub fn sub(&self, other: &Self) -> Result<Self>;

/// Scalar multiplication
pub fn scale(&self, scalar: T) -> Self;

/// Polynomial multiplication (equation 44)
pub fn multiply(&self, other: &Self) -> Self;

/// Polynomial division (equation 47)
/// Returns (quotient, remainder)
pub fn divide(&self, divisor: &Self) -> Result<(Self, Self)>;
```

**Calculus** (Farouki-Rajan Section 5.1):
```rust
/// Differentiation (equation 57)
pub fn derivative(&self) -> Self;

/// k-th derivative
pub fn derivative_nth(&self, k: usize) -> Self;

/// Indefinite integral (equation 58)
pub fn integral(&self) -> Self;

/// Definite integral on [0, 1] (equation 59)
pub fn definite_integral(&self) -> T;
```

**Evaluation** (Farouki-Rajan Section 3.4):
```rust
/// Evaluate at point using de Casteljau algorithm
pub fn evaluate(&self, x: T) -> T;

/// Evaluate derivative at point
pub fn evaluate_derivative(&self, x: T, order: usize) -> T;

/// Evaluate function and first k derivatives
pub fn evaluate_with_derivatives(&self, x: T, k: usize) -> Vec<T>;
```

**Utility**:
```rust
/// Greatest common divisor with another polynomial
pub fn gcd(&self, other: &Self) -> Self;

/// Remove common factors
pub fn simplify(&mut self);

/// L2 norm on [0, 1] (equation 34)
pub fn norm_l2(&self) -> T;

/// Convert to different interval
pub fn transform_interval(&self, new_interval: (T, T)) -> Self;
```

### 1.3 Helper Functions

**Binomial Coefficients** (Farouki-Rajan Section 3.1):
```rust
/// Compute binomial coefficient using prime factorization
/// to avoid overflow
pub fn binomial(n: u32, k: u32) -> u64;

/// Precompute binomial coefficients up to degree n
pub struct BinomialCache {
    max_degree: usize,
    cache: Vec<Vec<u64>>,
}
```

---

## Phase 2: Rational Functions

**Goal**: Rational function type with proper handling of numerator/denominator

### 2.1 Core Type: `RationalFunction<T>`

```rust
pub struct RationalFunction<T> {
    numerator: BernsteinPolynomial<T>,
    denominator: BernsteinPolynomial<T>,
    /// Interval where defined
    interval: (T, T),
}
```

### 2.2 Operations

```rust
impl<T: Float> RationalFunction<T> {
    /// Create from numerator and denominator
    /// Automatically removes common factors
    pub fn new(
        numerator: BernsteinPolynomial<T>,
        denominator: BernsteinPolynomial<T>,
    ) -> Result<Self>;

    /// Create without simplification
    pub fn new_unsimplified(
        numerator: BernsteinPolynomial<T>,
        denominator: BernsteinPolynomial<T>,
    ) -> Result<Self>;

    /// Evaluate at point
    pub fn evaluate(&self, x: T) -> Result<T>;

    /// Evaluate derivative (quotient rule)
    pub fn evaluate_derivative(&self, x: T, order: usize) -> Result<T>;

    /// Find poles (zeros of denominator)
    pub fn find_poles(&self) -> Vec<T>;

    /// Check if denominator is zero at point
    pub fn is_defined_at(&self, x: T) -> bool;

    /// Simplify by removing common factors
    pub fn simplify(&mut self);

    /// Get degrees (n, m)
    pub fn degrees(&self) -> (usize, usize);
}
```

---

## Phase 3: Padé Approximants

**Goal**: Construct rational approximants from power series or derivative data

### 3.1 From Power Series (Gelfgren Section 1)

```rust
pub struct PadeApproximant<T> {
    rational: RationalFunction<T>,
    degrees: (usize, usize), // (n, m)
}

impl<T: Float> PadeApproximant<T> {
    /// Construct [n,m] Padé approximant from power series coefficients
    /// f(z) = a₀ + a₁z + a₂z² + ...
    ///
    /// Solves linear system (Gelfgren equation 1.4):
    /// Σⱼ aⱼ₋ₖ βₖ = { aⱼ if 0≤j≤n; 0 if n+1≤j≤n+m }
    pub fn from_series(
        series_coeffs: &[T],
        n: usize, // numerator degree
        m: usize, // denominator degree
    ) -> Result<Self>;

    /// Construct symmetric [n,m] Padé approximant at two points z₀, z₁
    /// using Newton series (Gelfgren Section 2)
    ///
    /// Newton basis: wₖ(z) = (z-z₀)^{k-⌊k/2⌋}(z-z₁)^{⌊k/2⌋}
    pub fn symmetric_at_two_points(
        z0: T,
        z1: T,
        n: usize,
        m: usize,
        function_evaluator: impl Fn(T) -> T,
    ) -> Result<Self>;
}
```

### 3.2 Linear System Solver

```rust
/// Solve the Padé linear system to find denominator coefficients βₖ
/// then compute numerator coefficients αₖ
fn solve_pade_system<T: Float>(
    series_coeffs: &[T],
    n: usize,
    m: usize,
) -> Result<(Vec<T>, Vec<T>)>; // (numerator, denominator)
```

---

## Phase 4: Traub's Lagrange-Hermite Interpolation

**Goal**: Interpolate with function and derivative values at multiple points

### 4.1 Bell Polynomials (Traub Section 2)

```rust
/// Bell polynomial Bₙ(ω; g₁, ..., gₙ)
/// Used in Traub's interpolation formula
pub struct BellPolynomial<T> {
    degree: usize,
    // Cached values or computed on demand
}

impl<T: Float> BellPolynomial<T> {
    /// Compute Bₙ(ω; g₁, ..., gₙ) using explicit formula (Traub equation 2.2)
    pub fn evaluate(
        n: usize,
        omega: T,
        g_values: &[T],
    ) -> T;

    /// Compute Uₙ,ₖ coefficients (Traub equation 2.2)
    pub fn compute_u_coefficients(
        n: usize,
        k: usize,
        g_values: &[T],
    ) -> T;
}
```

### 4.2 Lagrange-Hermite Interpolation (Traub Section 3)

```rust
/// Data for Lagrange-Hermite interpolation
pub struct HermiteData<T> {
    /// Interpolation points xᵢ
    points: Vec<T>,
    /// Values: yᵢ^(m) for i=0..n, m=0..p-1
    /// values[i][m] = yᵢ^(m)
    values: Vec<Vec<T>>,
}

/// Construct polynomial interpolating function and derivatives
/// Returns polynomial of degree p(n+1) - 1
///
/// Uses Traub's formula (equation 3.7):
/// Pₙ,ₚ(t) = Σᵢ Lᵢᵖ(t) Σₘ [(t-xᵢ)ᵐ/m!] yᵢ^(m) · Σᵣ [(t-xᵢ)ʳ/r!] Bᵣ(p; S₁,...,Sᵣ)
pub fn lagrange_hermite_interpolation<T: Float>(
    data: &HermiteData<T>,
    p: usize, // multiplicity at each point
) -> Result<BernsteinPolynomial<T>>;

/// Helper: Compute Sᵣ(xᵢ) = (-1)^r(r-1)! Σ_{v≠i} 1/(xᵢ - xᵥ)^r
fn compute_s_values<T: Float>(
    points: &[T],
    i: usize,
    max_r: usize,
) -> Vec<T>;

/// Helper: Compute Lagrange basis polynomials Lᵢ(t)
fn lagrange_basis<T: Float>(
    points: &[T],
    i: usize,
) -> BernsteinPolynomial<T>;
```

---

## Phase 5: Gelfgren's Piecewise Rational Interpolation

**Goal**: Main algorithm - construct piecewise rational approximants on mesh

### 5.1 Mesh Management

```rust
/// Partition of interval [a, b]
pub struct Mesh<T> {
    /// Interval endpoints: a = x₀ < x₁ < ... < xₙ = b
    points: Vec<T>,
}

impl<T: Float> Mesh<T> {
    /// Create uniform mesh with N subintervals
    pub fn uniform(a: T, b: T, n_intervals: usize) -> Self;

    /// Create from explicit points
    pub fn from_points(points: Vec<T>) -> Result<Self>;

    /// Get subinterval [xⱼ₋₁, xⱼ]
    pub fn subinterval(&self, j: usize) -> (T, T);

    /// Get all subintervals
    pub fn subintervals(&self) -> Vec<(T, T)>;

    /// Number of subintervals
    pub fn num_intervals(&self) -> usize;

    /// Number of mesh points
    pub fn num_points(&self) -> usize;

    /// Find which subinterval contains point x
    pub fn locate_interval(&self, x: T) -> Option<usize>;
}
```

### 5.2 Newton Series (Gelfgren equation 2.2)

```rust
/// Newton series basis functions for symmetric Padé approximant
/// wₖ(z) = (z-z₀)^{k-⌊k/2⌋}(z-z₁)^{⌊k/2⌋}
pub struct NewtonBasis<T> {
    z0: T,
    z1: T,
}

impl<T: Float> NewtonBasis<T> {
    pub fn new(z0: T, z1: T) -> Self;

    /// Evaluate k-th basis function
    pub fn evaluate(&self, k: usize, z: T) -> T;

    /// Compute Newton series coefficients from function
    pub fn compute_coefficients(
        &self,
        max_degree: usize,
        f: impl Fn(T) -> T,
    ) -> Vec<T>;
}
```

### 5.3 Piecewise Rational Function

```rust
/// Piecewise rational function Sₙ,ₘ(τ)(z)
pub struct PiecewiseRational<T> {
    /// Mesh partition
    mesh: Mesh<T>,
    /// Rational approximant on each subinterval
    /// approximants[j] is for interval [xⱼ₋₁, xⱼ]
    approximants: Vec<RationalFunction<T>>,
    /// Polynomial degrees (n, m)
    degrees: (usize, usize),
    /// Continuity: number of continuous derivatives (p-1)
    continuity: usize,
}

impl<T: Float> PiecewiseRational<T> {
    /// Evaluate at point x
    /// Finds appropriate subinterval and evaluates local approximant
    pub fn evaluate(&self, x: T) -> Result<T>;

    /// Evaluate k-th derivative
    pub fn evaluate_derivative(&self, x: T, k: usize) -> Result<T>;

    /// Check continuity at mesh points
    pub fn verify_continuity(&self) -> Result<()>;

    /// Get approximant for specific subinterval
    pub fn get_approximant(&self, interval_index: usize) -> &RationalFunction<T>;
}
```

### 5.4 Construction Algorithms

```rust
/// Construct piecewise rational approximant using Gelfgren's method
///
/// For each subinterval Δⱼ = [xⱼ₋₁, xⱼ]:
/// 1. Construct symmetric [n,m] Padé approximant at endpoints
/// 2. Use Newton series basis
/// 3. Match function values at endpoints
pub fn construct_piecewise_rational<T: Float>(
    mesh: &Mesh<T>,
    n: usize, // numerator degree
    m: usize, // denominator degree
    function: impl Fn(T) -> T,
) -> Result<PiecewiseRational<T>>;

/// Construct with derivative matching
/// Uses Traub's formulas to incorporate derivative information
pub fn construct_with_derivatives<T: Float>(
    mesh: &Mesh<T>,
    n: usize,
    m: usize,
    p: usize, // order of derivative matching
    function_data: &HermiteData<T>,
) -> Result<PiecewiseRational<T>>;
```

---

## Phase 6: Linear Constraint Systems for BVPs

**Goal**: API for constructing linear systems for boundary value problems

### 6.1 Constraint System

```rust
/// Linear constraint system for BVP formulation
pub struct LinearConstraintSystem<T> {
    /// Coefficient matrix A
    matrix: Vec<Vec<T>>,
    /// Right-hand side b
    rhs: Vec<T>,
    /// Variable names/indices for debugging
    variables: Vec<String>,
}

impl<T: Float> LinearConstraintSystem<T> {
    /// Create empty system
    pub fn new() -> Self;

    /// Add constraint (row in system)
    pub fn add_constraint(&mut self, coefficients: Vec<T>, rhs_value: T);

    /// Solve system Ax = b
    pub fn solve(&self) -> Result<Vec<T>>;

    /// Export to matrix form for external solver
    pub fn to_matrix(&self) -> (Vec<Vec<T>>, Vec<T>);
}
```

### 6.2 BVP Formulation

```rust
/// Boundary value problem specification
pub struct BoundaryValueProblem<T> {
    /// ODE: L[u] = f where L is differential operator
    /// Represented as linear constraints on polynomial coefficients
    operator: Box<dyn Fn(&BernsteinPolynomial<T>) -> BernsteinPolynomial<T>>,

    /// Right-hand side function
    rhs: Box<dyn Fn(T) -> T>,

    /// Boundary conditions at endpoints
    boundary_conditions: Vec<BoundaryCondition<T>>,

    /// Mesh for approximation
    mesh: Mesh<T>,
}

pub enum BoundaryCondition<T> {
    /// Value condition: u(x) = value
    Value { point: T, value: T },
    /// Derivative condition: u^(k)(x) = value
    Derivative { point: T, order: usize, value: T },
}

/// Construct linear constraint system from BVP
/// Returns system where solution gives polynomial coefficients
pub fn formulate_bvp_constraints<T: Float>(
    bvp: &BoundaryValueProblem<T>,
    n: usize, // polynomial degree
    m: usize, // denominator degree if rational
) -> Result<LinearConstraintSystem<T>>;
```

---

## Phase 7: API Entry Points

**Goal**: High-level API matching user requirements

### 7.1 Primary API

```rust
/// API Entry Point 1: Construct linear constraint system from mesh
///
/// Returns system Ax = b where solution x gives polynomial coefficients
/// for rational approximant on each subinterval
pub fn construct_constraint_system<T: Float>(
    mesh: &Mesh<T>,
    numerator_degree: usize,
    denominator_degree: usize,
) -> Result<LinearConstraintSystem<T>>;

/// API Entry Point 2: Construct approximant from function and derivatives
///
/// Given mesh points and function/derivative values, constructs
/// piecewise rational approximant directly
pub fn construct_approximant<T: Float>(
    mesh: &Mesh<T>,
    numerator_degree: usize,
    denominator_degree: usize,
    function_data: &HermiteData<T>,
) -> Result<PiecewiseRational<T>>;

/// API Entry Point 3: Solve boundary value problem
///
/// Uses piecewise rational approximation to solve ODE with boundary conditions
pub fn solve_boundary_value_problem<T: Float>(
    bvp: &BoundaryValueProblem<T>,
    mesh: &Mesh<T>,
    numerator_degree: usize,
    denominator_degree: usize,
) -> Result<PiecewiseRational<T>>;

/// Evaluate piecewise rational function at point
pub fn evaluate<T: Float>(
    approximant: &PiecewiseRational<T>,
    point: T,
) -> Result<T>;
```

### 7.2 Convenience Functions

```rust
/// Approximate special function (for testing/validation)
pub fn approximate_special_function<T: Float>(
    function: impl Fn(T) -> T,
    interval: (T, T),
    n_intervals: usize,
    numerator_degree: usize,
    denominator_degree: usize,
) -> Result<PiecewiseRational<T>>;

/// Compute error between approximant and true function
pub fn compute_approximation_error<T: Float>(
    approximant: &PiecewiseRational<T>,
    true_function: impl Fn(T) -> T,
    n_samples: usize,
) -> Result<T>; // Returns max absolute error
```

---

## Phase 8: FFI Layer (C Bindings)

**Goal**: Expose functionality to C and other languages

### 8.1 Opaque Handle Pattern

```rust
/// Opaque handle for piecewise rational function
pub struct GelfgrenApproximant {
    _private: [u8; 0],
}

/// Create approximant from function values
#[no_mangle]
pub extern "C" fn gelfgren_create_approximant(
    mesh_points: *const f64,
    n_points: usize,
    function_values: *const f64,
    numerator_degree: usize,
    denominator_degree: usize,
    out_handle: *mut *mut GelfgrenApproximant,
) -> i32;

/// Evaluate approximant at point
#[no_mangle]
pub extern "C" fn gelfgren_evaluate(
    handle: *const GelfgrenApproximant,
    x: f64,
    out_value: *mut f64,
) -> i32;

/// Free approximant
#[no_mangle]
pub extern "C" fn gelfgren_free(
    handle: *mut GelfgrenApproximant,
) -> i32;
```

### 8.2 Error Handling (Same as before)

Thread-local error storage, error codes, `gelfgren_last_error_message()` etc.

### 8.3 Matrix/Vector Interface

```rust
/// Create constraint system
#[no_mangle]
pub extern "C" fn gelfgren_construct_constraints(
    mesh_points: *const f64,
    n_points: usize,
    n_degree: usize,
    m_degree: usize,
    out_system: *mut *mut GelfgrenConstraintSystem,
) -> i32;

/// Get constraint matrix dimensions
#[no_mangle]
pub extern "C" fn gelfgren_system_dimensions(
    system: *const GelfgrenConstraintSystem,
    out_rows: *mut usize,
    out_cols: *mut usize,
) -> i32;

/// Get constraint matrix data (row-major)
#[no_mangle]
pub extern "C" fn gelfgren_get_matrix(
    system: *const GelfgrenConstraintSystem,
    out_data: *mut f64,
    capacity: usize,
) -> i32;
```

---

## Phase 9: Testing & Validation

### 9.1 Unit Tests

```rust
#[cfg(test)]
mod tests {
    // Bernstein polynomial tests
    #[test] fn test_degree_elevation();
    #[test] fn test_polynomial_multiplication();
    #[test] fn test_de_casteljau_evaluation();

    // Rational function tests
    #[test] fn test_rational_simplification();
    #[test] fn test_rational_evaluation();

    // Padé approximant tests
    #[test] fn test_pade_exp(); // e^x approximation
    #[test] fn test_pade_sin(); // sin(x) approximation

    // Lagrange-Hermite tests
    #[test] fn test_hermite_interpolation();
    #[test] fn test_bell_polynomial_evaluation();

    // Piecewise tests
    #[test] fn test_piecewise_continuity();
    #[test] fn test_piecewise_evaluation();
}
```

### 9.2 Integration Tests

```rust
// Test against known special functions
#[test]
fn test_approximate_erf() {
    // Approximate error function erf(x)
    // Compare against known values
}

#[test]
fn test_approximate_bessel() {
    // Approximate Bessel function J₀(x)
}

#[test]
fn test_simple_bvp() {
    // Solve u'' = -u with u(0)=0, u(π)=0
    // Should get u(x) = sin(x)
}
```

### 9.3 Validation Strategy

1. **Compare with Padé tables**: Test Padé approximants against published tables
2. **Special function approximation**: Validate against scipy, mpmath, etc.
3. **Numerical stability**: Verify Bernstein form gives better conditioning
4. **Convergence tests**: Verify error decreases as mesh is refined
5. **BVP validation**: Compare with established ODE solvers

---

## Implementation Order

### Stage 1: Foundation (Week 1)
1. ✅ Project structure, error handling (already done)
2. 🔲 `BernsteinPolynomial<T>` core type
3. 🔲 Degree elevation/reduction
4. 🔲 Arithmetic operations (add, multiply)
5. 🔲 de Casteljau evaluation
6. 🔲 Unit tests for polynomials

### Stage 2: Calculus & Rational Functions (Week 1-2)
7. 🔲 Differentiation and integration
8. 🔲 `RationalFunction<T>` type
9. 🔲 GCD and simplification
10. 🔲 Unit tests for rational functions

### Stage 3: Padé Approximants (Week 2)
11. 🔲 Linear system solver interface
12. 🔲 Padé from power series
13. 🔲 Symmetric Padé (Newton series basis)
14. 🔲 Tests against known Padé approximants

### Stage 4: Hermite Interpolation (Week 2-3)
15. 🔲 Bell polynomial computation
16. 🔲 Lagrange-Hermite interpolation
17. 🔲 Integration with Padé construction
18. 🔲 Interpolation tests

### Stage 5: Piecewise Construction (Week 3)
19. 🔲 `Mesh` type
20. 🔲 `PiecewiseRational<T>` type
21. 🔲 Gelfgren's algorithm implementation
22. 🔲 Continuity verification
23. 🔲 Integration tests

### Stage 6: BVP Support (Week 3-4)
24. 🔲 `LinearConstraintSystem` type
25. 🔲 BVP formulation
26. 🔲 High-level API functions
27. 🔲 BVP validation tests

### Stage 7: FFI Layer (Week 4)
28. 🔲 C header generation (cbindgen)
29. 🔲 Opaque handle pattern
30. 🔲 FFI functions for all APIs
31. 🔲 C integration tests

### Stage 8: Multi-Language Bindings (Week 5+)
32. 🔲 Python (PyO3)
33. 🔲 Java (JNI)
34. 🔲 R (extendr)
35. 🔲 Other languages as needed

---

## Dependencies

```toml
[dependencies]
num-traits = "0.2"       # Generic numeric traits
approx = "0.5"           # Floating-point comparisons

# Linear algebra (for solving systems)
nalgebra = "0.32"        # Or ndarray with ndarray-linalg

# Optional: arbitrary precision
rug = { version = "1.19", optional = true }

[dev-dependencies]
approx = "0.5"
proptest = "1.4"         # Property-based testing
criterion = "0.5"        # Benchmarking
```

---

## Success Criteria

### Mathematical Correctness
- ✓ Padé approximants match published tables
- ✓ Lagrange-Hermite matches analytic interpolants
- ✓ Piecewise rational has correct continuity
- ✓ BVP solutions match known solutions

### Numerical Stability
- ✓ Bernstein form more stable than power form
- ✓ Well-conditioned linear systems
- ✓ Accurate evaluation near boundaries

### API Usability
- ✓ Clear entry points for each use case
- ✓ Comprehensive error messages
- ✓ Examples for common scenarios

### Multi-Language Support
- ✓ C FFI with clean interface
- ✓ At least Python, Java, R bindings working
- ✓ Documentation for each language

---

## Open Questions

1. **Linear System Solver**: Use nalgebra, ndarray-linalg, or external (LAPACK)?
2. **Precision**: f32, f64, or also arbitrary precision (rug)?
3. **Parallel**: Parallelize per-subinterval construction?
4. **Interval Arithmetic**: Include for rigorous error bounds?
5. **Root Finding**: Need polynomial root finding for poles?

---

## References

- [1] Gelfgren, J. (1975). "Piecewise Rational Interpolation"
- [2] Traub, J.F. (1964). "On Lagrange-Hermite Interpolation"
- [3] Farouki, R.T. & Rajan, V.T. (1987). "Algorithms for Polynomials in Bernstein Form"
