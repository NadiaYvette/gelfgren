# Comprehensive Test Suite Design

## Overview

The test suite validates the Gelfgren library against competing methods and challenging test cases, measuring **accuracy**, **precision**, and **performance**.

---

## Test Categories

### 1. Special Function Approximation

Compare piecewise rational approximants against:
- Reference implementations (MPFR, Boost, scipy.special)
- Chebyshev polynomial approximations
- Minimax polynomial approximations
- Padé approximants (global, not piecewise)

**Test Functions**:

#### Elementary Functions
- `exp(x)` - Exponential
- `log(x)` - Natural logarithm
- `sin(x), cos(x), tan(x)` - Trigonometric
- `sinh(x), cosh(x), tanh(x)` - Hyperbolic
- `arcsin(x), arctan(x)` - Inverse trig

#### Special Functions
- `erf(x), erfc(x)` - Error functions
- `Γ(x)` - Gamma function
- `J₀(x), J₁(x), Jₙ(x)` - Bessel functions
- `I₀(x), K₀(x)` - Modified Bessel functions
- `Ei(x)` - Exponential integral

#### **Jacobi Elliptic Functions** (Primary Challenge Case)
- `sn(u, m)` - Jacobi sine amplitude
- `cn(u, m)` - Jacobi cosine amplitude
- `dn(u, m)` - Jacobi delta amplitude
- `am(u, m)` - Jacobi amplitude

**Why Jacobi elliptic functions are challenging**:
1. **Poles near real line**: `sn(u + iK', m)` has poles at `u = 2nK` (n ∈ ℤ)
2. **Large argument behavior**: For large `u`, Runge phenomenon can occur
3. **Periodicity**: Period `4K(m)` depends on elliptic modulus `m`
4. **Closed forms at rational multiples**: `sn(K/2, m) = 1/√(1+√m)`
5. **Coupling**: sn², cn², dn² satisfy `sn² + cn² = 1`, `dn² + m·sn² = 1`

**Test scenarios**:
```rust
#[test]
fn test_jacobi_sn_small_modulus() {
    // m → 0: sn(u,m) → sin(u)
    // Should match sine approximation quality
}

#[test]
fn test_jacobi_sn_large_argument() {
    // Large u near pole regions
    // Detect if Runge phenomenon occurs
    let u = 100.0;
    let m = 0.5;
    // Compare with reference implementation
}

#[test]
fn test_jacobi_closed_forms() {
    // Verify exact values at special points
    // sn(K, m) = 1
    // sn(0, m) = 0
    // sn(K/2, m) = 1/√(1+√m)
}

#[test]
fn test_jacobi_periodicity() {
    // sn(u + 4K, m) = sn(u, m)
    // Verify approximant respects periodicity
}

#[test]
fn test_jacobi_coupling_relations() {
    // sn² + cn² = 1
    // dn² + m·sn² = 1
    // These should hold to machine precision
}
```

**Metrics**:
- **Accuracy**: Max absolute error vs. reference
- **Precision**: Relative error for small values
- **Conditioning**: Error amplification near poles
- **Runge detection**: Error growth rate vs. polynomial degree

---

### 2. Boundary Value Problems (BVPs)

**Test Problems**:

#### Linear BVPs with Known Solutions

```rust
#[test]
fn test_bvp_harmonic_oscillator() {
    // u'' + u = 0, u(0) = 0, u(π) = 0
    // Exact: u(x) = C·sin(x)
    // Solution: u(x) = 0 or u(x) = sin(x)
}

#[test]
fn test_bvp_exponential() {
    // u'' - u = 0, u(0) = 1, u(1) = e
    // Exact: u(x) = e^x
}

#[test]
fn test_bvp_airy() {
    // u'' - x·u = 0, u(0) = Ai(0), u(1) = Ai(1)
    // Exact: u(x) = Ai(x) (Airy function)
}
```

#### Sturm-Liouville Problems

```rust
#[test]
fn test_sturm_liouville_eigenvalue() {
    // -u'' = λu, u(0) = 0, u(π) = 0
    // Eigenvalues: λₙ = n²
    // Eigenfunctions: uₙ(x) = sin(nx)
}
```

#### Nonlinear BVPs (if supported)

```rust
#[test]
fn test_bvp_bratu() {
    // u'' + λ·e^u = 0, u(0) = 0, u(1) = 0
    // Classic nonlinear BVP with bifurcation
}
```

**Comparison Methods**:
- Finite difference methods
- Shooting methods
- Galerkin methods with polynomial basis
- Spectral methods (Chebyshev collocation)

---

### 3. Initial Value Problems (IVPs)

**Test Problems**:

#### Scalar ODEs

```rust
#[test]
fn test_ivp_exponential() {
    // u' = u, u(0) = 1
    // Exact: u(t) = e^t
}

#[test]
fn test_ivp_logistic() {
    // u' = u(1 - u), u(0) = 0.5
    // Exact: u(t) = 1/(1 + e^(-t))
}

#[test]
fn test_ivp_stiff_problem() {
    // u' = -1000(u - cos(t)) - sin(t), u(0) = 0
    // Stiff ODE, tests numerical stability
}
```

#### Systems (Vector-Valued)

```rust
#[test]
fn test_ivp_harmonic_oscillator_system() {
    // u₁' = u₂, u₂' = -u₁
    // u₁(0) = 1, u₂(0) = 0
    // Exact: u₁(t) = cos(t), u₂(t) = -sin(t)
}

#[test]
fn test_ivp_lorenz_attractor() {
    // Chaotic system - tests robustness
    // σ=10, ρ=28, β=8/3
}

#[test]
fn test_ivp_planetary_orbit() {
    // Two-body problem
    // Tests conservation of energy/momentum
}
```

**Comparison Methods**:
- Runge-Kutta methods (RK4, RK45)
- Adams-Bashforth methods
- BDF methods (for stiff problems)

---

### 4. Space Curves (Parametric Curves)

**Definition**: A space curve is r(t) = (x(t), y(t), z(t)) where each component is approximated.

#### Test Curves

```rust
#[test]
fn test_helix() {
    // r(t) = (cos(t), sin(t), t)
    // Tests periodic + linear components
}

#[test]
fn test_trefoil_knot() {
    // r(t) = (sin(t) + 2sin(2t), cos(t) - 2cos(2t), -sin(3t))
    // Tests multiple frequencies
}

#[test]
fn test_viviani_curve() {
    // Intersection of sphere and cylinder
    // r(t) = (cos²(t), cos(t)sin(t), sin(t))
}

#[test]
fn test_lissajous_3d() {
    // r(t) = (sin(aωt + δ), sin(bωt), sin(cωt))
    // Tests frequency coupling
}
```

**Comparison with NURBS**:

```rust
#[test]
fn test_vs_nurbs_circle() {
    // Approximate circle with both methods
    // Compare control point count, error, evaluation speed
}

#[test]
fn test_vs_nurbs_bezier_curve() {
    // Cubic Bezier vs. piecewise rational
    // Measure approximation quality
}

#[test]
fn test_curvature_preservation() {
    // Check if curvature κ(t) = |r' × r''|/|r'|³ is preserved
}
```

**Metrics**:
- **Arc length error**: ∫|r'(t)| dt vs. true arc length
- **Curvature accuracy**: κ_approx(t) vs. κ_true(t)
- **Torsion accuracy** (for 3D): τ_approx(t) vs. τ_true(t)
- **Control point economy**: Number of control points for same accuracy

---

### 5. Coupling and Vector Functions

For coupled systems where components interact:

```rust
#[test]
fn test_coupled_oscillators() {
    // u₁'' + k₁₁u₁ + k₁₂u₂ = 0
    // u₂'' + k₂₁u₁ + k₂₂u₂ = 0
    // Normal modes depend on coupling
}

#[test]
fn test_jacobi_elliptic_coupling() {
    // sn(u,m), cn(u,m), dn(u,m) satisfy algebraic relations
    // Approximate all three and verify relations hold
}

#[test]
fn test_electromagnetic_field() {
    // E(t), B(t) satisfy Maxwell equations
    // Tests physical coupling constraints
}
```

---

### 6. Runge Phenomenon Detection

**Definition**: Oscillations at interval boundaries when interpolating with high-degree polynomials.

```rust
#[test]
fn test_runge_function() {
    // f(x) = 1/(1 + 25x²) on [-1, 1]
    // Classic case: high-degree polynomials fail
    // Piecewise rational should handle better
}

#[test]
fn test_runge_with_poles() {
    // f(x) = 1/(1 + (x-0.5i)²) on [-1, 1]
    // Poles at x = ±(1/5)i (near real axis)
    // Similar to Jacobi elliptic near poles
}

#[test]
fn test_runge_suppression() {
    // Measure oscillation amplitude near boundaries
    // Compare uniform vs. Chebyshev nodes
    // Compare polynomial vs. rational approximants
}
```

**Detection Method**:
```rust
fn detect_runge_phenomenon<T: Float>(
    approximant: &PiecewiseRational<T>,
    true_function: impl Fn(T) -> T,
    interval: (T, T),
) -> RungeMetrics {
    // Sample heavily near boundaries
    // Compute error ratio: max_boundary_error / max_interior_error
    // Ratio >> 1 indicates Runge phenomenon
}
```

---

## Test Suite Structure

```
tests/
├── unit/
│   ├── bernstein_polynomial.rs
│   ├── rational_function.rs
│   ├── pade_approximant.rs
│   └── mesh.rs
├── integration/
│   ├── special_functions/
│   │   ├── elementary.rs
│   │   ├── bessel.rs
│   │   ├── gamma.rs
│   │   └── jacobi_elliptic.rs    # Primary challenge
│   ├── ode/
│   │   ├── bvp_linear.rs
│   │   ├── bvp_sturm_liouville.rs
│   │   ├── ivp_scalar.rs
│   │   └── ivp_systems.rs
│   ├── space_curves/
│   │   ├── parametric.rs
│   │   └── vs_nurbs.rs
│   └── pathological/
│       ├── runge_phenomenon.rs
│       └── near_poles.rs
└── benchmarks/
    ├── evaluation_speed.rs
    ├── construction_time.rs
    └── memory_usage.rs
```

---

## Metrics and Reporting

### Accuracy Metrics

```rust
pub struct AccuracyMetrics {
    /// Maximum absolute error
    pub max_abs_error: f64,

    /// Maximum relative error
    pub max_rel_error: f64,

    /// RMS error
    pub rms_error: f64,

    /// Error at specific test points
    pub point_errors: Vec<(f64, f64)>, // (x, error)

    /// Condition number estimate
    pub condition_number: f64,
}
```

### Performance Metrics

```rust
pub struct PerformanceMetrics {
    /// Construction time
    pub construction_time: Duration,

    /// Single evaluation time
    pub evaluation_time: Duration,

    /// Memory usage (bytes)
    pub memory_usage: usize,

    /// Throughput (evaluations/sec)
    pub throughput: f64,
}
```

### Comparison Report

```rust
pub struct ComparisonReport {
    method_name: String,
    accuracy: AccuracyMetrics,
    performance: PerformanceMetrics,
}

pub fn generate_comparison_table(
    test_name: &str,
    results: Vec<ComparisonReport>,
) -> String {
    // Generate markdown table:
    // | Method          | Max Error | Time    | Memory  |
    // |-----------------|-----------|---------|---------|
    // | Gelfgren        | 1.2e-12   | 45 μs   | 8.4 KB  |
    // | Chebyshev       | 3.1e-11   | 12 μs   | 4.2 KB  |
    // | Global Padé     | 2.1e-8    | 8 μs    | 1.1 KB  |
    // | NURBS (if curve)| 5.4e-10   | 22 μs   | 12 KB   |
}
```

---

## Reference Implementations

### Libraries to Compare Against

**Rust**:
- `peroxide` - Numerical library with ODE solvers
- `numeric` - Special functions
- `rgsl` - Rust bindings to GSL

**C/C++** (via FFI):
- GNU Scientific Library (GSL) - Special functions, ODEs
- Boost - Math special functions
- MPFR - Arbitrary precision

**Python** (for validation):
- `scipy.special` - Special functions reference
- `scipy.integrate` - ODE solvers (solve_bvp, solve_ivp)
- `mpmath` - Arbitrary precision
- `jacobi` - Jacobi elliptic functions package

**Julia** (gold standard for numerical computing):
- `SpecialFunctions.jl`
- `Elliptic.jl` - Jacobi elliptic functions
- `DifferentialEquations.jl` - ODE/BVP solvers

---

## Automated Test Execution

```rust
// tests/integration/runner.rs

pub fn run_all_comparison_tests() -> TestSuiteReport {
    let mut report = TestSuiteReport::new();

    // Special functions
    report.add_category("Special Functions", vec![
        run_test("exp", test_exp_approximation),
        run_test("Bessel J0", test_bessel_j0),
        run_test("Jacobi sn (small m)", test_jacobi_sn_small),
        run_test("Jacobi sn (large u)", test_jacobi_sn_large),
        run_test("Jacobi sn (near pole)", test_jacobi_sn_pole),
    ]);

    // BVPs
    report.add_category("Boundary Value Problems", vec![
        run_test("Harmonic oscillator", test_bvp_harmonic),
        run_test("Exponential", test_bvp_exponential),
        run_test("Airy function", test_bvp_airy),
    ]);

    // IVPs
    report.add_category("Initial Value Problems", vec![
        run_test("Scalar exponential", test_ivp_exp),
        run_test("Logistic equation", test_ivp_logistic),
        run_test("Harmonic system", test_ivp_harmonic_system),
    ]);

    // Space curves
    report.add_category("Space Curves", vec![
        run_test("Helix", test_curve_helix),
        run_test("vs NURBS circle", test_vs_nurbs_circle),
    ]);

    // Pathological cases
    report.add_category("Pathological Cases", vec![
        run_test("Runge function", test_runge_function),
        run_test("Near pole behavior", test_near_pole),
    ]);

    report.generate_html("test_report.html");
    report
}
```

---

## CI/CD Integration

```yaml
# .github/workflows/comprehensive-tests.yml
name: Comprehensive Test Suite

on: [push, pull_request]

jobs:
  test-suite:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Install dependencies
        run: |
          sudo apt-get install libgsl-dev  # For reference comparisons
          pip install scipy mpmath jacobi   # Python references

      - name: Run unit tests
        run: cargo test --lib

      - name: Run integration tests
        run: cargo test --test '*' -- --nocapture

      - name: Run benchmarks
        run: cargo bench --no-fail-fast

      - name: Generate comparison report
        run: cargo run --bin generate_report

      - name: Upload test report
        uses: actions/upload-artifact@v3
        with:
          name: test-report
          path: test_report.html
```

---

## Success Criteria

For the library to be considered production-ready:

### Accuracy
- ✓ **Special functions**: Error < 10⁻¹² for double precision
- ✓ **Jacobi elliptic**: Closed forms exact to machine precision
- ✓ **Jacobi coupling**: Relations satisfied to < 10⁻¹⁴
- ✓ **ODEs (BVP/IVP)**: Comparable to established solvers
- ✓ **Runge suppression**: No oscillations for degree < mesh size

### Performance
- ✓ **Construction**: Faster than solving equivalent linear system directly
- ✓ **Evaluation**: Comparable to Horner's method (within 2x)
- ✓ **Memory**: O(n·m) per subinterval where n,m are degrees

### Robustness
- ✓ **Near poles**: Graceful degradation, not catastrophic failure
- ✓ **Large arguments**: No overflow in Jacobi elliptic functions
- ✓ **Stiff ODEs**: Stable for moderately stiff problems
- ✓ **Ill-conditioned**: Better than power basis (proven by Farouki-Rajan)

---

## Example Test Output

```
Running Gelfgren Comprehensive Test Suite
==========================================

Special Functions
-----------------
  exp(x) on [-1,1]
    Gelfgren (n=4,m=4):  max_err=1.2e-12  time=45μs  ✓
    Chebyshev (deg 8):   max_err=3.1e-11  time=12μs  ✓
    Global Padé [8/8]:   max_err=2.1e-8   time=8μs   ⚠

  Jacobi sn(u, 0.5) on [0, 4K]
    Gelfgren (n=6,m=6):  max_err=2.4e-13  time=112μs ✓
    scipy.special:       (reference)
    Near pole (u=K):     rel_err=3.1e-12            ✓
    Coupling |sn²+cn²-1|: 4.2e-15                   ✓
    Runge check:         boundary_err/interior_err=1.02 ✓

Boundary Value Problems
-----------------------
  Harmonic oscillator u''+u=0
    Gelfgren:            max_err=5.1e-11  time=2.3ms ✓
    Finite difference:   max_err=1.2e-8   time=0.8ms ✓
    Spectral (Cheby):    max_err=8.4e-13  time=3.1ms ✓

Space Curves
------------
  Helix r(t)=(cos t, sin t, t)
    Gelfgren (per coord):     arc_err=1.3e-11 ✓
    NURBS (cubic, 12 ctrl):   arc_err=2.8e-10 ✓
    Control point ratio: 0.67 (Gelfgren uses fewer)

PASSED: 47/50 tests
WARNINGS: 3 (acceptable error levels)
FAILED: 0
```

---

## Future Additions

- **Multivariate extensions**: Surface approximation (beyond scope initially)
- **Adaptive mesh refinement**: Automatic subinterval subdivision
- **GPU acceleration**: Batch evaluation on GPU
- **Uncertainty quantification**: Rigorous error bounds via interval arithmetic
