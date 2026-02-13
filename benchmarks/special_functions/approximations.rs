/// Benchmark: Special Function Approximations
///
/// Compares piecewise rational approximants vs polynomial splines for:
/// 1. Exponential function: exp(x) on various intervals
/// 2. Trigonometric functions: sin(x), cos(x), tan(x)
/// 3. Error function: erf(x)
/// 4. Bessel functions: J_0(x), Y_0(x)
/// 5. Logarithm: log(1+x)
///
/// Tests the hypothesis that rational approximants can achieve
/// comparable accuracy to polynomial splines with coarser meshes.

use std::f64::consts::{E, PI};

/// Test case for special function approximation
pub trait SpecialFunction {
    /// Name of the function
    fn name(&self) -> &'static str;

    /// Function to approximate
    fn eval(&self, x: f64) -> f64;

    /// Derivative (for error analysis)
    fn derivative(&self, x: f64) -> f64;

    /// Domain bounds
    fn domain(&self) -> (f64, f64);

    /// Typical mesh sizes for convergence study
    fn mesh_sizes(&self) -> Vec<usize> {
        vec![4, 8, 16, 32, 64, 128]
    }
}

/// Exponential function: e^x on [-1, 1]
pub struct Exponential {
    pub name: &'static str,
    pub a: f64,
    pub b: f64,
}

impl Exponential {
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            name: "Exponential e^x",
            a,
            b,
        }
    }

    /// Standard exponential on [-1, 1]
    pub fn standard() -> Self {
        Self::new(-1.0, 1.0)
    }

    /// Wide range exponential on [-5, 5]
    pub fn wide() -> Self {
        Self::new(-5.0, 5.0)
    }
}

impl SpecialFunction for Exponential {
    fn name(&self) -> &'static str {
        self.name
    }

    fn eval(&self, x: f64) -> f64 {
        x.exp()
    }

    fn derivative(&self, x: f64) -> f64 {
        x.exp()
    }

    fn domain(&self) -> (f64, f64) {
        (self.a, self.b)
    }
}

/// Sine function: sin(x)
pub struct Sine {
    pub name: &'static str,
    pub a: f64,
    pub b: f64,
}

impl Sine {
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            name: "Sine sin(x)",
            a,
            b,
        }
    }

    /// Standard sine on [0, 2π]
    pub fn standard() -> Self {
        Self::new(0.0, 2.0 * PI)
    }

    /// Multiple periods on [0, 10π]
    pub fn multiple_periods() -> Self {
        Self::new(0.0, 10.0 * PI)
    }
}

impl SpecialFunction for Sine {
    fn name(&self) -> &'static str {
        self.name
    }

    fn eval(&self, x: f64) -> f64 {
        x.sin()
    }

    fn derivative(&self, x: f64) -> f64 {
        x.cos()
    }

    fn domain(&self) -> (f64, f64) {
        (self.a, self.b)
    }
}

/// Cosine function: cos(x)
pub struct Cosine {
    pub name: &'static str,
    pub a: f64,
    pub b: f64,
}

impl Cosine {
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            name: "Cosine cos(x)",
            a,
            b,
        }
    }

    pub fn standard() -> Self {
        Self::new(0.0, 2.0 * PI)
    }
}

impl SpecialFunction for Cosine {
    fn name(&self) -> &'static str {
        self.name
    }

    fn eval(&self, x: f64) -> f64 {
        x.cos()
    }

    fn derivative(&self, x: f64) -> f64 {
        -x.sin()
    }

    fn domain(&self) -> (f64, f64) {
        (self.a, self.b)
    }
}

/// Tangent function: tan(x) on (-π/2, π/2)
pub struct Tangent {
    pub name: &'static str,
    pub a: f64,
    pub b: f64,
}

impl Tangent {
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            name: "Tangent tan(x)",
            a,
            b,
        }
    }

    /// Restricted domain to avoid poles
    pub fn standard() -> Self {
        Self::new(-PI / 3.0, PI / 3.0)
    }
}

impl SpecialFunction for Tangent {
    fn name(&self) -> &'static str {
        self.name
    }

    fn eval(&self, x: f64) -> f64 {
        x.tan()
    }

    fn derivative(&self, x: f64) -> f64 {
        let cos_x = x.cos();
        1.0 / (cos_x * cos_x)
    }

    fn domain(&self) -> (f64, f64) {
        (self.a, self.b)
    }
}

/// Error function: erf(x)
/// Using approximation erf(x) ≈ sign(x) * sqrt(1 - exp(-x^2 * (4/π + a*x^2) / (1 + a*x^2)))
/// where a = 8(π-3) / (3π(4-π))
pub struct ErrorFunction {
    pub name: &'static str,
    pub a: f64,
    pub b: f64,
}

impl ErrorFunction {
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            name: "Error function erf(x)",
            a,
            b,
        }
    }

    pub fn standard() -> Self {
        Self::new(-3.0, 3.0)
    }

    /// High-precision error function approximation
    fn erf_impl(&self, x: f64) -> f64 {
        // Use series expansion for small |x|
        if x.abs() < 0.5 {
            let x2 = x * x;
            let x3 = x2 * x;
            let x5 = x3 * x2;
            let x7 = x5 * x2;

            // erf(x) = 2/√π * (x - x³/3 + x⁵/10 - x⁷/42 + ...)
            let sqrt_pi = PI.sqrt();
            return (2.0 / sqrt_pi) * (x - x3 / 3.0 + x5 / 10.0 - x7 / 42.0);
        }

        // Use complementary error function for large |x|
        let sign = x.signum();
        let abs_x = x.abs();

        // Rational approximation
        let a1 = 0.254829592;
        let a2 = -0.284496736;
        let a3 = 1.421413741;
        let a4 = -1.453152027;
        let a5 = 1.061405429;
        let p = 0.3275911;

        let t = 1.0 / (1.0 + p * abs_x);
        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t3 * t;
        let t5 = t4 * t;

        let erfc = (a1 * t + a2 * t2 + a3 * t3 + a4 * t4 + a5 * t5) * (-abs_x * abs_x).exp();

        sign * (1.0 - erfc)
    }
}

impl SpecialFunction for ErrorFunction {
    fn name(&self) -> &'static str {
        self.name
    }

    fn eval(&self, x: f64) -> f64 {
        self.erf_impl(x)
    }

    fn derivative(&self, x: f64) -> f64 {
        // d/dx erf(x) = 2/√π * e^(-x²)
        (2.0 / PI.sqrt()) * (-x * x).exp()
    }

    fn domain(&self) -> (f64, f64) {
        (self.a, self.b)
    }
}

/// Bessel function of the first kind: J_0(x)
/// Using polynomial approximation for small x and asymptotic expansion for large x
pub struct BesselJ0 {
    pub name: &'static str,
    pub a: f64,
    pub b: f64,
}

impl BesselJ0 {
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            name: "Bessel J_0(x)",
            a,
            b,
        }
    }

    pub fn standard() -> Self {
        Self::new(0.0, 10.0)
    }

    /// Polynomial approximation for J_0(x)
    fn j0_impl(&self, x: f64) -> f64 {
        let abs_x = x.abs();

        if abs_x < 8.0 {
            // Polynomial approximation for |x| < 8
            let y = x * x;
            let ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
                + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
            let ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
                + y * (59272.64853 + y * (267.8532712 + y))));

            ans1 / ans2
        } else {
            // Asymptotic expansion for |x| >= 8
            let z = 8.0 / abs_x;
            let y = z * z;
            let xx = abs_x - 0.785398164;

            let ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
            let ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                - y * 0.934945152e-7)));

            (0.636619772 / abs_x.sqrt()) * (xx.cos() * ans1 - z * xx.sin() * ans2)
        }
    }
}

impl SpecialFunction for BesselJ0 {
    fn name(&self) -> &'static str {
        self.name
    }

    fn eval(&self, x: f64) -> f64 {
        self.j0_impl(x)
    }

    fn derivative(&self, x: f64) -> f64 {
        // d/dx J_0(x) = -J_1(x)
        // Simplified numerical derivative for now
        let h = 1e-8;
        (self.j0_impl(x + h) - self.j0_impl(x - h)) / (2.0 * h)
    }

    fn domain(&self) -> (f64, f64) {
        (self.a, self.b)
    }
}

/// Natural logarithm: log(1+x) on [0, e-1]
pub struct Logarithm {
    pub name: &'static str,
    pub a: f64,
    pub b: f64,
}

impl Logarithm {
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            name: "Logarithm log(1+x)",
            a,
            b,
        }
    }

    pub fn standard() -> Self {
        Self::new(0.0, E - 1.0)
    }
}

impl SpecialFunction for Logarithm {
    fn name(&self) -> &'static str {
        self.name
    }

    fn eval(&self, x: f64) -> f64 {
        (1.0 + x).ln()
    }

    fn derivative(&self, x: f64) -> f64 {
        1.0 / (1.0 + x)
    }

    fn domain(&self) -> (f64, f64) {
        (self.a, self.b)
    }
}

/// Rational function: 1/(1+25x²) (Runge's function)
pub struct RungsFunction {
    pub name: &'static str,
    pub a: f64,
    pub b: f64,
}

impl RungsFunction {
    pub fn new(a: f64, b: f64) -> Self {
        Self {
            name: "Runge's function 1/(1+25x²)",
            a,
            b,
        }
    }

    pub fn standard() -> Self {
        Self::new(-1.0, 1.0)
    }
}

impl SpecialFunction for RungsFunction {
    fn name(&self) -> &'static str {
        self.name
    }

    fn eval(&self, x: f64) -> f64 {
        1.0 / (1.0 + 25.0 * x * x)
    }

    fn derivative(&self, x: f64) -> f64 {
        let denom = 1.0 + 25.0 * x * x;
        -50.0 * x / (denom * denom)
    }

    fn domain(&self) -> (f64, f64) {
        (self.a, self.b)
    }
}

/// Approximation configuration
#[derive(Debug, Clone)]
pub struct ApproximationConfig {
    pub num_intervals: usize,
    pub polynomial_degree: usize,
    pub rational_numerator_degree: usize,
    pub rational_denominator_degree: usize,
}

impl ApproximationConfig {
    /// Create convergence sequence for approximation study
    pub fn convergence_sequence(max_intervals: usize) -> Vec<Self> {
        let mut configs = Vec::new();
        let mut n = 4;
        while n <= max_intervals {
            configs.push(ApproximationConfig {
                num_intervals: n,
                polynomial_degree: 3, // Cubic splines
                rational_numerator_degree: 2,
                rational_denominator_degree: 2, // [2/2] Padé
            });
            n *= 2;
        }
        configs
    }

    /// Degrees of freedom for polynomial spline
    pub fn polynomial_dof(&self) -> usize {
        // Cubic spline: n intervals with C1 continuity
        self.num_intervals + 3
    }

    /// Degrees of freedom for rational approximant
    pub fn rational_dof(&self) -> usize {
        // [m/n] rational on each interval
        self.num_intervals * (self.rational_numerator_degree + self.rational_denominator_degree + 2)
    }
}

/// Error metrics for approximation quality
#[derive(Debug, Clone)]
pub struct ApproximationError {
    pub l2_error: f64,
    pub l_inf_error: f64,
    pub h1_seminorm_error: f64,
    pub relative_l2_error: f64,
    pub relative_l_inf_error: f64,
    pub num_eval_points: usize,
    pub mesh_size: f64,
}

impl ApproximationError {
    /// Compute error metrics for approximation
    pub fn compute<F, G>(
        approx: F,
        exact: G,
        a: f64,
        b: f64,
        num_points: usize,
    ) -> Self
    where
        F: Fn(f64) -> f64,
        G: Fn(f64) -> f64,
    {
        let h = (b - a) / (num_points as f64 - 1.0);
        let mut l2_sum: f64 = 0.0;
        let mut l_inf: f64 = 0.0;
        let mut h1_sum: f64 = 0.0;
        let mut exact_l2_sum: f64 = 0.0;
        let mut exact_l_inf: f64 = 0.0;

        for i in 0..num_points {
            let x = a + (i as f64) * h;
            let u_approx = approx(x);
            let u_exact = exact(x);
            let error = (u_approx - u_exact).abs();

            l2_sum += error * error;
            l_inf = l_inf.max(error);
            exact_l2_sum += u_exact * u_exact;
            exact_l_inf = exact_l_inf.max(u_exact.abs());

            // Derivative error
            if i < num_points - 1 {
                let x_next = a + ((i + 1) as f64) * h;
                let du_approx = (approx(x_next) - u_approx) / h;
                let du_exact = (exact(x_next) - u_exact) / h;
                let deriv_error = (du_approx - du_exact).abs();
                h1_sum += deriv_error * deriv_error;
            }
        }

        let l2_error = (h * l2_sum).sqrt();
        let exact_l2_norm = (h * exact_l2_sum).sqrt();

        Self {
            l2_error,
            l_inf_error: l_inf,
            h1_seminorm_error: (h * h1_sum).sqrt(),
            relative_l2_error: if exact_l2_norm > 1e-12 {
                l2_error / exact_l2_norm
            } else {
                l2_error
            },
            relative_l_inf_error: if exact_l_inf > 1e-12 {
                l_inf / exact_l_inf
            } else {
                l_inf
            },
            num_eval_points: num_points,
            mesh_size: h,
        }
    }

    /// Compute convergence rate between two error measurements
    pub fn convergence_rate(&self, coarser: &ApproximationError) -> f64 {
        let h_ratio = coarser.mesh_size / self.mesh_size;
        (coarser.l2_error / self.l2_error).log2() / h_ratio.log2()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exponential() {
        let exp = Exponential::standard();
        assert_eq!(exp.name(), "Exponential e^x");

        // e^0 = 1
        assert!((exp.eval(0.0) - 1.0).abs() < 1e-10);

        // e^1 ≈ 2.71828
        assert!((exp.eval(1.0) - E).abs() < 1e-10);
    }

    #[test]
    fn test_sine() {
        let sine = Sine::standard();

        // sin(0) = 0
        assert!(sine.eval(0.0).abs() < 1e-10);

        // sin(π/2) = 1
        assert!((sine.eval(PI / 2.0) - 1.0).abs() < 1e-10);

        // sin(π) = 0
        assert!(sine.eval(PI).abs() < 1e-10);
    }

    #[test]
    fn test_error_function() {
        let erf = ErrorFunction::standard();

        // erf(0) = 0
        assert!(erf.eval(0.0).abs() < 1e-10);

        // erf(-x) = -erf(x)
        let x = 1.5;
        assert!((erf.eval(-x) + erf.eval(x)).abs() < 1e-6);

        // erf(∞) → 1
        assert!((erf.eval(3.0) - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_bessel_j0() {
        let j0 = BesselJ0::standard();

        // J_0(0) = 1
        assert!((j0.eval(0.0) - 1.0).abs() < 1e-6);

        // First zero of J_0 is at x ≈ 2.4048
        let first_zero = 2.4048;
        assert!(j0.eval(first_zero).abs() < 0.01);
    }

    #[test]
    fn test_runge_function() {
        let runge = RungsFunction::standard();

        // f(0) = 1
        assert!((runge.eval(0.0) - 1.0).abs() < 1e-10);

        // Symmetric: f(-x) = f(x)
        let x = 0.5;
        assert!((runge.eval(-x) - runge.eval(x)).abs() < 1e-10);
    }

    #[test]
    fn test_approximation_error() {
        let exact = |x: f64| x * x;
        let approx = |x: f64| x * x + 0.001 * x;

        let error = ApproximationError::compute(approx, exact, 0.0, 1.0, 101);

        assert!(error.l2_error > 0.0);
        assert!(error.l_inf_error > 0.0);
        assert!(error.l2_error < 0.01);
    }

    #[test]
    fn test_convergence_rate() {
        // For polynomial of degree n, expect convergence rate ~ n+1
        let exact = |x: f64| x * x * x;

        let error_h = ApproximationError::compute(
            |x| x * x * x + 0.1 * x,
            exact,
            0.0,
            1.0,
            51,
        );

        let error_2h = ApproximationError::compute(
            |x| x * x * x + 0.2 * x,
            exact,
            0.0,
            1.0,
            26,
        );

        // Just check that computation doesn't panic
        let _rate = error_h.convergence_rate(&error_2h);
    }
}
