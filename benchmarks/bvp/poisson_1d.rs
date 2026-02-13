/// Benchmark: 1D Poisson Equation
///
/// Problem: -u''(x) = f(x) for x ∈ [0, 1]
///          u(0) = 0, u(1) = 0
///
/// Test Cases:
/// 1. Smooth forcing: f(x) = π² sin(πx), exact solution: u(x) = sin(πx)
/// 2. Non-smooth forcing: f(x) = -2 for x ∈ [0.25, 0.75], 0 otherwise
/// 3. Oscillatory forcing: f(x) = 100π² sin(10πx)

use std::f64::consts::PI;

/// Test case 1: Smooth forcing with known exact solution
pub struct SmoothPoisson {
    pub name: &'static str,
}

impl SmoothPoisson {
    pub fn new() -> Self {
        Self {
            name: "Poisson 1D - Smooth (sin)",
        }
    }

    /// Forcing function f(x) = π² sin(πx)
    pub fn forcing(&self, x: f64) -> f64 {
        PI * PI * (PI * x).sin()
    }

    /// Exact solution u(x) = sin(πx)
    pub fn exact_solution(&self, x: f64) -> f64 {
        (PI * x).sin()
    }

    /// Boundary conditions: u(0) = 0, u(1) = 0
    pub fn left_bc(&self) -> f64 {
        0.0
    }

    pub fn right_bc(&self) -> f64 {
        0.0
    }
}

/// Test case 2: Piecewise constant forcing (discontinuous)
pub struct DiscontinuousPoisson {
    pub name: &'static str,
}

impl DiscontinuousPoisson {
    pub fn new() -> Self {
        Self {
            name: "Poisson 1D - Discontinuous forcing",
        }
    }

    /// Forcing function: f(x) = -2 for x ∈ [0.25, 0.75], 0 otherwise
    pub fn forcing(&self, x: f64) -> f64 {
        if x >= 0.25 && x <= 0.75 {
            -2.0
        } else {
            0.0
        }
    }

    /// Exact solution (piecewise polynomial)
    /// Solves -u'' = f with f = -2 on [0.25, 0.75], 0 otherwise
    /// Subject to u(0) = 0, u(1) = 0, with C¹ continuity
    pub fn exact_solution(&self, x: f64) -> f64 {
        if x < 0.25 {
            -0.5 * x
        } else if x <= 0.75 {
            x * x - x + 0.0625
        } else {
            0.5 * x - 0.5
        }
    }

    pub fn left_bc(&self) -> f64 {
        0.0
    }

    pub fn right_bc(&self) -> f64 {
        0.0
    }
}

/// Test case 3: Highly oscillatory forcing
pub struct OscillatoryPoisson {
    pub name: &'static str,
    pub frequency: f64,
}

impl OscillatoryPoisson {
    pub fn new(frequency: f64) -> Self {
        Self {
            name: "Poisson 1D - Oscillatory",
            frequency,
        }
    }

    /// Forcing function f(x) = (ωπ)² sin(ωπx)
    pub fn forcing(&self, x: f64) -> f64 {
        let omega = self.frequency;
        (omega * PI) * (omega * PI) * (omega * PI * x).sin()
    }

    /// Exact solution u(x) = sin(ωπx)
    pub fn exact_solution(&self, x: f64) -> f64 {
        (self.frequency * PI * x).sin()
    }

    pub fn left_bc(&self) -> f64 {
        0.0
    }

    pub fn right_bc(&self) -> f64 {
        0.0
    }
}

/// Mesh configuration for convergence study
#[derive(Debug, Clone)]
pub struct MeshConfig {
    pub num_intervals: usize,
    pub polynomial_degree: usize,
    pub rational_numerator_degree: usize,
    pub rational_denominator_degree: usize,
}

impl MeshConfig {
    /// Create uniform mesh configurations for convergence study
    pub fn convergence_sequence(max_intervals: usize) -> Vec<Self> {
        let mut configs = Vec::new();
        let mut n = 4;
        while n <= max_intervals {
            configs.push(MeshConfig {
                num_intervals: n,
                polynomial_degree: 3, // Cubic splines
                rational_numerator_degree: 2,
                rational_denominator_degree: 2, // [2/2] rational
            });
            n *= 2;
        }
        configs
    }

    /// Degrees of freedom for polynomial spline
    pub fn polynomial_dof(&self) -> usize {
        // For C1 cubic splines: n intervals, degree 3
        // DOF = n + 3 (n-1 interior nodes + 2 boundary + 2 for C1 continuity)
        self.num_intervals + 3
    }

    /// Degrees of freedom for rational approximant
    pub fn rational_dof(&self) -> usize {
        // For piecewise rational [m/n] on each interval
        self.num_intervals * (self.rational_numerator_degree + self.rational_denominator_degree + 2)
    }
}

/// Error metrics for solution quality
#[derive(Debug, Clone)]
pub struct ErrorMetrics {
    pub l2_error: f64,
    pub l_inf_error: f64,
    pub h1_seminorm_error: f64,
    pub num_eval_points: usize,
    pub mesh_size: f64,
}

impl ErrorMetrics {
    /// Compute error metrics given approximate and exact solutions
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

        for i in 0..num_points {
            let x = a + (i as f64) * h;
            let u_approx = approx(x);
            let u_exact = exact(x);
            let error = (u_approx - u_exact).abs();

            l2_sum += error * error;
            l_inf = l_inf.max(error);

            // Compute derivative error (finite difference approximation)
            if i < num_points - 1 {
                let x_next = a + ((i + 1) as f64) * h;
                let du_approx = (approx(x_next) - u_approx) / h;
                let du_exact = (exact(x_next) - u_exact) / h;
                let deriv_error = (du_approx - du_exact).abs();
                h1_sum += deriv_error * deriv_error;
            }
        }

        Self {
            l2_error: (h * l2_sum).sqrt(),
            l_inf_error: l_inf,
            h1_seminorm_error: (h * h1_sum).sqrt(),
            num_eval_points: num_points,
            mesh_size: h,
        }
    }

    /// Compute convergence rate between two error measurements
    pub fn convergence_rate(&self, coarser: &ErrorMetrics) -> f64 {
        let h_ratio = coarser.mesh_size / self.mesh_size;
        (coarser.l2_error / self.l2_error).log2() / h_ratio.log2()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smooth_poisson_exact_solution() {
        let problem = SmoothPoisson::new();

        // Test that forcing and exact solution are consistent
        // -u'' = f => -(sin(πx))'' = -(-π² sin(πx)) = π² sin(πx) ✓

        for &x in &[0.1, 0.3, 0.5, 0.7, 0.9] {
            let f = problem.forcing(x);
            let u = problem.exact_solution(x);

            // Verify second derivative numerically
            let h = 1e-5;  // Step size for finite difference
            let u_plus = problem.exact_solution(x + h);
            let u_minus = problem.exact_solution(x - h);
            let u_double_prime = (u_plus - 2.0 * u + u_minus) / (h * h);

            let residual = (f + u_double_prime).abs();
            assert!(
                residual < 1e-3,  // Tolerance for finite difference approximation
                "PDE not satisfied at x={}: f={}, u''={}, residual={}",
                x, f, u_double_prime, residual
            );
        }
    }

    #[test]
    fn test_discontinuous_poisson_exact_solution() {
        let problem = DiscontinuousPoisson::new();

        // Test boundary conditions
        assert_eq!(problem.exact_solution(0.0), 0.0);
        assert_eq!(problem.exact_solution(1.0), 0.0);

        // Test continuity at x=0.25 and x=0.75
        let eps = 1e-10;
        assert!((problem.exact_solution(0.25 - eps) - problem.exact_solution(0.25 + eps)).abs() < 1e-6);
        assert!((problem.exact_solution(0.75 - eps) - problem.exact_solution(0.75 + eps)).abs() < 1e-6);
    }

    #[test]
    fn test_error_metrics() {
        let exact = |x: f64| x * x;
        let approx = |x: f64| x * x + 0.01 * x;

        let metrics = ErrorMetrics::compute(approx, exact, 0.0, 1.0, 101);

        assert!(metrics.l2_error > 0.0);
        assert!(metrics.l_inf_error > 0.0);
        assert!(metrics.l2_error < 0.1); // Should be small
    }

    #[test]
    fn test_mesh_config_dof() {
        let config = MeshConfig {
            num_intervals: 10,
            polynomial_degree: 3,
            rational_numerator_degree: 2,
            rational_denominator_degree: 2,
        };

        let poly_dof = config.polynomial_dof();
        let rat_dof = config.rational_dof();

        assert_eq!(poly_dof, 13); // 10 + 3
        assert_eq!(rat_dof, 60);  // 10 * (2 + 2 + 2)
    }
}
