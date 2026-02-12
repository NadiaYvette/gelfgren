//! Numerical methods
//!
//! This module provides numerical integration, differentiation,
//! ODE solvers, root finding, and optimization algorithms.

use crate::Result;

/// Integrate a function using the trapezoidal rule
///
/// # Arguments
/// * `f` - Function to integrate
/// * `a` - Lower bound
/// * `b` - Upper bound
/// * `n` - Number of trapezoids
///
/// # Returns
/// Approximate integral value
pub fn integrate_trapezoid<F>(f: F, a: f64, b: f64, n: usize) -> Result<f64>
where
    F: Fn(f64) -> f64,
{
    if n == 0 {
        return Err(crate::GelfgrenError::InvalidArgument(
            "Number of trapezoids must be positive".to_string()
        ));
    }

    let h = (b - a) / n as f64;
    let mut sum = 0.5 * (f(a) + f(b));

    for i in 1..n {
        let x = a + i as f64 * h;
        sum += f(x);
    }

    Ok(sum * h)
}

/// Find a root using the bisection method
///
/// # Arguments
/// * `f` - Function to find root of
/// * `a` - Lower bound (f(a) and f(b) must have opposite signs)
/// * `b` - Upper bound
/// * `tol` - Tolerance for convergence
/// * `max_iter` - Maximum number of iterations
///
/// # Returns
/// Approximate root
pub fn bisection<F>(f: F, mut a: f64, mut b: f64, tol: f64, max_iter: usize) -> Result<f64>
where
    F: Fn(f64) -> f64,
{
    let mut fa = f(a);
    let fb = f(b);

    if fa * fb > 0.0 {
        return Err(crate::GelfgrenError::InvalidArgument(
            "f(a) and f(b) must have opposite signs".to_string()
        ));
    }

    for _ in 0..max_iter {
        let c = (a + b) / 2.0;
        let fc = f(c);

        if fc.abs() < tol || (b - a).abs() < tol {
            return Ok(c);
        }

        if fa * fc < 0.0 {
            b = c;
            // Note: fb = fc omitted as it's not used in subsequent iterations
        } else {
            a = c;
            fa = fc;
        }
    }

    Err(crate::GelfgrenError::ConvergenceFailure(
        format!("Bisection failed to converge after {} iterations", max_iter)
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_trapezoid_integration() {
        // Integrate x^2 from 0 to 1, should be 1/3
        let f = |x: f64| x * x;
        let result = integrate_trapezoid(f, 0.0, 1.0, 1000).unwrap();
        assert_relative_eq!(result, 1.0 / 3.0, epsilon = 0.001);
    }

    #[test]
    fn test_bisection() {
        // Find root of x^2 - 2 (should be sqrt(2) ≈ 1.414)
        let f = |x: f64| x * x - 2.0;
        let result = bisection(f, 0.0, 2.0, 1e-6, 100).unwrap();
        assert_relative_eq!(result, 2.0_f64.sqrt(), epsilon = 1e-6);
    }
}
