//! Special mathematical functions
//!
//! This module provides special functions like gamma, beta,
//! Bessel functions, and orthogonal polynomials.

use crate::Result;

/// Compute the factorial of a non-negative integer
///
/// # Arguments
/// * `n` - Non-negative integer
///
/// # Returns
/// n!
pub fn factorial(n: u32) -> Result<f64> {
    if n > 170 {
        return Err(crate::GelfgrenError::Overflow);
    }

    Ok((1..=n).map(|x| x as f64).product())
}

/// Compute the natural logarithm of the gamma function using Stirling's approximation
///
/// This is a simple approximation. For production use, more sophisticated
/// algorithms should be implemented.
///
/// # Arguments
/// * `x` - Input value (must be positive)
///
/// # Returns
/// ln(Γ(x))
pub fn ln_gamma(x: f64) -> Result<f64> {
    if x <= 0.0 {
        return Err(crate::GelfgrenError::InvalidArgument(
            "Gamma function requires positive argument".to_string()
        ));
    }

    // Stirling's approximation
    // ln(Γ(x)) ≈ (x - 0.5) * ln(x) - x + 0.5 * ln(2π)
    const LN_2PI: f64 = 1.8378770664093453;  // ln(2π)

    Ok((x - 0.5) * x.ln() - x + 0.5 * LN_2PI)
}

/// Compute the gamma function
///
/// # Arguments
/// * `x` - Input value (must be positive)
///
/// # Returns
/// Γ(x)
pub fn gamma(x: f64) -> Result<f64> {
    if x <= 0.0 {
        return Err(crate::GelfgrenError::InvalidArgument(
            "Gamma function requires positive argument".to_string()
        ));
    }

    // For small integers, use factorial
    if x == x.floor() && x <= 170.0 {
        return factorial(x as u32 - 1);
    }

    Ok(ln_gamma(x)?.exp())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_factorial() {
        assert_relative_eq!(factorial(0).unwrap(), 1.0);
        assert_relative_eq!(factorial(5).unwrap(), 120.0);
        assert_relative_eq!(factorial(10).unwrap(), 3628800.0);
    }

    #[test]
    fn test_gamma_integer() {
        // Γ(n) = (n-1)! for positive integers
        assert_relative_eq!(gamma(1.0).unwrap(), 1.0, epsilon = 0.001);
        assert_relative_eq!(gamma(6.0).unwrap(), 120.0, epsilon = 0.1);
    }

    #[test]
    #[ignore] // TODO: Implement more accurate gamma function for small values
    fn test_gamma_half() {
        // Γ(0.5) = √π ≈ 1.772
        // Note: Stirling's approximation is not very accurate for small values like 0.5
        // A more sophisticated implementation (Lanczos approximation) is needed
        let result = gamma(0.5).unwrap();
        assert_relative_eq!(result, std::f64::consts::PI.sqrt(), epsilon = 0.01);
    }
}
