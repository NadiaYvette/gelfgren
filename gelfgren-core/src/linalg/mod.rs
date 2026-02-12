//! Linear algebra operations
//!
//! This module provides matrix and vector operations, decompositions,
//! and linear system solvers.

mod vector;
mod matrix;

pub use vector::Vector;
pub use matrix::Matrix;

use crate::Result;

/// Compute the dot product of two vectors
///
/// # Arguments
/// * `a` - First vector
/// * `b` - Second vector
///
/// # Returns
/// The dot product as a scalar
///
/// # Errors
/// Returns `InvalidDimension` if vectors have different lengths
pub fn dot(a: &[f64], b: &[f64]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(crate::GelfgrenError::InvalidDimension(
            format!("Vector lengths must match: {} != {}", a.len(), b.len())
        ));
    }

    Ok(a.iter().zip(b.iter()).map(|(x, y)| x * y).sum())
}

/// Compute the Euclidean norm (L2 norm) of a vector
///
/// # Arguments
/// * `v` - Input vector
///
/// # Returns
/// The L2 norm
pub fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_dot_product() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];
        let result = dot(&a, &b).unwrap();
        assert_relative_eq!(result, 32.0);
    }

    #[test]
    fn test_dot_product_dimension_mismatch() {
        let a = vec![1.0, 2.0];
        let b = vec![1.0, 2.0, 3.0];
        assert!(dot(&a, &b).is_err());
    }

    #[test]
    fn test_norm() {
        let v = vec![3.0, 4.0];
        let result = norm(&v);
        assert_relative_eq!(result, 5.0);
    }
}
