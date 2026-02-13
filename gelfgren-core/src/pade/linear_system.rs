//! Linear system solver for Padé approximant construction.
//!
//! Implements Gaussian elimination with partial pivoting for solving
//! the linear systems that arise in rational approximant construction.

use crate::error::{GelfgrenError, Result};
use num_traits::{Float, FromPrimitive};

/// Solves the linear system Ax = b using Gaussian elimination with partial pivoting.
///
/// # Arguments
///
/// * `a` - Coefficient matrix (n×n), provided in row-major order as flattened vector
/// * `b` - Right-hand side vector (n elements)
/// * `n` - Size of the system
///
/// # Returns
///
/// Solution vector x such that Ax = b
///
/// # Errors
///
/// Returns error if the matrix is singular or nearly singular.
pub fn solve_linear_system<T: Float + FromPrimitive>(
    a: &[T],
    b: &[T],
    n: usize,
) -> Result<Vec<T>> {
    if a.len() != n * n {
        return Err(GelfgrenError::InvalidArgument(format!(
            "Matrix must be {}×{} (got {} elements)",
            n,
            n,
            a.len()
        )));
    }
    if b.len() != n {
        return Err(GelfgrenError::InvalidArgument(format!(
            "Right-hand side must have {} elements (got {})",
            n,
            b.len()
        )));
    }

    // Make mutable copies for Gaussian elimination
    let mut a_copy: Vec<T> = a.to_vec();
    let mut b_copy: Vec<T> = b.to_vec();

    let tol = T::from_f64(1e-14).unwrap();

    // Forward elimination with partial pivoting
    for k in 0..n {
        // Find pivot
        let mut max_val = a_copy[k * n + k].abs();
        let mut pivot_row = k;

        for i in (k + 1)..n {
            let val = a_copy[i * n + k].abs();
            if val > max_val {
                max_val = val;
                pivot_row = i;
            }
        }

        // Check for singularity
        if max_val < tol {
            return Err(GelfgrenError::SingularMatrix);
        }

        // Swap rows if needed
        if pivot_row != k {
            for j in 0..n {
                let temp = a_copy[k * n + j];
                a_copy[k * n + j] = a_copy[pivot_row * n + j];
                a_copy[pivot_row * n + j] = temp;
            }
            let temp = b_copy[k];
            b_copy[k] = b_copy[pivot_row];
            b_copy[pivot_row] = temp;
        }

        // Eliminate column below diagonal
        for i in (k + 1)..n {
            let factor = a_copy[i * n + k] / a_copy[k * n + k];
            for j in k..n {
                a_copy[i * n + j] = a_copy[i * n + j] - factor * a_copy[k * n + j];
            }
            b_copy[i] = b_copy[i] - factor * b_copy[k];
        }
    }

    // Back substitution
    let mut x = vec![T::zero(); n];
    for i in (0..n).rev() {
        let mut sum = T::zero();
        for j in (i + 1)..n {
            sum = sum + a_copy[i * n + j] * x[j];
        }
        x[i] = (b_copy[i] - sum) / a_copy[i * n + i];
    }

    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_solve_2x2_system() {
        // System: 2x + 3y = 8
        //         4x + 5y = 14
        // Solution: x = 1, y = 2
        let a = vec![2.0, 3.0, 4.0, 5.0];
        let b = vec![8.0, 14.0];

        let x = solve_linear_system(&a, &b, 2).unwrap();

        assert_relative_eq!(x[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_solve_3x3_system() {
        // System: x + 2y + 3z = 14
        //         2x + 5y + 6z = 30
        //         3x + 6y + 10z = 45
        // Solution: x = 1, y = 2, z = 3
        let a = vec![
            1.0, 2.0, 3.0,
            2.0, 5.0, 6.0,
            3.0, 6.0, 10.0
        ];
        let b = vec![14.0, 30.0, 45.0];

        let x = solve_linear_system(&a, &b, 3).unwrap();

        assert_relative_eq!(x[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 2.0, epsilon = 1e-10);
        assert_relative_eq!(x[2], 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_singular_matrix() {
        // Singular system: x + y = 1
        //                  2x + 2y = 3 (inconsistent)
        let a = vec![1.0, 1.0, 2.0, 2.0];
        let b = vec![1.0, 3.0];

        let result = solve_linear_system(&a, &b, 2);
        assert!(result.is_err());
    }
}
