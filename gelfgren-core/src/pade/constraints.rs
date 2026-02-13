//! Constraint interfaces for Hermite interpolation with rational approximants.
//!
//! This module exposes the constraint equations used in two-point Padé construction,
//! enabling users to:
//! - Inspect the linear system before solving
//! - Use custom solvers (Newton, optimization, etc.)
//! - Handle overdetermined/underdetermined systems
//! - Transform problems to optimization form

use crate::bernstein::BernsteinPolynomial;
use crate::error::{GelfgrenError, Result};
use crate::rational::RationalFunction;
use num_traits::{Float, FromPrimitive};

/// Hermite interpolation constraints for rational approximants.
///
/// Represents the matching conditions R^(k)(xᵢ) = f^(k)(xᵢ) that define
/// a rational approximant R(x) = P(x)/Q(x) matching derivatives at endpoints.
///
/// # Mathematical Background
///
/// For a [n/m] rational approximant matching p derivatives at each of 2 endpoints:
/// - Total constraints: 2p
/// - Total unknowns: (n+1) + m = n + m + 1 (with normalization b₀ + bₙ = 1)
/// - Constraint: n + m + 1 = 2p for exact determination
///
/// # Examples
///
/// ```rust
/// use gelfgren_core::pade::HermiteConstraints;
///
/// // Set up constraints for e^x on [0,1] with [2/1] rational
/// let left = vec![1.0, 1.0];          // [f(0), f'(0)]
/// let right = vec![2.718, 2.718];     // [f(1), f'(1)]
///
/// let constraints = HermiteConstraints::new(
///     left, right, 2, 1, 0.0, 1.0
/// ).unwrap();
///
/// // Get linear system Ax = b
/// let (a, b, n) = constraints.linear_system();
///
/// // Or evaluate residuals for given coefficients
/// let coeffs = vec![/* ... */];
/// let residuals = constraints.residuals(&coeffs);
/// ```
#[derive(Debug, Clone)]
pub struct HermiteConstraints<T: Float> {
    /// Derivative values at left endpoint [f(a), f'(a), ...]
    left_derivatives: Vec<T>,
    /// Derivative values at right endpoint [f(b), f'(b), ...]
    right_derivatives: Vec<T>,
    /// Numerator degree
    n: usize,
    /// Denominator degree
    m: usize,
    /// Left endpoint
    x0: T,
    /// Right endpoint
    x1: T,
    /// Spacing Δx = x₁ - x₀
    delta_x: T,
    /// Order p (number of derivatives at each endpoint)
    p: usize,
}

impl<T: Float + FromPrimitive + std::fmt::Debug + std::iter::Sum> HermiteConstraints<T> {
    /// Creates new Hermite interpolation constraints.
    ///
    /// # Arguments
    ///
    /// * `left_derivatives` - [f(a), f'(a), ..., f^(p-1)(a)]
    /// * `right_derivatives` - [f(b), f'(b), ..., f^(p-1)(b)]
    /// * `n` - Numerator degree
    /// * `m` - Denominator degree
    /// * `x0` - Left endpoint a
    /// * `x1` - Right endpoint b
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - n + m + 1 ≠ 2p (must be even for two-point Padé)
    /// - Derivative vectors have different lengths
    /// - x₀ ≥ x₁
    pub fn new(
        left_derivatives: Vec<T>,
        right_derivatives: Vec<T>,
        n: usize,
        m: usize,
        x0: T,
        x1: T,
    ) -> Result<Self> {
        if left_derivatives.len() != right_derivatives.len() {
            return Err(GelfgrenError::InvalidArgument(
                "Left and right derivative vectors must have same length".to_string(),
            ));
        }

        let p = left_derivatives.len();

        if (n + m + 1) != 2 * p {
            return Err(GelfgrenError::InvalidArgument(format!(
                "Two-point Padé requires n+m+1 = 2p (even). Got n={}, m={}, p={}, n+m+1={}",
                n,
                m,
                p,
                n + m + 1
            )));
        }

        if x0 >= x1 {
            return Err(GelfgrenError::InvalidArgument(
                "Left endpoint must be less than right endpoint".to_string(),
            ));
        }

        let delta_x = x1 - x0;

        Ok(Self {
            left_derivatives,
            right_derivatives,
            n,
            m,
            delta_x,
            x0,
            x1,
            p,
        })
    }

    /// Returns the linear system Ax = b representing the matching conditions.
    ///
    /// The system has 2p equations (p at each endpoint) and n+m+1 unknowns
    /// (numerator and denominator coefficients).
    ///
    /// # Returns
    ///
    /// Tuple (matrix_a, vector_b, num_unknowns) where:
    /// - `matrix_a`: Flattened constraint matrix (row-major order)
    /// - `vector_b`: Right-hand side vector
    /// - `num_unknowns`: Number of variables (n+m+1)
    ///
    /// # Note
    ///
    /// This exposes the same linear system that `construct_rational_pade`
    /// solves internally. You can:
    /// - Inspect the system structure
    /// - Solve with your own method
    /// - Check condition number
    /// - Handle special cases
    pub fn linear_system(&self) -> (Vec<T>, Vec<T>, usize) {
        super::two_point::build_hermite_linear_system(
            &self.left_derivatives,
            &self.right_derivatives,
            self.n,
            self.m,
            self.x0,
            self.x1,
            self.delta_x,
        )
    }

    /// Number of constraint equations (2p).
    pub fn num_constraints(&self) -> usize {
        2 * self.p
    }

    /// Number of unknowns (n + m + 1).
    pub fn num_unknowns(&self) -> usize {
        self.n + self.m + 1
    }

    /// Returns the numerator degree.
    pub fn numerator_degree(&self) -> usize {
        self.n
    }

    /// Returns the denominator degree.
    pub fn denominator_degree(&self) -> usize {
        self.m
    }

    /// Returns the order (number of derivatives at each endpoint).
    pub fn order(&self) -> usize {
        self.p
    }

    /// Returns the interval [a, b].
    pub fn interval(&self) -> (T, T) {
        (self.x0, self.x1)
    }

    /// Evaluates constraint residuals given coefficient vector.
    ///
    /// # Arguments
    ///
    /// * `coeffs` - Coefficient vector [a₀, ..., aₙ, b₁, ..., bₘ]
    ///              where aᵢ are numerator coefficients (Bernstein basis)
    ///              and bⱼ are denominator coefficients with b₀ = 1 - Σbⱼ
    ///
    /// # Returns
    ///
    /// Vector of residuals r where rᵢ = R^(k)(xⱼ) - f^(k)(xⱼ)
    /// Ideally all residuals should be zero for a valid solution.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// let constraints = HermiteConstraints::new(/* ... */);
    /// let coeffs = vec![1.0, 2.0, 3.0, 0.5]; // [a₀, a₁, a₂, b₁] for [2/1]
    /// let residuals = constraints.residuals(&coeffs);
    ///
    /// // Use for root-finding
    /// let solution = newton_method(|x| constraints.residuals(x), x0);
    /// ```
    pub fn residuals(&self, coeffs: &[T]) -> Result<Vec<T>> {
        if coeffs.len() != self.num_unknowns() {
            return Err(GelfgrenError::InvalidArgument(format!(
                "Expected {} coefficients, got {}",
                self.num_unknowns(),
                coeffs.len()
            )));
        }

        // Extract numerator and denominator coefficients
        let num_coeffs = &coeffs[0..=self.n];
        let den_coeffs_partial = &coeffs[self.n + 1..];

        // Reconstruct full denominator with normalization
        let mut den_coeffs = vec![T::zero(); self.m + 1];
        let sum: T = den_coeffs_partial.iter().copied().sum();
        den_coeffs[0] = T::one() - sum;
        for (i, &val) in den_coeffs_partial.iter().enumerate() {
            den_coeffs[i + 1] = val;
        }

        // Build polynomials
        let numerator = BernsteinPolynomial::new(num_coeffs.to_vec(), self.x0, self.x1)?;
        let denominator = BernsteinPolynomial::new(den_coeffs, self.x0, self.x1)?;

        // Build rational
        let rational = RationalFunction::new(numerator, denominator)?;

        // Compute residuals
        let mut residuals = Vec::with_capacity(self.num_constraints());

        // Left endpoint constraints
        for k in 0..self.p {
            let r_deriv = self.evaluate_derivative(&rational, self.x0, k)?;
            let target = self.left_derivatives[k];
            residuals.push(r_deriv - target);
        }

        // Right endpoint constraints
        for k in 0..self.p {
            let r_deriv = self.evaluate_derivative(&rational, self.x1, k)?;
            let target = self.right_derivatives[k];
            residuals.push(r_deriv - target);
        }

        Ok(residuals)
    }

    /// Objective function: sum of squared residuals.
    ///
    /// Returns ||r||² where r is the residual vector.
    /// Useful for least-squares optimization.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// let constraints = HermiteConstraints::new(/* ... */);
    ///
    /// // Minimize objective to find best-fit rational
    /// let solution = minimize(|x| constraints.objective(x), x0);
    /// ```
    pub fn objective(&self, coeffs: &[T]) -> Result<T> {
        let residuals = self.residuals(coeffs)?;
        Ok(residuals.iter().map(|&r| r * r).sum())
    }

    /// Computes the Jacobian matrix of residuals with respect to coefficients.
    ///
    /// Returns J where Jᵢⱼ = ∂rᵢ/∂cⱼ
    ///
    /// # Returns
    ///
    /// Flattened Jacobian matrix in row-major order
    /// (num_constraints × num_unknowns matrix)
    ///
    /// # Note
    ///
    /// Uses finite differences with step size h = √ε where ε is machine epsilon.
    /// For more accuracy, implement analytical Jacobian (future work).
    pub fn jacobian(&self, coeffs: &[T]) -> Result<Vec<T>> {
        let num_constraints = self.num_constraints();
        let num_unknowns = self.num_unknowns();

        let epsilon = T::epsilon();
        let h = epsilon.sqrt(); // Step size for finite differences

        let mut jacobian = vec![T::zero(); num_constraints * num_unknowns];

        let base_residuals = self.residuals(coeffs)?;

        for j in 0..num_unknowns {
            // Perturb coefficient j
            let mut coeffs_perturbed = coeffs.to_vec();
            coeffs_perturbed[j] = coeffs_perturbed[j] + h;

            let perturbed_residuals = self.residuals(&coeffs_perturbed)?;

            // Compute finite difference: (r(c+h) - r(c)) / h
            for i in 0..num_constraints {
                let derivative = (perturbed_residuals[i] - base_residuals[i]) / h;
                jacobian[i * num_unknowns + j] = derivative;
            }
        }

        Ok(jacobian)
    }

    /// Returns quadratic form matrices for least-squares optimization.
    ///
    /// For the linear least squares problem minimize ||Ax - b||², this
    /// is equivalent to minimize x^T H x + g^T x where:
    /// - H = A^T A (Hessian matrix)
    /// - g = -2 A^T b (gradient vector)
    ///
    /// # Returns
    ///
    /// Tuple (hessian, gradient) where:
    /// - `hessian`: Flattened H matrix (num_unknowns × num_unknowns)
    /// - `gradient`: g vector (num_unknowns)
    ///
    /// # Examples
    ///
    /// ```ignore
    /// let (h, g) = constraints.quadratic_form();
    ///
    /// // Use with quadratic programming solver
    /// let solution = quadratic_minimize(h, g, bounds);
    /// ```
    pub fn quadratic_form(&self) -> (Vec<T>, Vec<T>) {
        let (a, b, n) = self.linear_system();
        let num_rows = self.num_constraints();
        let num_cols = n;

        // Compute H = A^T A
        let mut hessian = vec![T::zero(); num_cols * num_cols];
        for i in 0..num_cols {
            for j in 0..num_cols {
                let mut sum = T::zero();
                for k in 0..num_rows {
                    sum = sum + a[k * num_cols + i] * a[k * num_cols + j];
                }
                hessian[i * num_cols + j] = sum;
            }
        }

        // Compute g = -2 A^T b
        let mut gradient = vec![T::zero(); num_cols];
        for i in 0..num_cols {
            let mut sum = T::zero();
            for k in 0..num_rows {
                sum = sum + a[k * num_cols + i] * b[k];
            }
            gradient[i] = sum * T::from(-2.0).unwrap();
        }

        (hessian, gradient)
    }

    /// Builds a rational function from coefficient vector.
    ///
    /// # Arguments
    ///
    /// * `coeffs` - Coefficient vector [a₀, ..., aₙ, b₁, ..., bₘ]
    ///
    /// # Returns
    ///
    /// RationalFunction R(x) = P(x)/Q(x) with given coefficients
    pub fn build_rational(&self, coeffs: &[T]) -> Result<RationalFunction<T>> {
        if coeffs.len() != self.num_unknowns() {
            return Err(GelfgrenError::InvalidArgument(format!(
                "Expected {} coefficients, got {}",
                self.num_unknowns(),
                coeffs.len()
            )));
        }

        // Extract and build numerator
        let num_coeffs = &coeffs[0..=self.n];
        let numerator = BernsteinPolynomial::new(num_coeffs.to_vec(), self.x0, self.x1)?;

        // Reconstruct denominator with normalization
        let den_coeffs_partial = &coeffs[self.n + 1..];
        let mut den_coeffs = vec![T::zero(); self.m + 1];
        let sum: T = den_coeffs_partial.iter().copied().sum();
        den_coeffs[0] = T::one() - sum;
        for (i, &val) in den_coeffs_partial.iter().enumerate() {
            den_coeffs[i + 1] = val;
        }

        let denominator = BernsteinPolynomial::new(den_coeffs, self.x0, self.x1)?;

        RationalFunction::new(numerator, denominator)
    }

    /// Helper: Evaluate k-th derivative of rational function at x.
    fn evaluate_derivative(&self, rational: &RationalFunction<T>, x: T, k: usize) -> Result<T> {
        if k == 0 {
            rational.evaluate(x)
        } else {
            // Use finite differences for higher derivatives
            // TODO: Implement analytical derivatives
            let h = T::from_f64(1e-7).unwrap();
            match k {
                1 => {
                    let f_plus = rational.evaluate(x + h)?;
                    let f_minus = rational.evaluate(x - h)?;
                    Ok((f_plus - f_minus) / (T::from(2.0).unwrap() * h))
                }
                2 => {
                    let f_center = rational.evaluate(x)?;
                    let f_plus = rational.evaluate(x + h)?;
                    let f_minus = rational.evaluate(x - h)?;
                    Ok((f_plus - T::from(2.0).unwrap() * f_center + f_minus) / (h * h))
                }
                _ => Err(GelfgrenError::NotImplemented(format!(
                    "Derivative order {} not yet supported",
                    k
                ))),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constraint_interface_basic() {
        // Test e^x on [0,1] with [2/1] rational approximant
        // e^0 = 1, e^1 ≈ 2.718, derivatives are the same
        let left = vec![1.0, 1.0];
        let right = vec![2.718281828, 2.718281828];

        let constraints = HermiteConstraints::new(left, right, 2, 1, 0.0, 1.0).unwrap();

        // Check dimensions
        assert_eq!(constraints.num_constraints(), 4); // 2*2
        assert_eq!(constraints.num_unknowns(), 4); // 2+1+1
        assert_eq!(constraints.numerator_degree(), 2);
        assert_eq!(constraints.denominator_degree(), 1);
        assert_eq!(constraints.order(), 2);
        assert_eq!(constraints.interval(), (0.0, 1.0));

        // Test linear system extraction
        let (a, b, n) = constraints.linear_system();
        assert_eq!(n, 4); // 4 unknowns
        assert_eq!(a.len(), 16); // 4x4 matrix
        assert_eq!(b.len(), 4); // 4 equations

        // The matrix should not be all zeros
        let sum: f64 = a.iter().map(|x| x.abs()).sum();
        assert!(sum > 0.0);
    }

    #[test]
    fn test_residuals_and_objective() {
        let left = vec![1.0, 1.0];
        let right = vec![2.718281828, 2.718281828];
        let constraints = HermiteConstraints::new(left, right, 2, 1, 0.0, 1.0).unwrap();

        // Test with some example coefficients
        let coeffs = vec![1.0, 2.0, 1.5, 0.5]; // [p0, p1, p2, q1]

        // Test residuals
        let residuals = constraints.residuals(&coeffs).unwrap();
        assert_eq!(residuals.len(), 4);

        // Test objective (should be non-negative)
        let obj = constraints.objective(&coeffs).unwrap();
        assert!(obj >= 0.0);

        // Objective should be sum of squared residuals
        let expected_obj: f64 = residuals.iter().map(|r| r * r).sum();
        assert!((obj - expected_obj).abs() < 1e-10);
    }

    #[test]
    fn test_jacobian() {
        let left = vec![1.0, 1.0];
        let right = vec![2.718281828, 2.718281828];
        let constraints = HermiteConstraints::new(left, right, 2, 1, 0.0, 1.0).unwrap();

        let coeffs = vec![1.0, 2.0, 1.5, 0.5];

        // Test Jacobian computation
        let jac = constraints.jacobian(&coeffs).unwrap();
        assert_eq!(jac.len(), 16); // 4x4 matrix

        // Jacobian should not be all zeros
        let sum: f64 = jac.iter().map(|x| x.abs()).sum();
        assert!(sum > 0.0);
    }

    #[test]
    fn test_quadratic_form() {
        let left = vec![1.0, 1.0];
        let right = vec![2.718281828, 2.718281828];
        let constraints = HermiteConstraints::new(left, right, 2, 1, 0.0, 1.0).unwrap();

        // Test quadratic form extraction
        let (h, g) = constraints.quadratic_form();
        assert_eq!(h.len(), 16); // 4x4 Hessian
        assert_eq!(g.len(), 4); // gradient

        // Hessian should be symmetric (H = A^T A)
        for i in 0..4 {
            for j in 0..4 {
                let h_ij = h[i * 4 + j];
                let h_ji = h[j * 4 + i];
                assert!((h_ij - h_ji).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_build_rational() {
        let left = vec![1.0, 1.0];
        let right = vec![2.718281828, 2.718281828];
        let constraints = HermiteConstraints::new(left, right, 2, 1, 0.0, 1.0).unwrap();

        // Get solution by solving linear system
        let (a, b, n) = constraints.linear_system();

        // Simple solution (not optimal, just for testing)
        let coeffs = vec![1.0, 1.5, 1.0, 0.3];

        // Build rational from coefficients
        let rational = constraints.build_rational(&coeffs).unwrap();

        // Verify rational can be evaluated
        let val = rational.evaluate(0.5).unwrap();
        assert!(val.is_finite());
    }

    #[test]
    fn test_validation_different_lengths() {
        let left = vec![1.0, 1.0];
        let right = vec![2.718281828]; // Different length

        let result = HermiteConstraints::new(left, right, 2, 1, 0.0, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_validation_invalid_interval() {
        let left = vec![1.0, 1.0];
        let right = vec![2.718281828, 2.718281828];

        // x0 >= x1
        let result = HermiteConstraints::new(left, right, 2, 1, 1.0, 0.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_validation_wrong_dof() {
        let left = vec![1.0, 1.0]; // p=2, so 2p=4
        let right = vec![2.718281828, 2.718281828];

        // n+m+1 = 3+2+1 = 6 ≠ 4
        let result = HermiteConstraints::new(left, right, 3, 2, 0.0, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_residuals_wrong_size() {
        let left = vec![1.0, 1.0];
        let right = vec![2.718281828, 2.718281828];
        let constraints = HermiteConstraints::new(left, right, 2, 1, 0.0, 1.0).unwrap();

        // Wrong number of coefficients
        let coeffs = vec![1.0, 2.0]; // Should be 4

        let result = constraints.residuals(&coeffs);
        assert!(result.is_err());
    }
}
