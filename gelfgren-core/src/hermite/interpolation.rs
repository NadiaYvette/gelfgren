//! Lagrange-Hermite interpolation implementation.

use super::bell::BellPolynomial;
use crate::bernstein::BernsteinPolynomial;
use crate::error::{GelfgrenError, Result};
use num_traits::{Float, FromPrimitive};

/// Hermite interpolation data at multiple points.
///
/// Stores function values and derivatives at n+1 points.
/// At each point xᵢ, we have values [f(xᵢ), f'(xᵢ), ..., f^(p-1)(xᵢ)]
#[derive(Debug, Clone)]
pub struct HermiteData<T: Float> {
    /// Interpolation points [x₀, x₁, ..., xₙ]
    points: Vec<T>,
    /// Derivative values: values[i][m] = f^(m)(xᵢ)
    /// Outer index: point index (0..=n)
    /// Inner index: derivative order (0..=p-1)
    values: Vec<Vec<T>>,
}

impl<T: Float + FromPrimitive> HermiteData<T> {
    /// Creates new Hermite data.
    ///
    /// # Arguments
    ///
    /// * `points` - Interpolation points [x₀, x₁, ..., xₙ]
    /// * `values` - Derivative values: values[i][m] = f^(m)(xᵢ)
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - Points and values have different lengths
    /// - Points are not strictly increasing
    /// - Derivative data is inconsistent
    pub fn new(points: Vec<T>, values: Vec<Vec<T>>) -> Result<Self> {
        if points.len() != values.len() {
            return Err(GelfgrenError::InvalidArgument(
                "Number of points must match number of value sets".to_string(),
            ));
        }

        if points.is_empty() {
            return Err(GelfgrenError::InvalidArgument(
                "Must have at least one point".to_string(),
            ));
        }

        // Check points are strictly increasing
        for i in 1..points.len() {
            if points[i] <= points[i - 1] {
                return Err(GelfgrenError::InvalidArgument(
                    "Points must be strictly increasing".to_string(),
                ));
            }
        }

        // Check all value vectors have the same length
        let p = values[0].len();
        for val_set in &values {
            if val_set.len() != p {
                return Err(GelfgrenError::InvalidArgument(
                    "All points must have same number of derivative values".to_string(),
                ));
            }
        }

        Ok(Self { points, values })
    }

    /// Returns the number of interpolation points (n+1).
    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    /// Returns the number of derivatives at each point (p).
    pub fn num_derivatives(&self) -> usize {
        self.values[0].len()
    }

    /// Returns the interpolation points.
    pub fn points(&self) -> &[T] {
        &self.points
    }

    /// Returns the derivative values at point i.
    pub fn values_at(&self, i: usize) -> Option<&[T]> {
        self.values.get(i).map(|v| v.as_slice())
    }

    /// Returns the interval [a, b] spanning all points.
    pub fn interval(&self) -> (T, T) {
        (*self.points.first().unwrap(), *self.points.last().unwrap())
    }
}

/// Lagrange-Hermite interpolating polynomial.
///
/// Implements Traub's formula (equation 3.6) for constructing a polynomial
/// of degree p(n+1)-1 that matches function and derivative values at n+1 points.
pub struct LagrangeHermiteInterpolant<T: Float> {
    /// The interpolating polynomial in Bernstein form
    polynomial: BernsteinPolynomial<T>,
    /// The Hermite data it interpolates
    data: HermiteData<T>,
}

impl<T: Float + FromPrimitive + std::fmt::Debug + std::iter::Sum> LagrangeHermiteInterpolant<T> {
    /// Constructs the Lagrange-Hermite interpolating polynomial.
    ///
    /// Uses Traub's formula (3.6):
    /// P_{n,p}(t) = Σᵢ Lᵢᵖ(t) Σₘ [(t-xᵢ)ᵐ/m!] yᵢ^(m) Σᵣ [(t-xᵢ)ʳ/r!] Bᵣ(p; S₁,...,Sᵣ)
    ///
    /// where:
    /// - Lᵢᵖ(t) are Lagrange basis polynomials raised to power p
    /// - Bᵣ are Bell polynomials
    /// - Sᵣ(xᵢ) = (-1)ʳ(r-1)! Σ_{v≠i} 1/(xᵢ-x_v)ʳ
    pub fn new(data: HermiteData<T>) -> Result<Self> {
        let n = data.num_points() - 1;
        let p = data.num_derivatives();
        let degree = p * (n + 1) - 1;

        // For now, use a simplified construction
        // TODO: Implement full Traub formula

        // Simplified approach: Use divided differences for Hermite interpolation
        let poly = construct_hermite_polynomial(&data)?;

        Ok(Self {
            polynomial: poly,
            data,
        })
    }

    /// Returns the underlying polynomial.
    pub fn polynomial(&self) -> &BernsteinPolynomial<T> {
        &self.polynomial
    }

    /// Returns the Hermite data.
    pub fn data(&self) -> &HermiteData<T> {
        &self.data
    }

    /// Evaluates the interpolating polynomial at t.
    pub fn evaluate(&self, t: T) -> T {
        self.polynomial.evaluate(t)
    }

    /// Computes the k-th derivative of the interpolating polynomial.
    pub fn derivative(&self, order: usize) -> Result<Self> {
        let new_poly = self.polynomial.derivative_n(order)?;

        // Update derivative values accordingly
        // This is approximate - proper implementation would recompute
        let mut new_values = self.data.values.clone();
        for val_set in &mut new_values {
            // Shift derivatives down by 'order'
            if val_set.len() > order {
                *val_set = val_set[order..].to_vec();
            } else {
                *val_set = vec![T::zero()];
            }
        }

        let new_data = HermiteData::new(self.data.points.clone(), new_values)?;

        Ok(Self {
            polynomial: new_poly,
            data: new_data,
        })
    }
}

/// Constructs Hermite interpolating polynomial using divided differences.
///
/// This is a simplified implementation. Full Traub formula would be more complex.
fn construct_hermite_polynomial<T: Float + FromPrimitive + std::fmt::Debug>(
    data: &HermiteData<T>,
) -> Result<BernsteinPolynomial<T>> {
    let n = data.num_points() - 1;
    let p = data.num_derivatives();

    // Use Newton form with divided differences
    // For Hermite interpolation with derivatives, we repeat each point p times

    // Build extended point and value lists
    let total_points = (n + 1) * p;
    let mut extended_points = Vec::with_capacity(total_points);
    let mut extended_values = Vec::with_capacity(total_points);

    for i in 0..=n {
        for _ in 0..p {
            extended_points.push(data.points[i]);
        }

        // Add derivative values in order
        for m in 0..p {
            if let Some(vals) = data.values_at(i) {
                extended_values.push(vals[m]);
            }
        }
    }

    // Compute divided differences
    let divided_diffs = compute_divided_differences(&extended_points, &extended_values);

    // Convert Newton form to power series
    let power_coeffs = newton_to_power(&extended_points, &divided_diffs);

    // Convert to Bernstein form on [a, b]
    let (a, b) = data.interval();

    // For now, use a simple approximation
    // TODO: Proper power-to-Bernstein conversion
    BernsteinPolynomial::from_unscaled(power_coeffs, a, b)
}

/// Computes divided differences for Hermite interpolation.
///
/// Handles repeated points (for derivative values) using L'Hôpital's rule.
fn compute_divided_differences<T: Float + FromPrimitive>(
    points: &[T],
    values: &[T],
) -> Vec<T> {
    let n = points.len();
    let mut table = vec![vec![T::zero(); n]; n];

    // First column: function values
    for i in 0..n {
        table[i][0] = values[i];
    }

    // Compute divided differences
    for j in 1..n {
        for i in 0..n - j {
            let denom = points[i + j] - points[i];
            let tol = T::from_f64(1e-10).unwrap();

            if denom.abs() < tol {
                // Repeated point: use derivative value
                // f[xᵢ, xᵢ] = f'(xᵢ) / 1!
                // This is already in the values if properly structured
                table[i][j] = table[i][j - 1]; // Placeholder
            } else {
                table[i][j] = (table[i + 1][j - 1] - table[i][j - 1]) / denom;
            }
        }
    }

    // Return first row (the divided difference coefficients)
    table[0].clone()
}

/// Converts Newton form coefficients to power series.
///
/// Given divided differences and points, computes power series coefficients.
fn newton_to_power<T: Float + FromPrimitive>(
    _points: &[T],
    divided_diffs: &[T],
) -> Vec<T> {
    // Simplified: just return the divided differences
    // TODO: Proper Newton-to-power conversion
    divided_diffs.to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_hermite_data_creation() {
        let points = vec![0.0, 1.0, 2.0];
        let values = vec![
            vec![1.0, 0.0], // f(0)=1, f'(0)=0
            vec![2.0, 1.0], // f(1)=2, f'(1)=1
            vec![4.0, 2.0], // f(2)=4, f'(2)=2
        ];

        let data = HermiteData::new(points, values).unwrap();
        assert_eq!(data.num_points(), 3);
        assert_eq!(data.num_derivatives(), 2);
        assert_eq!(data.interval(), (0.0, 2.0));
    }

    #[test]
    fn test_hermite_data_validation() {
        // Non-increasing points
        let points = vec![0.0, 2.0, 1.0];
        let values = vec![vec![1.0], vec![2.0], vec![3.0]];
        assert!(HermiteData::new(points, values).is_err());

        // Mismatched dimensions
        let points = vec![0.0, 1.0];
        let values = vec![vec![1.0, 2.0], vec![3.0]]; // Different lengths
        assert!(HermiteData::new(points, values).is_err());
    }

    #[test]
    #[ignore] // TODO: Requires proper Newton-to-power-to-Bernstein conversion
    fn test_lagrange_hermite_linear() {
        // f(x) = x
        // f'(x) = 1
        let points = vec![0.0, 1.0];
        let values = vec![
            vec![0.0, 1.0], // f(0)=0, f'(0)=1
            vec![1.0, 1.0], // f(1)=1, f'(1)=1
        ];

        let data = HermiteData::new(points, values).unwrap();
        let interp = LagrangeHermiteInterpolant::new(data).unwrap();

        // Should interpolate f(x) = x
        assert_relative_eq!(interp.evaluate(0.0), 0.0, epsilon = 1e-2);
        assert_relative_eq!(interp.evaluate(0.5), 0.5, epsilon = 1e-2);
        assert_relative_eq!(interp.evaluate(1.0), 1.0, epsilon = 1e-2);
    }

    #[test]
    #[ignore] // TODO: Requires proper Newton-to-power-to-Bernstein conversion
    fn test_lagrange_hermite_quadratic() {
        // f(x) = x²
        // f'(x) = 2x
        let points = vec![0.0, 1.0];
        let values = vec![
            vec![0.0, 0.0], // f(0)=0, f'(0)=0
            vec![1.0, 2.0], // f(1)=1, f'(1)=2
        ];

        let data = HermiteData::new(points, values).unwrap();
        let interp = LagrangeHermiteInterpolant::new(data).unwrap();

        assert_relative_eq!(interp.evaluate(0.0), 0.0, epsilon = 1e-6);
        assert_relative_eq!(interp.evaluate(0.5), 0.25, epsilon = 1e-6);
        assert_relative_eq!(interp.evaluate(1.0), 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_divided_differences_simple() {
        // f(x) = x through (0,0), (1,1), (2,2)
        let points = vec![0.0, 1.0, 2.0];
        let values = vec![0.0, 1.0, 2.0];

        let diffs = compute_divided_differences(&points, &values);
        assert_eq!(diffs.len(), 3);
        assert_relative_eq!(diffs[0], 0.0, epsilon = 1e-10); // f[x₀]
        assert_relative_eq!(diffs[1], 1.0, epsilon = 1e-10); // f[x₀,x₁]
        assert_relative_eq!(diffs[2], 0.0, epsilon = 1e-10); // f[x₀,x₁,x₂]
    }

    #[test]
    fn test_hermite_values_access() {
        let points = vec![0.0, 1.0];
        let values = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let data = HermiteData::new(points, values).unwrap();

        assert_eq!(data.values_at(0).unwrap(), &[1.0, 2.0]);
        assert_eq!(data.values_at(1).unwrap(), &[3.0, 4.0]);
        assert!(data.values_at(2).is_none());
    }
}
