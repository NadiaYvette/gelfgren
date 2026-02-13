//! R bindings for Gelfgren library using extendr.

use extendr_api::prelude::*;
use gelfgren_core::bernstein::BernsteinPolynomial as CoreBernstein;
use gelfgren_core::rational::RationalFunction as CoreRational;
use gelfgren_core::pade::PadeApproximant as CorePade;
use gelfgren_core::piecewise::Mesh as CoreMesh;

/// Bernstein polynomial class for R.
///
/// @description
/// Represents a polynomial using the numerically stable Bernstein basis.
///
/// @details
/// Bernstein polynomials are defined on an interval [a, b] and provide
/// numerically stable operations for polynomial arithmetic and calculus.
///
/// @examples
/// # Create polynomial P(x) = 1 + 2x + 3x²
/// p <- bernstein_polynomial(c(1, 2, 3), 0, 1)
///
/// # Evaluate at points
/// x <- seq(0, 1, by = 0.25)
/// y <- p$evaluate(x)
///
/// # Compute derivative
/// p_prime <- p$derivative()
///
/// # Polynomial arithmetic
/// q <- bernstein_polynomial(c(1, 1), 0, 1)
/// sum_poly <- p$add(q)
///
/// @export
#[derive(Debug, Clone)]
struct BernsteinPolynomial {
    inner: CoreBernstein<f64>,
}

#[extendr]
impl BernsteinPolynomial {
    /// Create a new Bernstein polynomial
    ///
    /// @param coeffs Numeric vector of unscaled coefficients
    /// @param a Left endpoint of interval
    /// @param b Right endpoint of interval
    /// @return BernsteinPolynomial object
    fn new(coeffs: Vec<f64>, a: f64, b: f64) -> Result<Self> {
        let inner = CoreBernstein::new(coeffs, a, b)
            .map_err(|e| format!("Failed to create polynomial: {:?}", e))?;
        Ok(Self { inner })
    }

    /// Evaluate polynomial at points
    ///
    /// @param x Numeric vector of evaluation points
    /// @return Numeric vector of P(x) values
    fn evaluate(&self, x: Vec<f64>) -> Vec<f64> {
        x.iter().map(|&xi| self.inner.evaluate(xi)).collect()
    }

    /// Get polynomial degree
    ///
    /// @return Integer degree
    fn degree(&self) -> i32 {
        self.inner.degree() as i32
    }

    /// Get interval [a, b]
    ///
    /// @return Numeric vector c(a, b)
    fn interval(&self) -> Vec<f64> {
        let (a, b) = self.inner.interval();
        vec![a, b]
    }

    /// Compute derivative polynomial
    ///
    /// @return BernsteinPolynomial representing P'(x)
    fn derivative(&self) -> Result<Self> {
        let inner = self.inner.derivative()
            .map_err(|e| format!("Failed to compute derivative: {:?}", e))?;
        Ok(Self { inner })
    }

    /// Compute integral (antiderivative) polynomial
    ///
    /// @return BernsteinPolynomial representing ∫P(x)dx
    fn integral(&self) -> Result<Self> {
        let inner = self.inner.integral()
            .map_err(|e| format!("Failed to compute integral: {:?}", e))?;
        Ok(Self { inner })
    }

    /// Elevate polynomial degree by 1
    ///
    /// @return BernsteinPolynomial with degree n+1
    fn elevate(&self) -> Self {
        let inner = self.inner.degree_elevate();
        Self { inner }
    }

    /// Add two polynomials
    ///
    /// @param other Another BernsteinPolynomial
    /// @return BernsteinPolynomial representing P + Q
    fn add(&self, other: &BernsteinPolynomial) -> Result<Self> {
        let inner = (self.inner.clone() + other.inner.clone())
            .map_err(|e| format!("Failed to add polynomials: {:?}", e))?;
        Ok(Self { inner })
    }

    /// Subtract two polynomials
    ///
    /// @param other Another BernsteinPolynomial
    /// @return BernsteinPolynomial representing P - Q
    fn subtract(&self, other: &BernsteinPolynomial) -> Result<Self> {
        let inner = (self.inner.clone() - other.inner.clone())
            .map_err(|e| format!("Failed to subtract polynomials: {:?}", e))?;
        Ok(Self { inner })
    }

    /// Multiply two polynomials
    ///
    /// @param other Another BernsteinPolynomial
    /// @return BernsteinPolynomial representing P * Q
    fn multiply(&self, other: &BernsteinPolynomial) -> Result<Self> {
        let inner = (self.inner.clone() * other.inner.clone())
            .map_err(|e| format!("Failed to multiply polynomials: {:?}", e))?;
        Ok(Self { inner })
    }

    /// Print method for R
    fn print(&self) {
        let (a, b) = self.inner.interval();
        rprintln!("BernsteinPolynomial(degree={}, interval=[{}, {}])",
                  self.inner.degree(), a, b);
    }
}

/// Rational function class for R.
///
/// @description
/// Represents a rational function P(x)/Q(x) as ratio of Bernstein polynomials.
///
/// @export
#[derive(Debug, Clone)]
struct RationalFunction {
    inner: CoreRational<f64>,
}

#[extendr]
impl RationalFunction {
    /// Create a new rational function
    ///
    /// @param numerator BernsteinPolynomial for P(x)
    /// @param denominator BernsteinPolynomial for Q(x)
    /// @return RationalFunction object
    fn new(numerator: &BernsteinPolynomial, denominator: &BernsteinPolynomial) -> Result<Self> {
        let inner = CoreRational::new(
            numerator.inner.clone(),
            denominator.inner.clone(),
        ).map_err(|e| format!("Failed to create rational function: {:?}", e))?;
        Ok(Self { inner })
    }

    /// Evaluate rational function at points
    ///
    /// @param x Numeric vector of evaluation points
    /// @return Numeric vector of R(x) values
    /// @details Throws error if denominator is zero at any point
    fn evaluate(&self, x: Vec<f64>) -> Result<Vec<f64>> {
        x.iter()
            .map(|&xi| self.inner.evaluate(xi)
                .map_err(|e| format!("Evaluation failed at x={}: {:?}", xi, e)))
            .collect()
    }

    /// Compute derivative using quotient rule
    ///
    /// @return RationalFunction representing R'(x)
    fn derivative(&self) -> Result<Self> {
        let inner = self.inner.derivative()
            .map_err(|e| format!("Failed to compute derivative: {:?}", e))?;
        Ok(Self { inner })
    }

    fn print(&self) {
        rprintln!("RationalFunction(P/Q)");
    }
}

/// Padé approximant class for R.
///
/// @export
#[derive(Debug, Clone)]
struct PadeApproximant {
    inner: CorePade<f64>,
}

#[extendr]
impl PadeApproximant {
    /// Create Padé approximant from power series
    ///
    /// @param coeffs Numeric vector of power series coefficients
    /// @param n Numerator degree
    /// @param m Denominator degree
    /// @param center Expansion center
    /// @param a Left endpoint
    /// @param b Right endpoint
    /// @return PadeApproximant object
    fn new(coeffs: Vec<f64>, n: i32, m: i32, center: f64, a: f64, b: f64) -> Result<Self> {
        let inner = CorePade::from_series(&coeffs, n as usize, m as usize, center, (a, b))
            .map_err(|e| format!("Failed to create Padé approximant: {:?}", e))?;
        Ok(Self { inner })
    }

    /// Evaluate Padé approximant at points
    ///
    /// @param x Numeric vector of evaluation points
    /// @return Numeric vector of approximant values
    fn evaluate(&self, x: Vec<f64>) -> Result<Vec<f64>> {
        x.iter()
            .map(|&xi| self.inner.evaluate(xi)
                .map_err(|e| format!("Evaluation failed: {:?}", e)))
            .collect()
    }

    fn print(&self) {
        rprintln!("PadeApproximant([n/m])");
    }
}

/// Mesh partition class for R.
///
/// @export
#[derive(Debug, Clone)]
struct Mesh {
    inner: CoreMesh<f64>,
}

#[extendr]
impl Mesh {
    /// Create uniform mesh
    ///
    /// @param a Left endpoint
    /// @param b Right endpoint
    /// @param n Number of subintervals
    /// @return Mesh object
    fn uniform(a: f64, b: f64, n: i32) -> Result<Self> {
        let inner = CoreMesh::uniform(a, b, n as usize)
            .map_err(|e| format!("Failed to create mesh: {:?}", e))?;
        Ok(Self { inner })
    }

    /// Create Chebyshev mesh with boundary clustering
    ///
    /// @param a Left endpoint
    /// @param b Right endpoint
    /// @param n Number of subintervals
    /// @return Mesh object
    fn chebyshev(a: f64, b: f64, n: i32) -> Result<Self> {
        let inner = CoreMesh::chebyshev(a, b, n as usize)
            .map_err(|e| format!("Failed to create mesh: {:?}", e))?;
        Ok(Self { inner })
    }

    /// Get mesh interval
    fn interval(&self) -> Vec<f64> {
        let (a, b) = self.inner.interval();
        vec![a, b]
    }

    /// Get number of mesh points
    fn num_points(&self) -> i32 {
        self.inner.num_points() as i32
    }

    fn print(&self) {
        let (a, b) = self.inner.interval();
        rprintln!("Mesh(interval=[{}, {}], points={})", a, b, self.inner.num_points());
    }
}

// Macro to generate exports
extendr_module! {
    mod gelfgren;
    impl BernsteinPolynomial;
    impl RationalFunction;
    impl PadeApproximant;
    impl Mesh;
}
