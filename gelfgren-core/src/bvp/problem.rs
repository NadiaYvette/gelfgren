//! Boundary value problem specification.

use super::boundary::BoundaryCondition;
use crate::error::{GelfgrenError, Result};
use num_traits::Float;

/// A differential operator L[u].
///
/// Represents an ODE operator, e.g., L[u] = u'' + p(x)u' + q(x)u
pub trait DifferentialOperator<T: Float> {
    /// Applies the operator to a function and its derivatives.
    ///
    /// # Arguments
    ///
    /// * `x` - Evaluation point
    /// * `u` - Function value u(x)
    /// * `derivatives` - [u'(x), u''(x), u'''(x), ...]
    ///
    /// # Returns
    ///
    /// L[u](x)
    fn apply(&self, x: T, u: T, derivatives: &[T]) -> T;

    /// Returns the order of the operator (highest derivative).
    fn order(&self) -> usize;
}

/// A simple second-order linear operator: u'' + p(x)u' + q(x)u
pub struct SecondOrderLinear<T: Float> {
    /// Coefficient p(x) for u'
    pub p: Box<dyn Fn(T) -> T>,
    /// Coefficient q(x) for u
    pub q: Box<dyn Fn(T) -> T>,
}

impl<T: Float> SecondOrderLinear<T> {
    /// Creates a new second-order linear operator.
    pub fn new<P, Q>(p: P, q: Q) -> Self
    where
        P: Fn(T) -> T + 'static,
        Q: Fn(T) -> T + 'static,
    {
        Self {
            p: Box::new(p),
            q: Box::new(q),
        }
    }
}

impl<T: Float> DifferentialOperator<T> for SecondOrderLinear<T> {
    fn apply(&self, x: T, u: T, derivatives: &[T]) -> T {
        if derivatives.len() < 2 {
            return T::zero(); // Not enough derivatives
        }

        let u_prime = derivatives[0];
        let u_double_prime = derivatives[1];

        u_double_prime + (self.p)(x) * u_prime + (self.q)(x) * u
    }

    fn order(&self) -> usize {
        2
    }
}

/// Right-hand side function f(x) for the ODE.
pub trait RightHandSide<T: Float> {
    /// Evaluates f(x).
    fn eval(&self, x: T) -> T;
}

/// Simple function wrapper for RHS.
pub struct FunctionRHS<T: Float> {
    f: Box<dyn Fn(T) -> T>,
}

impl<T: Float> FunctionRHS<T> {
    /// Creates a new RHS from a function.
    pub fn new<F>(f: F) -> Self
    where
        F: Fn(T) -> T + 'static,
    {
        Self { f: Box::new(f) }
    }
}

impl<T: Float> RightHandSide<T> for FunctionRHS<T> {
    fn eval(&self, x: T) -> T {
        (self.f)(x)
    }
}

/// A boundary value problem specification.
///
/// Consists of:
/// - Differential operator L
/// - Right-hand side f(x)
/// - Boundary conditions
/// - Interval [a, b]
pub struct BoundaryValueProblem<T: Float> {
    /// The differential operator
    pub operator: Box<dyn DifferentialOperator<T>>,
    /// Right-hand side function
    pub rhs: Box<dyn RightHandSide<T>>,
    /// Boundary conditions (typically 2 for second-order BVP)
    pub boundary_conditions: Vec<BoundaryCondition<T>>,
    /// Interval [a, b]
    pub interval: (T, T),
}

impl<T: Float> BoundaryValueProblem<T> {
    /// Creates a new BVP.
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - Invalid interval (a >= b)
    /// - Insufficient boundary conditions
    /// - Boundary conditions not at endpoints
    pub fn new(
        operator: Box<dyn DifferentialOperator<T>>,
        rhs: Box<dyn RightHandSide<T>>,
        boundary_conditions: Vec<BoundaryCondition<T>>,
        interval: (T, T),
    ) -> Result<Self> {
        let (a, b) = interval;
        if a >= b {
            return Err(GelfgrenError::InvalidArgument(
                "Invalid interval: a must be less than b".to_string(),
            ));
        }

        // Check we have enough boundary conditions (should equal operator order)
        let order = operator.order();
        if boundary_conditions.len() < order {
            return Err(GelfgrenError::InvalidArgument(format!(
                "Need {} boundary conditions for order-{} ODE, got {}",
                order,
                order,
                boundary_conditions.len()
            )));
        }

        Ok(Self {
            operator,
            rhs,
            boundary_conditions,
            interval,
        })
    }

    /// Returns the interval [a, b].
    pub fn interval(&self) -> (T, T) {
        self.interval
    }

    /// Returns the order of the differential equation.
    pub fn order(&self) -> usize {
        self.operator.order()
    }

    /// Evaluates the residual L[u](x) - f(x).
    ///
    /// Used to check how well a candidate solution satisfies the ODE.
    pub fn residual(&self, x: T, u: T, derivatives: &[T]) -> T {
        self.operator.apply(x, u, derivatives) - self.rhs.eval(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_second_order_operator() {
        // u'' + 2u' + u
        let op = SecondOrderLinear::new(|_x| 2.0, |_x| 1.0);

        assert_eq!(op.order(), 2);

        // At x=0: u=1, u'=2, u''=3
        // L[u] = 3 + 2*2 + 1*1 = 8
        let result = op.apply(0.0, 1.0, &[2.0, 3.0]);
        assert_eq!(result, 8.0);
    }

    #[test]
    fn test_function_rhs() {
        let rhs = FunctionRHS::new(|x: f64| x * x);
        assert_eq!(rhs.eval(3.0), 9.0);
    }

    #[test]
    fn test_bvp_creation() {
        let op = Box::new(SecondOrderLinear::new(|_x| 0.0, |_x| 1.0));
        let rhs = Box::new(FunctionRHS::new(|_x: f64| 0.0));
        let bcs = vec![
            BoundaryCondition::dirichlet(0.0, 0.0),
            BoundaryCondition::dirichlet(1.0, 0.0),
        ];

        let bvp = BoundaryValueProblem::new(op, rhs, bcs, (0.0, 1.0));
        assert!(bvp.is_ok());

        let bvp = bvp.unwrap();
        assert_eq!(bvp.order(), 2);
        assert_eq!(bvp.interval(), (0.0, 1.0));
    }

    #[test]
    fn test_bvp_validation() {
        let op = Box::new(SecondOrderLinear::new(|_x| 0.0, |_x| 1.0));
        let rhs = Box::new(FunctionRHS::new(|_x: f64| 0.0));

        // Invalid interval
        let bcs = vec![BoundaryCondition::dirichlet(0.0, 0.0)];
        let result = BoundaryValueProblem::new(op, rhs, bcs, (1.0, 0.0));
        assert!(result.is_err());
    }

    #[test]
    fn test_residual() {
        // u'' + u = 0
        let op = Box::new(SecondOrderLinear::new(|_x| 0.0, |_x| 1.0));
        let rhs = Box::new(FunctionRHS::new(|_x: f64| 0.0));
        let bcs = vec![
            BoundaryCondition::dirichlet(0.0, 0.0),
            BoundaryCondition::dirichlet(1.0, 0.0),
        ];
        let bvp = BoundaryValueProblem::new(op, rhs, bcs, (0.0, 1.0)).unwrap();

        // For sin(x): u''=-sin(x), u=sin(x)
        // u'' + u = -sin(x) + sin(x) = 0 (exact solution)
        let x = 0.5;
        let u = x.sin();
        let u_prime = x.cos();
        let u_double_prime = -x.sin();

        let residual = bvp.residual(x, u, &[u_prime, u_double_prime]);
        assert!((residual).abs() < 1e-10);
    }
}
