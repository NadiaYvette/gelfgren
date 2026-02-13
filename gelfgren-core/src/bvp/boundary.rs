//! Boundary condition types and specifications.

use num_traits::Float;

/// Type of boundary condition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryConditionType {
    /// Dirichlet: u(x) = value
    Dirichlet,
    /// Neumann: u'(x) = value
    Neumann,
    /// Robin: c₀·u(x) + c₁·u'(x) = value
    Robin,
}

/// A boundary condition at one endpoint.
///
/// Specifies constraints on the solution and/or its derivatives.
#[derive(Debug, Clone)]
pub struct BoundaryCondition<T: Float> {
    /// Type of boundary condition
    pub bc_type: BoundaryConditionType,
    /// Location (typically a or b)
    pub location: T,
    /// Value (right-hand side)
    pub value: T,
    /// Coefficients for Robin BC: [c₀, c₁]
    /// For Dirichlet/Neumann, this is ignored
    pub coefficients: Vec<T>,
}

impl<T: Float> BoundaryCondition<T> {
    /// Creates a Dirichlet boundary condition: u(x) = value
    pub fn dirichlet(location: T, value: T) -> Self {
        Self {
            bc_type: BoundaryConditionType::Dirichlet,
            location,
            value,
            coefficients: Vec::new(),
        }
    }

    /// Creates a Neumann boundary condition: u'(x) = value
    pub fn neumann(location: T, value: T) -> Self {
        Self {
            bc_type: BoundaryConditionType::Neumann,
            location,
            value,
            coefficients: Vec::new(),
        }
    }

    /// Creates a Robin boundary condition: c₀·u(x) + c₁·u'(x) = value
    pub fn robin(location: T, c0: T, c1: T, value: T) -> Self {
        Self {
            bc_type: BoundaryConditionType::Robin,
            location,
            value,
            coefficients: vec![c0, c1],
        }
    }

    /// Checks if this boundary condition is at the left endpoint.
    pub fn is_left(&self, a: T) -> bool {
        (self.location - a).abs() < T::from(1e-10).unwrap()
    }

    /// Checks if this boundary condition is at the right endpoint.
    pub fn is_right(&self, b: T) -> bool {
        (self.location - b).abs() < T::from(1e-10).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dirichlet_bc() {
        let bc = BoundaryCondition::dirichlet(0.0, 1.0);
        assert_eq!(bc.bc_type, BoundaryConditionType::Dirichlet);
        assert_eq!(bc.value, 1.0);
        assert!(bc.is_left(0.0));
    }

    #[test]
    fn test_neumann_bc() {
        let bc = BoundaryCondition::neumann(1.0, 2.0);
        assert_eq!(bc.bc_type, BoundaryConditionType::Neumann);
        assert_eq!(bc.value, 2.0);
        assert!(bc.is_right(1.0));
    }

    #[test]
    fn test_robin_bc() {
        let bc = BoundaryCondition::robin(0.0, 1.0, 2.0, 3.0);
        assert_eq!(bc.bc_type, BoundaryConditionType::Robin);
        assert_eq!(bc.coefficients, vec![1.0, 2.0]);
        assert_eq!(bc.value, 3.0);
    }
}
