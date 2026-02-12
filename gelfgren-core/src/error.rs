use std::fmt;

/// Error types for Gelfgren numerical computing operations
#[derive(Debug, Clone, PartialEq)]
pub enum GelfgrenError {
    /// Invalid dimensions for matrix/vector operations
    InvalidDimension(String),

    /// Matrix is singular and cannot be inverted
    SingularMatrix,

    /// Numerical method failed to converge
    ConvergenceFailure(String),

    /// Invalid argument provided to a function
    InvalidArgument(String),

    /// Operation would result in numerical overflow
    Overflow,

    /// Operation would result in numerical underflow
    Underflow,

    /// Division by zero
    DivisionByZero,

    /// Feature not yet implemented
    NotImplemented(String),
}

impl fmt::Display for GelfgrenError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDimension(msg) => write!(f, "Invalid dimension: {}", msg),
            Self::SingularMatrix => write!(f, "Matrix is singular"),
            Self::ConvergenceFailure(msg) => write!(f, "Convergence failure: {}", msg),
            Self::InvalidArgument(msg) => write!(f, "Invalid argument: {}", msg),
            Self::Overflow => write!(f, "Numerical overflow"),
            Self::Underflow => write!(f, "Numerical underflow"),
            Self::DivisionByZero => write!(f, "Division by zero"),
            Self::NotImplemented(msg) => write!(f, "Not implemented: {}", msg),
        }
    }
}

impl std::error::Error for GelfgrenError {}

/// Result type for Gelfgren operations
pub type Result<T> = std::result::Result<T, GelfgrenError>;
