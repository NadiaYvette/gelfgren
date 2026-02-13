//! Boundary value problem solver.
//!
//! Solves ordinary differential equations (ODEs) with boundary conditions
//! using piecewise rational approximation.
//!
//! # Mathematical Background
//!
//! A boundary value problem (BVP) consists of:
//! - An ODE: L[u](x) = f(x) on [a, b]
//! - Boundary conditions at endpoints
//!
//! where L is a differential operator, e.g.:
//! - L[u] = u''(x) + p(x)u'(x) + q(x)u(x)  (second-order linear)
//! - L[u] = u^(n)(x) + ...                  (n-th order)
//!
//! # Boundary Conditions
//!
//! - **Dirichlet**: u(a) = α, u(b) = β
//! - **Neumann**: u'(a) = α, u'(b) = β
//! - **Robin**: c₀u(a) + c₁u'(a) = α
//! - **Periodic**: u(a) = u(b), u'(a) = u'(b)
//!
//! # Solution Method
//!
//! 1. Choose mesh τ on [a, b]
//! 2. Represent solution as piecewise rational S_{n,m}(τ)
//! 3. Enforce ODE at collocation points (e.g., mesh points)
//! 4. Enforce boundary conditions
//! 5. Solve resulting linear system for coefficients
//!
//! # Algorithm Components
//!
//! - **Collocation**: Enforce L[S] = f at chosen points
//! - **Galerkin**: Weak formulation with test functions
//! - **Shooting**: Convert BVP to initial value problems
//! - **Finite Differences**: Approximate derivatives on mesh
//!
//! # Example
//!
//! Solve u''(x) + u(x) = 0 with u(0) = 0, u(π) = 0:
//! ```text
//! 1. Create mesh on [0, π]
//! 2. Set up collocation at mesh points
//! 3. Build constraint matrix from u'' + u
//! 4. Add boundary condition rows
//! 5. Solve linear system
//! 6. Construct piecewise rational solution
//! ```
//!
//! # References
//!
//! - Gelfgren (1975): Piecewise rational for BVPs
//! - Ascher, Mattheij & Russell (1995): "Numerical Solution of BVPs for ODEs"
//! - Keller (1968): Numerical methods for two-point BVPs

mod problem;
mod boundary;
mod solver;

pub use problem::{BoundaryValueProblem, DifferentialOperator, RightHandSide, SecondOrderLinear, FunctionRHS};
pub use boundary::{BoundaryCondition, BoundaryConditionType};
pub use solver::BVPSolver;
