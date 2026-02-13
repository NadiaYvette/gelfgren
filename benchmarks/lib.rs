/// Gelfgren Benchmarks
///
/// This library provides comprehensive benchmarks for the Gelfgren piecewise
/// rational interpolation library.
///
/// # Benchmark Categories
///
/// 1. **Boundary Value Problems (BVP)**: Solving differential equations
/// 2. **Special Functions**: Approximating standard mathematical functions
/// 3. **Convergence Analysis**: Comparing polynomial vs rational approximants
///
/// # Structure
///
/// - `bvp/`: Boundary value problem benchmarks
/// - `special_functions/`: Special function approximation benchmarks
/// - `convergence/`: Convergence analysis tools
/// - `python/`: Python scripts for running studies and generating reports
/// - `data/`: Generated benchmark data (JSON)
/// - `reports/`: Generated reports (LaTeX, PDF)

pub mod bvp;
pub mod special_functions;

pub use bvp::*;
pub use special_functions::*;
