/// Boundary Value Problem benchmarks module
///
/// Re-exports from the bvp/ directory

// Include the Poisson 1D module
#[path = "bvp/poisson_1d.rs"]
pub mod poisson_1d;

pub use poisson_1d::*;
