/// Special function approximation benchmarks module
///
/// Re-exports from the special_functions/ directory

// Include the approximations module
#[path = "special_functions/approximations.rs"]
pub mod approximations;

pub use approximations::*;
