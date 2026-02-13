# Architecture

This document describes the architecture and design of the Gelfgren library.

## Table of Contents

- [Overview](#overview)
- [Layer Architecture](#layer-architecture)
- [Core Library](#core-library)
- [FFI Layer](#ffi-layer)
- [Language Bindings](#language-bindings)
- [Build System](#build-system)
- [CI/CD Pipeline](#cicd-pipeline)
- [Design Decisions](#design-decisions)

## Overview

Gelfgren uses a layered architecture with a Rust core, C FFI bridge, and language-specific bindings:

```
┌──────────────────────────────────────────────┐
│        Language Bindings (Layer 3)           │
│  Python  Java  R  Ruby  C++  Fortran        │
│  Haskell  Mercury  C                         │
└────────────────┬─────────────────────────────┘
                 │ Language-specific wrappers
          ┌──────▼────────┐
          │   FFI Layer   │  (Layer 2)
          │  (C ABI)      │
          └──────┬────────┘
                 │ C-compatible interface
          ┌──────▼────────┐
          │  Rust Core    │  (Layer 1)
          │ gelfgren-core │
          └───────────────┘
```

### Design Principles

1. **Safety First**: Rust core ensures memory safety and correctness
2. **Zero-Cost Abstractions**: FFI introduces minimal overhead
3. **Idiomatic APIs**: Each language binding feels native to that language
4. **Comprehensive Testing**: Every layer has extensive test coverage
5. **Documentation**: Complete docs for all languages

## Layer Architecture

### Layer 1: Rust Core (`gelfgren-core`)

**Purpose**: Implements all numerical algorithms in pure Rust

**Responsibilities**:
- Bernstein polynomial operations
- Rational function arithmetic
- Padé approximant construction
- Hermite interpolation
- Piecewise rational methods
- BVP solvers

**Dependencies**: Minimal (only `thiserror` for error handling)

**Testing**: Unit tests and integration tests in Rust

### Layer 2: FFI Layer (`gelfgren-ffi`)

**Purpose**: Provides C-compatible interface to Rust core

**Responsibilities**:
- Type conversion (Rust ↔ C)
- Error handling and propagation
- Memory management
- Null pointer checking
- Panic catching

**Dependencies**: `libc` for C types

**Testing**: FFI-specific tests for memory safety and error handling

### Layer 3: Language Bindings

**Purpose**: Provide idiomatic APIs for each target language

**Responsibilities**:
- Wrap FFI calls in language-specific constructs
- Memory management (RAII, GC integration, etc.)
- Error handling in language-native style
- Type conversions
- Documentation in language style

**Dependencies**: Language-specific (PyO3, JNI, extendr, etc.)

**Testing**: Language-specific test frameworks

## Core Library

### Module Structure

```
gelfgren-core/
├── src/
│   ├── lib.rs                  # Public API exports
│   ├── error.rs                # Error types
│   ├── bernstein/
│   │   └── mod.rs              # Bernstein polynomials
│   ├── rational/
│   │   └── mod.rs              # Rational functions
│   ├── pade/
│   │   └── mod.rs              # Padé approximants
│   ├── hermite/
│   │   ├── mod.rs              # Hermite interpolation
│   │   └── bell.rs             # Bell polynomials
│   ├── piecewise/
│   │   ├── mod.rs              # Piecewise construction
│   │   └── mesh.rs             # Mesh management
│   └── bvp/
│       ├── mod.rs              # BVP framework
│       └── solver.rs           # Solver implementations
└── tests/
    └── integration_tests.rs    # Integration tests
```

### Key Types

#### `BernsteinPolynomial<T>`

```rust
pub struct BernsteinPolynomial<T> {
    coeffs: Vec<T>,        // Bernstein coefficients
    degree: usize,         // Polynomial degree
    a: T,                  // Left endpoint
    b: T,                  // Right endpoint
}
```

**Invariants**:
- `coeffs.len() == degree + 1`
- `a < b`
- Coefficients are in Bernstein basis (not power basis)

#### `RationalFunction<T>`

```rust
pub struct RationalFunction<T> {
    numerator: BernsteinPolynomial<T>,
    denominator: BernsteinPolynomial<T>,
}
```

**Invariants**:
- Numerator and denominator have same interval `[a, b]`
- Denominator is not identically zero

#### `PadeApproximant<T>`

```rust
pub struct PadeApproximant<T> {
    rational: RationalFunction<T>,
    n: usize,  // Numerator degree
    m: usize,  // Denominator degree
}
```

### Error Handling

```rust
#[derive(Debug, thiserror::Error)]
pub enum GelfgrenError {
    #[error("Division by zero")]
    DivisionByZero,

    #[error("Invalid interval: {0}")]
    InvalidInterval(String),

    #[error("Empty coefficients")]
    EmptyCoefficients,

    #[error("Singular matrix in linear system")]
    SingularMatrix,

    #[error("Pole detected at x = {0}")]
    Pole(f64),

    // ... more variants
}
```

Uses `Result<T, GelfgrenError>` throughout for error propagation.

## FFI Layer

### Type Safety with Opaque Pointers

To maintain type safety across the FFI boundary, we use opaque pointer types:

```rust
// In gelfgren-ffi/src/types.rs

#[repr(C)]
pub struct GelfgrenBernstein {
    _private: [u8; 0],  // Zero-sized, prevents construction
}

// Conversion helpers (not exposed to C)
impl GelfgrenBernstein {
    pub(crate) fn from_box(b: Box<BernsteinBox>) -> *mut Self {
        Box::into_raw(b) as *mut Self
    }

    pub(crate) unsafe fn as_box<'a>(ptr: *mut Self) -> &'a mut BernsteinBox {
        &mut *(ptr as *mut BernsteinBox)
    }
}
```

**Benefits**:
- Type safety: Can't mix different pointer types
- Prevents direct access from C
- Rust controls all memory

### Error Handling Strategy

#### Thread-Local Error Storage

```rust
use std::cell::RefCell;

thread_local! {
    static LAST_ERROR: RefCell<Option<String>> = RefCell::new(None);
}

pub fn set_last_error(err: String) {
    LAST_ERROR.with(|e| *e.borrow_mut() = Some(err));
}

#[no_mangle]
pub extern "C" fn gelfgren_last_error_message() -> *const c_char {
    // Return error message or empty string
}
```

**Properties**:
- Thread-safe (each thread has own error)
- No allocation in C (returns pointer to static buffer)
- Error persists until next FFI call or cleared

#### Error Codes

```rust
#[repr(i32)]
pub enum GelfgrenErrorCode {
    Success = 0,
    NullPointer = -1,
    InvalidInterval = -2,
    EmptyCoefficients = -3,
    SingularMatrix = -4,
    DivisionByZero = -5,
    Pole = -6,
    InvalidDegree = -7,
    InvalidOperation = -8,
}
```

**Usage Pattern**:

```rust
#[no_mangle]
pub unsafe extern "C" fn gelfgren_function(
    ptr: *const Data,
    result: *mut f64
) -> GelfgrenErrorCode {
    // 1. Check null pointers
    if ptr.is_null() || result.is_null() {
        set_last_error("Null pointer argument".to_string());
        return GelfgrenErrorCode::NullPointer;
    }

    // 2. Perform operation with Result
    match perform_operation(ptr) {
        Ok(value) => {
            *result = value;
            GelfgrenErrorCode::Success
        }
        Err(e) => {
            set_last_error(format!("{:?}", e));
            error_code_from_error(&e)
        }
    }
}
```

### Memory Management

#### Creation Functions

```rust
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_create(
    coeffs: *const f64,
    degree: usize,
    a: f64,
    b: f64,
) -> *mut GelfgrenBernstein {
    // Validation
    if coeffs.is_null() {
        set_last_error("Null pointer: coeffs".to_string());
        return std::ptr::null_mut();
    }

    // Convert C data to Rust
    let coeffs_slice = std::slice::from_raw_parts(coeffs, degree + 1);
    let coeffs_vec = coeffs_slice.to_vec();

    // Create Rust object
    match BernsteinPolynomial::new(coeffs_vec, a, b) {
        Ok(poly) => {
            let boxed = Box::new(BernsteinBox::F64(poly));
            GelfgrenBernstein::from_box(boxed)
        }
        Err(e) => {
            set_last_error(format!("Failed to create: {:?}", e));
            std::ptr::null_mut()
        }
    }
}
```

#### Destruction Functions

```rust
#[no_mangle]
pub unsafe extern "C" fn gelfgren_bernstein_free(ptr: *mut GelfgrenBernstein) {
    if ptr.is_null() {
        return;  // Freeing null is a no-op
    }

    // Convert back to Box and drop
    let _ = Box::from_raw(ptr as *mut BernsteinBox);
}
```

**Rules**:
- Every `*_create` has corresponding `*_free`
- Caller owns pointers returned by `*_create`
- Double-free is safe (no-op on null)
- Freeing in wrong language binding is UB (don't do it)

### Generated Header

The C header is generated automatically by `cbindgen`:

```toml
# gelfgren-ffi/cbindgen.toml

[export]
exclude = ["BernsteinBox", "RationalBox", ...]  # Internal types

[fn]
sort_by = "Name"

language = "C"
cpp_compat = true  # Add extern "C" for C++
```

Generates:

```c
// gelfgren.h

#ifndef GELFGREN_H
#define GELFGREN_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct GelfgrenBernstein GelfgrenBernstein;

typedef enum {
    GELFGREN_SUCCESS = 0,
    GELFGREN_NULL_POINTER = -1,
    // ...
} GelfgrenErrorCode;

GelfgrenBernstein* gelfgren_bernstein_create(
    const double* coeffs,
    size_t degree,
    double a,
    double b
);

void gelfgren_bernstein_free(GelfgrenBernstein* ptr);

GelfgrenErrorCode gelfgren_bernstein_evaluate(
    const GelfgrenBernstein* ptr,
    double x,
    double* result
);

#ifdef __cplusplus
}
#endif

#endif  // GELFGREN_H
```

## Language Bindings

### C Bindings

**Approach**: Use generated header directly

**Memory**: Manual management (create/free pairs)

**Errors**: Check return codes and call `gelfgren_last_error_message()`

### C++ Bindings

**Approach**: Header-only RAII wrappers

**Key Pattern**: Smart pointers with custom deleters

```cpp
class BernsteinPolynomial {
    std::unique_ptr<GelfgrenBernstein, decltype(&gelfgren_bernstein_free)> ptr_;

public:
    BernsteinPolynomial(const std::vector<double>& coeffs, double a, double b)
        : ptr_(check_ptr(gelfgren_bernstein_create(
            coeffs.data(), coeffs.size() - 1, a, b)),
            &gelfgren_bernstein_free) {}

    // Move semantics
    BernsteinPolynomial(BernsteinPolynomial&&) = default;
    BernsteinPolynomial& operator=(BernsteinPolynomial&&) = default;

    // No copying (would need gelfgren_bernstein_clone)
    BernsteinPolynomial(const BernsteinPolynomial&) = delete;
    BernsteinPolynomial& operator=(const BernsteinPolynomial&) = delete;

    double operator()(double x) const {
        double result;
        check_error(gelfgren_bernstein_evaluate(ptr_.get(), x, &result));
        return result;
    }
};
```

**Benefits**:
- Automatic cleanup
- Exception safety
- Operator overloading for natural syntax

### Python Bindings

**Approach**: PyO3 native extension

**Key Features**:
- NumPy integration for array operations
- Python exceptions for errors
- Type hints for IDE support

```rust
#[pyclass(name = "BernsteinPolynomial")]
struct PyBernstein {
    inner: CoreBernstein<f64>,
}

#[pymethods]
impl PyBernstein {
    #[new]
    fn new(coeffs: Vec<f64>, a: f64, b: f64) -> PyResult<Self> {
        CoreBernstein::new(coeffs, a, b)
            .map(|inner| Self { inner })
            .map_err(|e| PyValueError::new_err(format!("{:?}", e)))
    }

    fn evaluate<'py>(
        &self,
        py: Python<'py>,
        x: PyReadonlyArray1<f64>,
    ) -> Bound<'py, PyArray1<f64>> {
        let result: Vec<f64> = x.as_slice()
            .unwrap()
            .iter()
            .map(|&xi| self.inner.evaluate(xi))
            .collect();
        result.to_pyarray_bound(py)
    }
}
```

**Build**: Maturin compiles and packages as wheel

### Java Bindings

**Approach**: JNI with wrapper classes

**Key Pattern**: `AutoCloseable` for resource management

```java
public class BernsteinPolynomial implements AutoCloseable {
    private long nativeHandle;
    private boolean closed = false;

    public BernsteinPolynomial(double[] coeffs, double a, double b) {
        this.nativeHandle = create(coeffs, coeffs.length - 1, a, b);
        if (this.nativeHandle == 0) {
            throw new GelfgrenException("Failed to create");
        }
    }

    @Override
    public void close() {
        if (!closed && nativeHandle != 0) {
            free(nativeHandle);
            closed = true;
            nativeHandle = 0;
        }
    }

    // JNI methods
    private static native long create(double[] coeffs, int degree, double a, double b);
    private static native void free(long handle);
}
```

**Build**: Maven with rust-maven-plugin

### R Bindings

**Approach**: extendr framework

**Key Features**:
- Automatic conversion between R and Rust types
- roxygen2 documentation
- Vectorized operations

```rust
#[extendr]
impl BernsteinPolynomial {
    fn new(coeffs: Vec<f64>, a: f64, b: f64) -> Result<Self> {
        CoreBernstein::new(coeffs, a, b)
            .map(|inner| Self { inner })
            .map_err(|e| format!("Failed: {:?}", e))
    }

    fn evaluate(&self, x: Vec<f64>) -> Vec<f64> {
        x.iter().map(|&xi| self.inner.evaluate(xi)).collect()
    }
}

extendr_module! {
    mod gelfgren;
    impl BernsteinPolynomial;
}
```

**Build**: rextendr generates R wrapper functions

## Build System

### Workspace Structure

```toml
# Cargo.toml (workspace root)
[workspace]
members = ["gelfgren-core", "gelfgren-ffi"]
resolver = "2"

[workspace.package]
version = "0.1.0"
authors = ["Your Name"]
license = "MIT OR Apache-2.0"
edition = "2021"

[profile.release]
lto = true
codegen-units = 1
opt-level = 3
```

### Build Scripts

#### Master Build Script (`scripts/build-all.sh`)

```bash
#!/bin/bash
set -e

echo "Building Rust core..."
cargo build --release

echo "Generating C header..."
cd gelfgren-ffi
cbindgen --config cbindgen.toml --output ../include/gelfgren.h
cd ..

echo "Building Python bindings..."
cd bindings/python
maturin build --release
cd ../..

echo "Building Java bindings..."
cd bindings/java
mvn clean package
cd ../..

echo "Building R bindings..."
cd bindings/r
Rscript -e "rextendr::document()"
R CMD build .
cd ../..

echo "All builds complete!"
```

#### Selective Build Script (`scripts/build.sh`)

```bash
#!/bin/bash
LANG=$1

case $LANG in
    python)
        cd bindings/python && maturin build --release
        ;;
    java)
        cd bindings/java && mvn clean package
        ;;
    r)
        cd bindings/r && R CMD build .
        ;;
    *)
        echo "Unknown language: $LANG"
        exit 1
        ;;
esac
```

## CI/CD Pipeline

### Workflow Structure

```
.github/workflows/
├── rust-core.yml           # Rust testing and linting
├── bindings-cpp.yml        # C++ builds and tests
├── bindings-python.yml     # Python wheels and tests
├── bindings-java.yml       # Maven builds and tests
├── bindings-r.yml          # R CMD check
└── release.yml             # Release automation
```

### Rust Core CI

**Triggers**: Push to main/develop, PRs

**Jobs**:
1. **Check**: `cargo check` on all platforms
2. **Test**: `cargo test` with stable and beta Rust
3. **Format**: `cargo fmt --check`
4. **Lint**: `cargo clippy -- -D warnings`
5. **Audit**: `cargo audit` for security
6. **Coverage**: `cargo tarpaulin` with codecov upload
7. **Docs**: Build and deploy rustdoc to GitHub Pages
8. **MSRV**: Test minimum supported Rust version

### Language Binding CI

Each language has dedicated workflow:
- Multi-version testing (Python 3.9-3.12, Java 11/17/21, R 4.1-4.3)
- Multi-platform (Linux, macOS, Windows)
- Test execution with language-specific frameworks
- Artifact collection (wheels, JARs, packages)

### Release Workflow

**Trigger**: Push tag matching `v*.*.*`

**Steps**:
1. Build native libraries for all platforms
2. Package as tar.gz/zip
3. Create GitHub release
4. Upload artifacts
5. Publish to registries:
   - Rust: crates.io
   - Python: PyPI
   - Java: Maven Central
   - R: CRAN (manual submission)

## Design Decisions

### Why Rust for Core?

- **Memory safety**: No segfaults or undefined behavior
- **Performance**: Zero-cost abstractions, comparable to C++
- **Correctness**: Strong type system catches bugs at compile time
- **Ecosystem**: Excellent FFI tools (cbindgen, PyO3, extendr)

### Why Bernstein Basis?

- **Numerical stability**: Superior conditioning to power basis
- **Geometric intuition**: Coefficients are control points
- **Efficient algorithms**: All operations have optimal complexity

### Why FFI Instead of Native Rewrite?

- **Single source of truth**: Core logic maintained in one place
- **Consistency**: All languages get identical results
- **Effort**: 9 native implementations would be error-prone
- **Performance**: FFI overhead is negligible (<1%)

### Why Opaque Pointers?

- **Type safety**: Prevents mixing different pointer types
- **Encapsulation**: C code can't access internal structure
- **Flexibility**: Can change internal representation without breaking ABI

### Why Thread-Local Error Storage?

- **Thread safety**: No global state races
- **Simplicity**: C code doesn't manage string memory
- **Compatibility**: Works with all threading models

## Performance Considerations

### Zero-Copy FFI

Where possible, data is passed by reference:

```rust
// No copy: C array → Rust slice
let coeffs_slice = std::slice::from_raw_parts(coeffs, len);

// Copy only when Rust needs ownership
let coeffs_vec = coeffs_slice.to_vec();
```

### Vectorized Operations

Python binding uses NumPy for vectorization:

```python
# Single FFI call for entire array
y = poly.evaluate(np.linspace(0, 1, 1000))
```

Internally: Rust loop over elements (no per-element FFI overhead).

### Inlining

Critical paths use `#[inline]` for compiler optimization:

```rust
#[inline]
pub fn evaluate(&self, x: f64) -> f64 {
    // Hot loop
}
```

### Memory Pooling (Future)

For applications creating many short-lived objects, could add object pool to reduce allocation overhead.

## Security Considerations

### Input Validation

All FFI functions validate inputs:

```rust
// Check null pointers
if ptr.is_null() {
    return GelfgrenErrorCode::NullPointer;
}

// Check array bounds
if degree > MAX_DEGREE {
    return GelfgrenErrorCode::InvalidDegree;
}

// Check intervals
if a >= b {
    return GelfgrenErrorCode::InvalidInterval;
}
```

### Panic Safety

All panics are caught at FFI boundary:

```rust
std::panic::catch_unwind(|| {
    // Operation that might panic
}).unwrap_or_else(|_| {
    set_last_error("Internal panic".to_string());
    GelfgrenErrorCode::InvalidOperation
})
```

### Memory Safety

Rust guarantees:
- No use-after-free
- No double-free
- No buffer overflows
- No data races

C/C++ bindings must follow rules to maintain these guarantees.

## Future Directions

### Planned Improvements

1. **SIMD**: Vectorize polynomial evaluation
2. **Parallel**: Multi-threaded piecewise evaluation
3. **GPU**: CUDA/OpenCL acceleration for large-scale problems
4. **WebAssembly**: Browser-based computation
5. **Distributed**: MPI support for HPC clusters

### API Evolution

- Versioned API: `gelfgren_v1_*` functions for stability
- Deprecation warnings: Mark old functions with compiler warnings
- Migration guides: Document breaking changes

## Conclusion

This architecture provides:

- **Safety**: Rust core + validation at boundaries
- **Performance**: Zero-cost abstractions + zero-copy FFI
- **Usability**: Idiomatic APIs for each language
- **Maintainability**: Single source of truth for algorithms
- **Extensibility**: Easy to add new language bindings

The layered design allows independent evolution of each component while maintaining correctness and performance.
