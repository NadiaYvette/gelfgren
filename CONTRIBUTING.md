# Contributing to Gelfgren

Thank you for your interest in contributing to Gelfgren! This document provides guidelines and instructions for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Code Style](#code-style)
- [Testing](#testing)
- [Documentation](#documentation)
- [Pull Request Process](#pull-request-process)
- [Areas for Contribution](#areas-for-contribution)

## Code of Conduct

We are committed to providing a welcoming and inclusive environment. Please:

- Be respectful and considerate of others
- Welcome newcomers and help them get started
- Focus on constructive criticism
- Respect differing viewpoints and experiences

## Getting Started

### Prerequisites

Install the required development tools:

```bash
# Rust toolchain
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Additional tools
cargo install cbindgen cargo-audit cargo-tarpaulin

# For Python bindings
pip install maturin pytest numpy

# For Java bindings
# Install JDK 11+ and Maven 3.6+

# For R bindings
# Install R 4.0+ and run:
# install.packages(c("rextendr", "devtools", "testthat"))
```

### Clone and Build

```bash
# Fork the repository on GitHub, then clone your fork
git clone https://github.com/YOUR_USERNAME/gelfgren.git
cd gelfgren

# Build the core library
cargo build

# Run tests
cargo test --workspace

# Build language bindings (optional)
./scripts/build.sh python
```

## Development Workflow

### 1. Create a Branch

```bash
git checkout -b feature/my-new-feature
# or
git checkout -b fix/issue-123
```

Use descriptive branch names:
- `feature/` for new features
- `fix/` for bug fixes
- `docs/` for documentation
- `refactor/` for refactoring
- `test/` for test improvements

### 2. Make Changes

- Write clear, focused commits
- Include tests for new functionality
- Update documentation as needed
- Ensure code follows style guidelines

### 3. Test Your Changes

```bash
# Run all Rust tests
cargo test --workspace

# Check formatting
cargo fmt --all -- --check

# Run linter
cargo clippy --all-targets --all-features -- -D warnings

# Security audit
cargo audit

# Test specific language bindings
cd bindings/python && pytest tests/
cd bindings/java && mvn test
cd bindings/r && R CMD check .
```

### 4. Commit Your Changes

Write clear commit messages following the [Conventional Commits](https://www.conventionalcommits.org/) format:

```
type(scope): short description

Longer explanation if needed.

Fixes #123
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `test`: Test improvements
- `refactor`: Code refactoring
- `perf`: Performance improvements
- `chore`: Maintenance tasks

Examples:
```bash
git commit -m "feat(bernstein): add degree reduction algorithm"
git commit -m "fix(ffi): handle null pointers in error path"
git commit -m "docs(python): add NumPy integration examples"
```

### 5. Push and Create Pull Request

```bash
git push origin feature/my-new-feature
```

Then create a pull request on GitHub with:
- Clear description of changes
- Reference to related issues
- Screenshots/examples if applicable

## Code Style

### Rust

Follow the [Rust Style Guide](https://rust-lang.github.io/api-guidelines/):

```rust
// Use rustfmt (automatic formatting)
cargo fmt --all

// Good style examples:

// 1. Clear function names and documentation
/// Evaluates the Bernstein polynomial at the given point.
///
/// # Arguments
///
/// * `x` - Point at which to evaluate (must be in [a, b])
///
/// # Returns
///
/// The value of the polynomial at `x`
pub fn evaluate(&self, x: f64) -> f64 {
    // Implementation
}

// 2. Error handling with Result
pub fn new(coeffs: Vec<f64>, a: f64, b: f64) -> Result<Self, GelfgrenError> {
    if coeffs.is_empty() {
        return Err(GelfgrenError::EmptyCoefficients);
    }
    // ...
}

// 3. Clear variable names
let num_points = data.len();
let degree = coeffs.len() - 1;
```

### FFI Code

Follow C FFI best practices:

```rust
// 1. Always check for null pointers
#[no_mangle]
pub unsafe extern "C" fn gelfgren_function(ptr: *const Data) -> GelfgrenErrorCode {
    if ptr.is_null() {
        set_last_error("Null pointer argument".to_string());
        return GelfgrenErrorCode::NullPointer;
    }
    // ...
}

// 2. Use opaque pointer types
#[repr(C)]
pub struct GelfgrenObject {
    _private: [u8; 0],
}

// 3. Provide clear error messages
match result {
    Ok(val) => val,
    Err(e) => {
        set_last_error(format!("Operation failed: {:?}", e));
        return GelfgrenErrorCode::InvalidOperation;
    }
}
```

### Python

Follow [PEP 8](https://peps.python.org/pep-0008/):

```python
# Clear docstrings
def evaluate(self, x: np.ndarray) -> np.ndarray:
    """Evaluate the polynomial at given points.

    Args:
        x: NumPy array of evaluation points

    Returns:
        NumPy array of function values

    Raises:
        ValueError: If x contains values outside the domain
    """
    pass
```

### Other Languages

- **Java**: Follow [Google Java Style Guide](https://google.github.io/styleguide/javaguide.html)
- **C++**: Follow [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/)
- **R**: Follow [Tidyverse Style Guide](https://style.tidyverse.org/)

## Testing

### Rust Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_evaluation() {
        let poly = BernsteinPolynomial::new(vec![1.0, 2.0, 3.0], 0.0, 1.0).unwrap();
        let result = poly.evaluate(0.5);
        assert!((result - 2.0).abs() < 1e-10);
    }

    #[test]
    #[should_panic(expected = "EmptyCoefficients")]
    fn test_empty_coefficients() {
        BernsteinPolynomial::new(vec![], 0.0, 1.0).unwrap();
    }
}
```

### Coverage

We aim for >80% code coverage:

```bash
cargo tarpaulin --workspace --out Html --output-dir coverage/
```

### Integration Tests

Place integration tests in `tests/` directories:

```
bindings/python/tests/test_bernstein.py
bindings/java/src/test/java/org/gelfgren/BernsteinTest.java
bindings/r/tests/testthat/test-bernstein.R
```

## Documentation

### Rust Documentation

Use rustdoc comments with examples:

```rust
/// Computes the derivative of the Bernstein polynomial.
///
/// The derivative is also a Bernstein polynomial of degree n-1.
///
/// # Examples
///
/// ```
/// use gelfgren_core::bernstein::BernsteinPolynomial;
///
/// let poly = BernsteinPolynomial::new(vec![0.0, 1.0, 4.0], 0.0, 1.0).unwrap();
/// let dpoly = poly.derivative();
/// assert_eq!(dpoly.degree(), 1);
/// ```
///
/// # Returns
///
/// A new `BernsteinPolynomial` representing the derivative
pub fn derivative(&self) -> Self {
    // Implementation
}
```

Generate and check docs:

```bash
cargo doc --no-deps --open
```

### Language-Specific Docs

- **Python**: Use Sphinx with autodoc
- **Java**: Use Javadoc
- **R**: Use roxygen2
- **C/C++**: Doxygen comments in headers

### README and Guides

- Update README.md for user-facing changes
- Add tutorials to `docs/` directory
- Include runnable examples in `examples/`

## Pull Request Process

### Before Submitting

- [ ] Tests pass locally: `cargo test --workspace`
- [ ] Code is formatted: `cargo fmt --all`
- [ ] No clippy warnings: `cargo clippy --all-targets -- -D warnings`
- [ ] Documentation is updated
- [ ] CHANGELOG.md is updated (for significant changes)
- [ ] Commit messages follow conventional format

### PR Description Template

```markdown
## Description

Brief description of the changes.

## Type of Change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update

## Testing

Describe the tests you ran and how to reproduce them.

## Checklist

- [ ] My code follows the style guidelines
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] I have made corresponding changes to the documentation
- [ ] My changes generate no new warnings
- [ ] I have added tests that prove my fix is effective or that my feature works
- [ ] New and existing unit tests pass locally with my changes

## Related Issues

Fixes #123
Related to #456
```

### Review Process

1. Automated CI checks must pass
2. At least one maintainer review required
3. Address review comments
4. Squash commits if requested
5. Maintainer will merge when approved

## Areas for Contribution

### High Priority

- **Additional language bindings**: Ruby, Fortran, Haskell, Mercury
- **Performance optimizations**: SIMD instructions, parallel evaluation
- **Documentation**: Tutorials, examples, API reference improvements
- **Testing**: Increase coverage, add benchmarks

### Medium Priority

- **Error handling improvements**: Better error messages, validation
- **API enhancements**: Convenience methods, builder patterns
- **Build system**: Improve cross-compilation, reduce dependencies
- **CI/CD**: Add more platform coverage, automated releases

### Good First Issues

Look for issues labeled `good-first-issue` in the issue tracker. These are suitable for newcomers and include:

- Documentation improvements
- Test additions
- Example programs
- Minor bug fixes

## Getting Help

- **Questions**: Use [GitHub Discussions](https://github.com/yourusername/gelfgren/discussions)
- **Bugs**: Open an [issue](https://github.com/yourusername/gelfgren/issues)
- **Chat**: Join our community (link to chat platform)

## Recognition

Contributors are acknowledged in:
- CONTRIBUTORS.md file
- Release notes
- Project documentation

Thank you for contributing to Gelfgren!
