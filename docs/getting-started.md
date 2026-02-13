# Getting Started with Gelfgren

This guide will help you get started with Gelfgren across all supported programming languages.

## Table of Contents

- [What is Gelfgren?](#what-is-gelfgren)
- [Installation](#installation)
- [Basic Concepts](#basic-concepts)
- [Language-Specific Guides](#language-specific-guides)
  - [Rust](#rust)
  - [Python](#python)
  - [C++](#cpp)
  - [Java](#java)
  - [R](#r)
  - [C](#c)
- [Common Tasks](#common-tasks)
- [Next Steps](#next-steps)

## What is Gelfgren?

Gelfgren is a numerical computing library that implements piecewise rational interpolation methods. It provides:

- **Bernstein Polynomials**: Numerically stable polynomial representation
- **Rational Functions**: Ratio of polynomials (P(x)/Q(x))
- **Padé Approximants**: Rational approximation of transcendental functions
- **Hermite Interpolation**: Interpolation matching function values and derivatives
- **Piecewise Methods**: Construct approximations on mesh intervals
- **BVP Solvers**: Boundary value problems for differential equations

The library is written in Rust with a C FFI layer that enables bindings to many languages.

## Installation

### Rust

Add to your `Cargo.toml`:

```toml
[dependencies]
gelfgren-core = "0.1.0"
```

Or use cargo:

```bash
cargo add gelfgren-core
```

### Python

```bash
pip install gelfgren
```

For development from source:

```bash
cd bindings/python
pip install maturin
maturin develop --release
```

### C++

1. Download the appropriate native library from [releases](https://github.com/yourusername/gelfgren/releases)
2. Copy `bindings/cpp/gelfgren.hpp` to your include path
3. Link against the native library

Using CMake:

```cmake
find_package(gelfgren REQUIRED)
target_link_libraries(your_target gelfgren::gelfgren)
```

### Java

Maven `pom.xml`:

```xml
<dependencies>
    <dependency>
        <groupId>org.gelfgren</groupId>
        <artifactId>gelfgren</artifactId>
        <version>0.1.0</version>
    </dependency>
</dependencies>
```

Gradle `build.gradle`:

```gradle
dependencies {
    implementation 'org.gelfgren:gelfgren:0.1.0'
}
```

### R

From CRAN (once published):

```r
install.packages("gelfgren")
```

From source:

```r
# Install dependencies
install.packages(c("rextendr", "devtools"))

# Build and install
setwd("bindings/r")
rextendr::document()
devtools::install()
```

### C

Download the native library from [releases](https://github.com/yourusername/gelfgren/releases) and use with your build system:

```bash
# With pkg-config
gcc myprogram.c $(pkg-config --cflags --libs gelfgren) -o myprogram

# Or manually
gcc myprogram.c -I/path/to/include -L/path/to/lib -lgelfgren -o myprogram
```

## Basic Concepts

### Bernstein Polynomials

Bernstein polynomials provide a numerically stable representation for polynomials on an interval [a, b]. A Bernstein polynomial of degree n is defined by (n+1) coefficients.

**Key properties:**
- Better numerical conditioning than power basis
- Coefficients have geometric meaning (control points)
- Efficient degree elevation and evaluation
- Natural bounds on function values

### Rational Functions

A rational function is the ratio of two polynomials: R(x) = P(x)/Q(x). Gelfgren represents both P and Q as Bernstein polynomials.

**Key properties:**
- Can approximate functions with poles
- More flexible than polynomials alone
- Requires careful handling of singularities

### Padé Approximants

A Padé approximant [n/m] is a rational function constructed to match a power series or derivative values up to order n+m.

**Uses:**
- Approximate transcendental functions (exp, sin, log)
- Improve convergence of Taylor series
- Analytic continuation

## Language-Specific Guides

### Rust

#### Creating a Bernstein Polynomial

```rust
use gelfgren_core::bernstein::BernsteinPolynomial;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create polynomial: 2x^2 + 3x + 1 on [0, 1]
    // Bernstein coefficients: [1, 2, 6]
    let coeffs = vec![1.0, 2.0, 6.0];
    let poly = BernsteinPolynomial::new(coeffs, 0.0, 1.0)?;

    // Evaluate at a point
    let y = poly.evaluate(0.5);
    println!("f(0.5) = {}", y);

    // Compute derivative
    let dpoly = poly.derivative();
    let dy = dpoly.evaluate(0.5);
    println!("f'(0.5) = {}", dy);

    // Integrate from 0 to 1
    let integral = poly.integral();
    let area = integral.evaluate(1.0) - integral.evaluate(0.0);
    println!("∫f(x)dx from 0 to 1 = {}", area);

    Ok(())
}
```

#### Creating a Rational Function

```rust
use gelfgren_core::{bernstein::BernsteinPolynomial, rational::RationalFunction};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Numerator: x
    let num = BernsteinPolynomial::new(vec![0.0, 1.0], 0.0, 1.0)?;

    // Denominator: 1 + x
    let den = BernsteinPolynomial::new(vec![1.0, 2.0], 0.0, 1.0)?;

    // Create rational function: x / (1 + x)
    let rat = RationalFunction::new(num, den)?;

    // Evaluate
    let y = rat.evaluate(0.5)?;
    println!("f(0.5) = {}", y);

    // Differentiate
    let drat = rat.derivative();

    Ok(())
}
```

#### Creating a Padé Approximant

```rust
use gelfgren_core::pade::PadeApproximant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Approximate exp(x) near x=0 with [2/2] Padé approximant
    // Power series: 1 + x + x^2/2 + x^3/6 + x^4/24
    let coeffs = vec![1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0];

    let pade = PadeApproximant::from_power_series(coeffs, 2, 2, -1.0, 1.0)?;

    // Evaluate
    let x = 0.5;
    let approx = pade.evaluate(x)?;
    let exact = x.exp();
    println!("Padé: {}, Exact: {}, Error: {}", approx, exact, (approx - exact).abs());

    Ok(())
}
```

### Python

#### Basic Usage

```python
import gelfgren as gf
import numpy as np
import matplotlib.pyplot as plt

# Create a Bernstein polynomial
coeffs = [1.0, 2.0, 6.0]
poly = gf.BernsteinPolynomial(coeffs, a=0.0, b=1.0)

# Evaluate at a single point
y = poly.evaluate(np.array([0.5]))
print(f"f(0.5) = {y[0]}")

# Vectorized evaluation
x = np.linspace(0, 1, 100)
y = poly.evaluate(x)

# Plot
plt.plot(x, y, label='Polynomial')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.show()

# Derivative
dpoly = poly.derivative()
dy = dpoly.evaluate(x)

plt.plot(x, y, label='f(x)')
plt.plot(x, dy, label="f'(x)")
plt.legend()
plt.show()
```

#### Rational Functions

```python
import gelfgren as gf
import numpy as np

# Create rational function x / (1 + x)
num = gf.BernsteinPolynomial([0.0, 1.0], a=0.0, b=1.0)
den = gf.BernsteinPolynomial([1.0, 2.0], a=0.0, b=1.0)
rat = gf.RationalFunction(num, den)

# Evaluate
x = np.linspace(0, 1, 100)
y = rat.evaluate(x)

# Handle potential poles
try:
    y_at_pole = rat.evaluate(np.array([-1.0]))
except ValueError as e:
    print(f"Pole detected: {e}")
```

### C++

#### Basic Usage with RAII

```cpp
#include "gelfgren.hpp"
#include <iostream>
#include <vector>

int main() {
    try {
        // Create polynomial (automatic memory management)
        std::vector<double> coeffs = {1.0, 2.0, 6.0};
        gelfgren::BernsteinPolynomial poly(coeffs, 0.0, 1.0);

        // Evaluate using operator()
        double y = poly(0.5);
        std::cout << "f(0.5) = " << y << std::endl;

        // Derivative
        auto dpoly = poly.derivative();
        double dy = dpoly(0.5);
        std::cout << "f'(0.5) = " << dy << std::endl;

        // Integral
        auto integral = poly.integral();
        double area = integral(1.0) - integral(0.0);
        std::cout << "Area = " << area << std::endl;

    } catch (const gelfgren::GelfgrenException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    // Automatic cleanup on scope exit
    return 0;
}
```

#### Rational Functions

```cpp
#include "gelfgren.hpp"
#include <iostream>

int main() {
    try {
        gelfgren::BernsteinPolynomial num({0.0, 1.0}, 0.0, 1.0);
        gelfgren::BernsteinPolynomial den({1.0, 2.0}, 0.0, 1.0);

        gelfgren::RationalFunction rat(num, den);

        double y = rat(0.5);
        std::cout << "f(0.5) = " << y << std::endl;

    } catch (const gelfgren::GelfgrenException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
```

### Java

#### Basic Usage with Try-With-Resources

```java
import org.gelfgren.BernsteinPolynomial;
import org.gelfgren.GelfgrenException;

public class Example {
    public static void main(String[] args) {
        // Try-with-resources ensures cleanup
        try (BernsteinPolynomial poly =
                new BernsteinPolynomial(new double[]{1.0, 2.0, 6.0}, 0.0, 1.0)) {

            // Evaluate
            double y = poly.evaluate(0.5);
            System.out.printf("f(0.5) = %.6f%n", y);

            // Derivative
            try (BernsteinPolynomial dpoly = poly.derivative()) {
                double dy = dpoly.evaluate(0.5);
                System.out.printf("f'(0.5) = %.6f%n", dy);
            }

            // Integral
            try (BernsteinPolynomial integral = poly.integral()) {
                double area = integral.evaluate(1.0) - integral.evaluate(0.0);
                System.out.printf("Area = %.6f%n", area);
            }

        } catch (GelfgrenException e) {
            System.err.println("Error: " + e.getMessage());
            System.exit(1);
        }
    }
}
```

#### Working with Arrays

```java
import org.gelfgren.BernsteinPolynomial;

public class VectorExample {
    public static void main(String[] args) {
        try (BernsteinPolynomial poly =
                new BernsteinPolynomial(new double[]{1.0, 2.0, 6.0}, 0.0, 1.0)) {

            // Create array of evaluation points
            double[] x = new double[100];
            for (int i = 0; i < x.length; i++) {
                x[i] = i / 99.0;
            }

            // Evaluate at all points
            for (double xi : x) {
                double y = poly.evaluate(xi);
                System.out.printf("f(%.2f) = %.6f%n", xi, y);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
```

### R

#### Basic Usage

```r
library(gelfgren)

# Create a Bernstein polynomial
poly <- BernsteinPolynomial$new(c(1.0, 2.0, 6.0), a = 0.0, b = 1.0)

# Evaluate at a point
y <- poly$evaluate(0.5)
print(paste("f(0.5) =", y))

# Vectorized evaluation
x <- seq(0, 1, length.out = 100)
y <- poly$evaluate(x)

# Plot
plot(x, y, type = "l", xlab = "x", ylab = "f(x)", main = "Bernstein Polynomial")

# Derivative
dpoly <- poly$derivative()
dy <- dpoly$evaluate(x)

# Plot both
plot(x, y, type = "l", col = "blue", xlab = "x", ylab = "y")
lines(x, dy, col = "red")
legend("topleft", legend = c("f(x)", "f'(x)"), col = c("blue", "red"), lty = 1)

# Integral
integral_poly <- poly$integral()
area <- integral_poly$evaluate(1.0) - integral_poly$evaluate(0.0)
print(paste("Area =", area))
```

#### Rational Functions

```r
library(gelfgren)

# Create rational function x / (1 + x)
num <- BernsteinPolynomial$new(c(0.0, 1.0), a = 0.0, b = 1.0)
den <- BernsteinPolynomial$new(c(1.0, 2.0), a = 0.0, b = 1.0)
rat <- RationalFunction$new(num, den)

# Evaluate
x <- seq(0, 1, length.out = 100)
y <- rat$evaluate(x)

# Plot
plot(x, y, type = "l", xlab = "x", ylab = "f(x)", main = "Rational Function")
```

### C

#### Basic Usage

```c
#include "gelfgren.h"
#include <stdio.h>
#include <stdlib.h>

int main() {
    // Create polynomial
    double coeffs[] = {1.0, 2.0, 6.0};
    GelfgrenBernstein* poly = gelfgren_bernstein_create(coeffs, 2, 0.0, 1.0);

    if (poly == NULL) {
        fprintf(stderr, "Error: %s\n", gelfgren_last_error_message());
        return 1;
    }

    // Evaluate
    double result;
    GelfgrenErrorCode err = gelfgren_bernstein_evaluate(poly, 0.5, &result);

    if (err != GELFGREN_SUCCESS) {
        fprintf(stderr, "Error: %s\n", gelfgren_last_error_message());
        gelfgren_bernstein_free(poly);
        return 1;
    }

    printf("f(0.5) = %f\n", result);

    // Derivative
    GelfgrenBernstein* dpoly = gelfgren_bernstein_derivative(poly);
    if (dpoly != NULL) {
        gelfgren_bernstein_evaluate(dpoly, 0.5, &result);
        printf("f'(0.5) = %f\n", result);
        gelfgren_bernstein_free(dpoly);
    }

    // Clean up
    gelfgren_bernstein_free(poly);

    return 0;
}
```

## Common Tasks

### Approximating a Function

**Python:**

```python
import gelfgren as gf
import numpy as np

# Sample a function
def f(x):
    return np.sin(x)

# Create mesh points
x_data = np.linspace(0, np.pi, 10)
y_data = f(x_data)

# Fit Bernstein polynomial (convert from samples)
# Note: This requires Hermite or least squares fitting (future feature)
# For now, manually provide Bernstein coefficients
```

### Computing Derivatives

All language bindings provide a `derivative()` method that returns a new polynomial/rational representing the derivative.

### Computing Integrals

Use the `integral()` method to get the antiderivative, then evaluate at bounds:

```python
# Python
area = poly.integral().evaluate(b) - poly.integral().evaluate(a)
```

```cpp
// C++
double area = poly.integral()(b) - poly.integral()(a);
```

### Checking for Poles

When working with rational functions, always check for potential poles:

```python
# Python
try:
    y = rat.evaluate(x)
except ValueError:
    print("Pole detected at x =", x)
```

```cpp
// C++
try {
    double y = rat(x);
} catch (const gelfgren::GelfgrenException& e) {
    std::cerr << "Pole at x = " << x << std::endl;
}
```

## Next Steps

- **[Mathematical Background](mathematics.md)**: Learn the theory behind the algorithms
- **[Architecture Guide](architecture.md)**: Understand the library structure
- **[API Reference](https://yourusername.github.io/gelfgren)**: Complete API documentation
- **[Examples](../examples/)**: Browse complete working examples
- **[Contributing](../CONTRIBUTING.md)**: Help improve Gelfgren

## Further Reading

- Gelfgren, J. (1975). "Piecewise Rational Interpolation". *BIT Numerical Mathematics*, 15, 382-393.
- Farouki, R. T., & Rajan, V. T. (1987). "Algorithms for Polynomials in Bernstein Form". *Computer Aided Geometric Design*, 5(1), 1-26.
- Traub, J. F. (1964). "On Lagrange-Hermite Interpolation". *SIAM Journal on Numerical Analysis*, 1(1), 1-15.
