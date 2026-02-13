# gelfgren: Piecewise Rational Interpolation for R

R bindings for the Gelfgren piecewise rational interpolation library.

## Installation

### From source

```r
# Install rextendr
install.packages("rextendr")

# Build and install
rextendr::document()
devtools::install()
```

## Quick Start

```r
library(gelfgren)

# Create a Bernstein polynomial
p <- bernstein_polynomial(c(1, 2, 3), 0, 1)

# Evaluate at points
x <- seq(0, 1, by = 0.1)
y <- p$evaluate(x)

# Visualize
plot(x, y, type = "b", main = "Bernstein Polynomial")

# Compute derivative
p_prime <- p$derivative()
y_prime <- p_prime$evaluate(x)
lines(x, y_prime, col = "red")

# Polynomial arithmetic
q <- bernstein_polynomial(c(1, 1), 0, 1)
sum_poly <- p$add(q)
```

## Features

- **BernsteinPolynomial**: Numerically stable polynomial operations
- **RationalFunction**: P(x)/Q(x) with automatic pole detection
- **PadeApproximant**: Rational approximations from power series
- **Mesh**: Uniform and Chebyshev mesh generation

## Examples

### Rational Functions

```r
# Create R(x) = (1 + x) / (1 + 2x)
num <- bernstein_polynomial(c(1, 1), 0, 1)
den <- bernstein_polynomial(c(1, 2), 0, 1)
r <- rational_function(num, den)

x <- seq(0, 1, length.out = 50)
y <- r$evaluate(x)
plot(x, y, type = "l", main = "Rational Function")
```

### Padé Approximants

```r
# Approximate exp(x) with [2/2] Padé
coeffs <- c(1, 1, 0.5, 1/6, 1/24)
pade <- pade_approximant(coeffs, 2, 2, 0, -1, 1)

x <- seq(-1, 1, by = 0.1)
approx <- pade$evaluate(x)
exact <- exp(x)

plot(x, exact, type = "l", col = "blue", lwd = 2,
     main = "Padé Approximant for exp(x)")
points(x, approx, col = "red", pch = 19)
legend("topleft", c("exp(x)", "Padé [2/2]"),
       col = c("blue", "red"), lty = c(1, NA), pch = c(NA, 19))
```

## Testing

```r
devtools::test()
```

## Documentation

```r
# View package documentation
?gelfgren

# View function help
?bernstein_polynomial
?rational_function
?pade_approximant
```

## System Requirements

- R (≥ 3.6.0)
- Rust (≥ 1.70)
- Cargo (Rust package manager)

## License

MIT OR Apache-2.0
