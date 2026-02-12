# Gelfgren.jl

Julia bindings for the Gelfgren numerical computing library, providing high-performance piecewise rational interpolation methods with a native Julia interface.

## Features

- **Native Julia Interface**: Idiomatic Julia API with multiple dispatch
- **Vectorized Operations**: Broadcast evaluation over arrays
- **Automatic Memory Management**: Finalizers handle cleanup automatically
- **Type-Safe**: Leverages Julia's type system for safety
- **Operator Overloading**: Natural mathematical syntax (+, -, *, /)
- **Zero-Cost Abstractions**: Direct ccall to C library for performance
- **Comprehensive Documentation**: Docstrings for all exported functions

## Requirements

- Julia 1.9 or later
- Gelfgren C library (libgelfgren.so/dylib)

## Installation

### From Source

```bash
# Build the Rust library first
cd ../..
cargo build --release

# Start Julia in the bindings directory
cd bindings/julia
julia

# In Julia REPL
julia> ] dev .
julia> using Gelfgren
```

### Setting Library Path

If the library is not in a standard location, set the environment variable:

```bash
export GELFGREN_LIB=/path/to/libgelfgren.so
```

Or in Julia before loading:

```julia
ENV["GELFGREN_LIB"] = "/path/to/libgelfgren.so"
using Gelfgren
```

## Quick Start

```julia
using Gelfgren

# Create a Bernstein polynomial
poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)

# Evaluate at a point
y = evaluate(poly, 0.5)  # Returns 3.0

# Vectorized evaluation
xs = [0.0, 0.25, 0.5, 0.75, 1.0]
ys = evaluate(poly, xs)

# Compute derivative
dpoly = derivative(poly)
dy = evaluate(dpoly, 0.5)

# Arithmetic operations
p1 = BernsteinPolynomial([1.0, 2.0, 3.0], 0.0, 1.0)
p2 = BernsteinPolynomial([2.0, 1.0, 1.0], 0.0, 1.0)
sum_poly = p1 + p2
scaled = 2.0 * p1
```

## Usage

### Bernstein Polynomials

```julia
# Create from coefficients on interval [a, b]
poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)

# Get degree
deg = degree(poly)  # Returns 2

# Evaluate at single point
y = evaluate(poly, 0.5)

# Evaluate at multiple points (vectorized)
xs = range(0, 1, length=100)
ys = evaluate(poly, collect(xs))

# Calculus operations
dpoly = derivative(poly)     # Derivative
ipoly = integral(poly)        # Antiderivative

# Arithmetic operations
p1 + p2        # Addition
p1 - p2        # Subtraction
p1 * p2        # Multiplication
2.0 * p1       # Scalar multiplication
p1 * 2.0       # Also works
```

### Rational Functions

```julia
# Create rational function P(x)/Q(x)
num = BernsteinPolynomial([0.0, 1.0], 0.0, 1.0)
den = BernsteinPolynomial([1.0, 2.0], 0.0, 1.0)
rat = RationalFunction(num, den)

# Evaluate
y = evaluate(rat, 0.5)

# Vectorized
ys = evaluate(rat, [0.0, 0.5, 1.0])

# Derivative (quotient rule)
drat = derivative(rat)
```

### Padé Approximants

```julia
# Approximate exp(x) with [2/2] Padé approximant
coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]
pade = PadeApproximant(coeffs, 2, 2, -1.0, 1.0)

# Evaluate
y = evaluate(pade, 0.5)

# Compare with exact value
exact = exp(0.5)
error = abs(y - exact)
```

### Working with Arrays

```julia
# Using ranges
xs = 0.0:0.1:1.0
ys = evaluate(poly, collect(xs))

# Using arrays
xs = [0.0, 0.25, 0.5, 0.75, 1.0]
ys = evaluate(poly, xs)

# Statistical operations
using Statistics
mean_y = mean(ys)
std_y = std(ys)
```

### Integration with Julia Ecosystem

#### With Plots.jl

```julia
using Plots
using Gelfgren

poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)
xs = range(0, 1, length=100)
ys = evaluate(poly, collect(xs))

plot(xs, ys, label="f(x)", linewidth=2)
title!("Bernstein Polynomial")
xlabel!("x")
ylabel!("f(x)")
```

#### With DataFrames.jl

```julia
using DataFrames
using Gelfgren

poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)
xs = 0.0:0.1:1.0
ys = evaluate(poly, collect(xs))

df = DataFrame(x = collect(xs), y = ys)
```

#### Numerical Integration

```julia
using QuadGK

poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)

# Integrate using QuadGK
integral_value, error = quadgk(x -> evaluate(poly, x), 0.0, 1.0)

# Compare with built-in integral
ipoly = integral(poly)
y1 = evaluate(ipoly, 1.0)
y0 = evaluate(ipoly, 0.0)
area = y1 - y0
```

## API Reference

### Types

#### `BernsteinPolynomial`

Represents a Bernstein polynomial on interval [a, b].

**Constructor:**
```julia
BernsteinPolynomial(coeffs::AbstractVector{Float64}, a::Float64, b::Float64)
```

**Methods:**
- `evaluate(poly, x)` - Evaluate at point x
- `evaluate(poly, xs)` - Vectorized evaluation
- `derivative(poly)` - Compute derivative
- `integral(poly)` - Compute antiderivative
- `degree(poly)` - Get polynomial degree
- `+`, `-`, `*` - Arithmetic operations

#### `RationalFunction`

Represents rational function P(x)/Q(x).

**Constructor:**
```julia
RationalFunction(num::BernsteinPolynomial, den::BernsteinPolynomial)
```

**Methods:**
- `evaluate(rat, x)` - Evaluate at point x
- `evaluate(rat, xs)` - Vectorized evaluation
- `derivative(rat)` - Compute derivative using quotient rule

#### `PadeApproximant`

Represents Padé approximant [n/m].

**Constructor:**
```julia
PadeApproximant(coeffs::AbstractVector{Float64}, n::Int, m::Int, a::Float64, b::Float64)
```

**Methods:**
- `evaluate(pade, x)` - Evaluate at point x
- `evaluate(pade, xs)` - Vectorized evaluation

### Exception

#### `GelfgrenError`

Exception type for all Gelfgren errors.

```julia
struct GelfgrenError <: Exception
    msg::String
end
```

## Memory Management

The Julia bindings use finalizers for automatic memory management:

```julia
# Memory is automatically freed when object goes out of scope
function foo()
    poly = BernsteinPolynomial([1.0, 2.0], 0.0, 1.0)
    y = evaluate(poly, 0.5)
    # poly is automatically freed when function returns
end
```

Force garbage collection if needed:

```julia
GC.gc()
```

## Performance

The Julia bindings offer excellent performance:

- **Direct C calls**: Uses `ccall` for zero-overhead FFI
- **Type inference**: Julia's compiler optimizes type-stable code
- **Vectorization**: Efficient array operations
- **SIMD**: Julia can vectorize loops automatically

### Benchmarking

```julia
using BenchmarkTools

poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)

# Single evaluation
@btime evaluate($poly, 0.5)

# Vectorized evaluation
xs = rand(1000)
@btime evaluate($poly, $xs)
```

### Performance Tips

1. Use vectorized `evaluate(poly, xs)` instead of loops
2. Avoid global variables in performance-critical code
3. Use type annotations for better inference
4. Profile with `@profile` to find bottlenecks

## Testing

Run the test suite:

```julia
using Pkg
Pkg.test("Gelfgren")
```

Or from the command line:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## Examples

Complete examples are in the `examples/` directory:

```bash
julia examples/basic_usage.jl
```

## Troubleshooting

### Library Not Found

If you get "library not found" errors:

```julia
# Check library path
ENV["GELFGREN_LIB"] = "/full/path/to/libgelfgren.so"

# On Linux
ENV["LD_LIBRARY_PATH"] = "/path/to/gelfgren/target/release:" * get(ENV, "LD_LIBRARY_PATH", "")

# On macOS
ENV["DYLD_LIBRARY_PATH"] = "/path/to/gelfgren/target/release:" * get(ENV, "DYLD_LIBRARY_PATH", "")
```

### Julia Version

Check your Julia version:

```julia
versioninfo()
```

Requires Julia 1.9 or later.

### Precompilation Issues

Force recompilation:

```julia
using Pkg
Pkg.build("Gelfgren")
```

## Contributing

Contributions are welcome! See the main [CONTRIBUTING.md](../../CONTRIBUTING.md).

## License

Licensed under either of MIT or Apache-2.0 at your option.

## See Also

- [Main Gelfgren Documentation](../../README.md)
- [Julia Documentation](https://docs.julialang.org/)
- [Julia C Interface](https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/)
- [C API Header](../../include/gelfgren.h)
