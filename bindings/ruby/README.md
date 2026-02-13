# Gelfgren Ruby Bindings

Ruby bindings for the Gelfgren numerical computing library, providing high-performance piecewise rational interpolation methods.

## Installation

Add this line to your application's Gemfile:

```ruby
gem 'gelfgren'
```

And then execute:

```bash
bundle install
```

Or install it yourself as:

```bash
gem install gelfgren
```

## Requirements

- Ruby 3.0 or later
- Rust toolchain (for building from source)

## Usage

### Basic Example

```ruby
require 'gelfgren'

# Create a Bernstein polynomial
poly = Gelfgren::BernsteinPolynomial.new([1.0, 2.0, 6.0], 0.0, 1.0)

# Evaluate at a point
y = poly.evaluate(0.5)
puts "f(0.5) = #{y}"

# Vectorized evaluation
x_array = [0.0, 0.25, 0.5, 0.75, 1.0]
y_array = poly.evaluate(x_array)

# Compute derivative
dpoly = poly.derivative
```

### Rational Functions

```ruby
# Create rational function: x / (1 + x)
num = Gelfgren::BernsteinPolynomial.new([0.0, 1.0], 0.0, 1.0)
den = Gelfgren::BernsteinPolynomial.new([1.0, 2.0], 0.0, 1.0)
rat = Gelfgren::RationalFunction.new(num, den)

# Evaluate
y = rat.evaluate(0.5)

# Derivative
drat = rat.derivative
```

### Padé Approximants

```ruby
# Approximate exp(x) with [2/2] Padé approximant
coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]
pade = Gelfgren::PadeApproximant.from_power_series(coeffs, 2, 2, -1.0, 1.0)

# Evaluate
x = 0.5
approx = pade.evaluate(x)
exact = Math.exp(x)
error = (approx - exact).abs
puts "Error: #{error}"
```

### Polynomial Arithmetic

```ruby
p1 = Gelfgren::BernsteinPolynomial.new([1.0, 2.0, 3.0], 0.0, 1.0)
p2 = Gelfgren::BernsteinPolynomial.new([2.0, 1.0, 1.0], 0.0, 1.0)

# Addition
p_sum = p1 + p2

# Subtraction
p_diff = p1 - p2

# Multiplication
p_prod = p1 * p2

# Scaling
p_scaled = p1.scale(2.0)
```

### Convenience Methods

```ruby
# Use module-level methods for concise syntax
poly = Gelfgren.bernstein([1.0, 2.0, 6.0], a: 0.0, b: 1.0)

num = Gelfgren.bernstein([0.0, 1.0], a: 0.0, b: 1.0)
den = Gelfgren.bernstein([1.0, 2.0], a: 0.0, b: 1.0)
rat = Gelfgren.rational(num, den)

pade = Gelfgren.pade_from_series([1.0, 1.0, 0.5], 1, 1)
```

## API Reference

### `Gelfgren::BernsteinPolynomial`

#### Constructor

```ruby
BernsteinPolynomial.new(coeffs, a, b)
```

- `coeffs` (Array<Float>): Bernstein coefficients
- `a` (Float): Left endpoint of interval
- `b` (Float): Right endpoint of interval

#### Methods

- `evaluate(x)`: Evaluate at point(s). Accepts Float or Array<Float>
- `derivative()`: Compute derivative (returns new BernsteinPolynomial)
- `integral()`: Compute antiderivative (returns new BernsteinPolynomial)
- `degree()`: Get polynomial degree (Integer)
- `interval()`: Get interval endpoints (Array<Float>)
- `+(other)`: Add polynomials
- `-(other)`: Subtract polynomials
- `*(other)`: Multiply polynomials
- `scale(scalar)`: Multiply by scalar
- `elevate_degree(r)`: Elevate degree by r
- `to_s()`: String representation

### `Gelfgren::RationalFunction`

#### Constructor

```ruby
RationalFunction.new(numerator, denominator)
```

- `numerator` (BernsteinPolynomial): Numerator polynomial
- `denominator` (BernsteinPolynomial): Denominator polynomial

#### Methods

- `evaluate(x)`: Evaluate at point(s)
- `derivative()`: Compute derivative
- `numerator()`: Get numerator polynomial
- `denominator()`: Get denominator polynomial
- `+(other)`: Add rational functions
- `-(other)`: Subtract rational functions
- `*(other)`: Multiply rational functions
- `/(other)`: Divide rational functions
- `reciprocal()`: Compute reciprocal
- `to_s()`: String representation

### `Gelfgren::PadeApproximant`

#### Constructors

```ruby
PadeApproximant.from_power_series(coeffs, n, m, a, b)
```

- `coeffs` (Array<Float>): Power series coefficients
- `n` (Integer): Numerator degree
- `m` (Integer): Denominator degree
- `a`, `b` (Float): Interval endpoints

```ruby
PadeApproximant.from_derivatives(derivs, n, m, x0, a, b)
```

- `derivs` (Array<Float>): Derivative values at x0
- `n`, `m` (Integer): Numerator and denominator degrees
- `x0` (Float): Expansion point
- `a`, `b` (Float): Interval endpoints

#### Methods

- `evaluate(x)`: Evaluate at point(s)
- `rational()`: Get underlying RationalFunction
- `orders()`: Get [n, m] array
- `to_s()`: String representation

### Module Methods

- `Gelfgren.version`: Get version string
- `Gelfgren.bernstein(coeffs, a:, b:)`: Create BernsteinPolynomial
- `Gelfgren.rational(num, den)`: Create RationalFunction
- `Gelfgren.pade_from_series(coeffs, n, m, a:, b:)`: Create PadeApproximant from series
- `Gelfgren.pade_from_derivatives(derivs, n, m, x0:, a:, b:)`: Create PadeApproximant from derivatives

## Development

After checking out the repo, run:

```bash
bundle install
rake compile
```

Run tests:

```bash
rake spec
```

Run examples:

```bash
rake examples
```

Build gem:

```bash
rake build
```

Install locally:

```bash
rake install_local
```

## Performance

The Ruby bindings use the high-performance Rust core via FFI with Magnus, providing:

- Native performance for numerical operations
- Minimal overhead for vectorized operations
- Automatic memory management

Benchmarks show 5-10x speedup compared to pure Ruby implementations for polynomial evaluation.

## Contributing

Bug reports and pull requests are welcome on GitHub at https://github.com/yourusername/gelfgren.

## License

The gem is available as open source under the terms of the MIT License or Apache License 2.0.

## See Also

- [Main Gelfgren Documentation](https://yourusername.github.io/gelfgren)
- [Rust API Documentation](https://docs.rs/gelfgren-core)
- [Mathematical Background](https://yourusername.github.io/gelfgren/mathematics.html)
