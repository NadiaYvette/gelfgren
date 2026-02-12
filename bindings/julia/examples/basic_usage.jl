#!/usr/bin/env julia

# Basic usage example for Gelfgren Julia bindings

using Gelfgren

println("Gelfgren Julia Bindings Example")
println("=" ^ 60)
println()

# Example 1: Bernstein Polynomial
println("1. Bernstein Polynomial")
println("-" ^ 60)

# Create a quadratic polynomial: 2x^2 + 3x + 1 on [0, 1]
poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)

# Get degree
deg = degree(poly)
println("Created polynomial with degree: $deg")

# Evaluate at a point
y = evaluate(poly, 0.5)
println("f(0.5) = $y")

# Vectorized evaluation
println("\nVectorized evaluation:")
xs = [0.0, 0.25, 0.5, 0.75, 1.0]
ys = evaluate(poly, xs)
for (x, y) in zip(xs, ys)
    println("  f($x) = $y")
end

# Derivative
println()
dpoly = derivative(poly)
dy = evaluate(dpoly, 0.5)
println("Derivative at x=0.5: $dy")

# Integral
ipoly = integral(poly)
y1 = evaluate(ipoly, 1.0)
y0 = evaluate(ipoly, 0.0)
area = y1 - y0
println("Integral from 0 to 1: $area")

println()

# Example 2: Polynomial Arithmetic
println("2. Polynomial Arithmetic")
println("-" ^ 60)

p1 = BernsteinPolynomial([1.0, 2.0, 3.0], 0.0, 1.0)
p2 = BernsteinPolynomial([2.0, 1.0, 1.0], 0.0, 1.0)

y1 = evaluate(p1, 0.5)
y2 = evaluate(p2, 0.5)
println("p1(0.5) = $y1")
println("p2(0.5) = $y2")

# Addition
sum_poly = p1 + p2
y_sum = evaluate(sum_poly, 0.5)
println("p1 + p2 at x=0.5: $y_sum")

# Scaling
scaled = 2.0 * p1
y_scaled = evaluate(scaled, 0.5)
println("2 * p1 at x=0.5: $y_scaled")

println()

# Example 3: Rational Function
println("3. Rational Function")
println("-" ^ 60)

# Create rational function: x / (1 + x) on [0, 1]
num = BernsteinPolynomial([0.0, 1.0], 0.0, 1.0)
den = BernsteinPolynomial([1.0, 2.0], 0.0, 1.0)

rat = RationalFunction(num, den)

x = 0.5
y = evaluate(rat, x)
expected = x / (1.0 + x)
println("R($x) = $y")
println("Expected: $expected")

# Vectorized evaluation
println("\nVectorized evaluation:")
xs = [0.0, 0.25, 0.5, 0.75, 1.0]
ys = evaluate(rat, xs)
for (x, y) in zip(xs, ys)
    expected = x / (1.0 + x)
    println("  R($x) = $y (expected: $expected)")
end

println()

# Example 4: Padé Approximant
println("4. Padé Approximant")
println("-" ^ 60)

# Approximate exp(x) near x=0 with [2/2] Padé approximant
# Power series: 1 + x + x^2/2 + x^3/6 + x^4/24
coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]
pade = PadeApproximant(coeffs, 2, 2, -1.0, 1.0)

println("Approximation of exp(x):")
xs = [-0.5, -0.25, 0.0, 0.25, 0.5]
ys = evaluate(pade, xs)

for (x, y) in zip(xs, ys)
    exact = exp(x)
    error = abs(y - exact)
    @printf("  x=%.2f: Padé=%.6f, exp(x)=%.6f, error=%.8f\n", x, y, exact, error)
end

println()

# Example 5: Broadcasting and Ranges
println("5. Broadcasting with Ranges")
println("-" ^ 60)

poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)

# Evaluate over a range
xs = 0.0:0.1:1.0
ys = evaluate(poly, collect(xs))

println("Evaluated at $(length(ys)) points")
println("  min: $(minimum(ys))")
println("  max: $(maximum(ys))")
println("  mean: $(sum(ys) / length(ys))")

println()

# Example 6: Error Handling
println("6. Error Handling")
println("-" ^ 60)

try
    # This should fail - empty coefficients
    bad_poly = BernsteinPolynomial(Float64[], 0.0, 1.0)
catch e
    if e isa GelfgrenError
        println("Caught expected error: $(e.msg)")
    else
        rethrow(e)
    end
end

println()
println("All examples completed successfully!")
