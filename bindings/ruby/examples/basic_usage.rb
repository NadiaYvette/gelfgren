#!/usr/bin/env ruby
# frozen_string_literal: true

require "bundler/setup"
require "gelfgren"

puts "Gelfgren Ruby Bindings Example"
puts "=" * 40

# Example 1: Bernstein Polynomial
puts "\n1. Bernstein Polynomial"
puts "-" * 40

# Create a quadratic polynomial: 2x^2 + 3x + 1 on [0, 1]
# Bernstein coefficients: [1, 2, 6]
poly = Gelfgren::BernsteinPolynomial.new([1.0, 2.0, 6.0], 0.0, 1.0)

puts "Created: #{poly}"
puts "Degree: #{poly.degree}"
puts "Interval: #{poly.interval.inspect}"

# Evaluate at a point
x = 0.5
y = poly.evaluate(x)
puts "f(#{x}) = #{y}"

# Vectorized evaluation
x_array = [0.0, 0.25, 0.5, 0.75, 1.0]
y_array = poly.evaluate(x_array)
puts "\nVectorized evaluation:"
x_array.zip(y_array).each do |xi, yi|
  puts "  f(#{xi}) = #{yi}"
end

# Derivative
dpoly = poly.derivative
puts "\nDerivative: #{dpoly}"
dy = dpoly.evaluate(0.5)
puts "f'(0.5) = #{dy}"

# Integral
integral = poly.integral
puts "\nIntegral: #{integral}"
area = integral.evaluate(1.0) - integral.evaluate(0.0)
puts "∫f(x)dx from 0 to 1 = #{area}"

# Example 2: Rational Function
puts "\n2. Rational Function"
puts "-" * 40

# Create rational function: x / (1 + x) on [0, 1]
num = Gelfgren::BernsteinPolynomial.new([0.0, 1.0], 0.0, 1.0)
den = Gelfgren::BernsteinPolynomial.new([1.0, 2.0], 0.0, 1.0)
rat = Gelfgren::RationalFunction.new(num, den)

puts "Created: #{rat}"

# Evaluate
x = 0.5
y = rat.evaluate(x)
puts "R(#{x}) = #{y}"
puts "Expected: #{x / (1 + x)}"

# Vectorized evaluation
x_array = [0.0, 0.25, 0.5, 0.75, 1.0]
y_array = rat.evaluate(x_array)
puts "\nVectorized evaluation:"
x_array.zip(y_array).each do |xi, yi|
  expected = xi / (1 + xi)
  puts "  R(#{xi}) = #{yi} (expected: #{expected})"
end

# Derivative
drat = rat.derivative
puts "\nDerivative: #{drat}"
dy = drat.evaluate(0.5)
puts "R'(0.5) = #{dy}"

# Example 3: Padé Approximant
puts "\n3. Padé Approximant"
puts "-" * 40

# Approximate exp(x) near x=0 with [2/2] Padé approximant
# Power series: 1 + x + x^2/2 + x^3/6 + x^4/24
coeffs = [1.0, 1.0, 0.5, 1.0 / 6.0, 1.0 / 24.0]
pade = Gelfgren::PadeApproximant.from_power_series(coeffs, 2, 2, -1.0, 1.0)

puts "Created: #{pade}"
puts "Orders: #{pade.orders.inspect}"

# Evaluate and compare with exact
x_array = [-0.5, -0.25, 0.0, 0.25, 0.5]
puts "\nApproximation of exp(x):"
x_array.each do |x|
  approx = pade.evaluate(x)
  exact = Math.exp(x)
  error = (approx - exact).abs
  puts "  x=#{x}: Padé=#{approx.round(6)}, exp(x)=#{exact.round(6)}, error=#{error.round(8)}"
end

# Example 4: Polynomial Arithmetic
puts "\n4. Polynomial Arithmetic"
puts "-" * 40

p1 = Gelfgren::BernsteinPolynomial.new([1.0, 2.0, 3.0], 0.0, 1.0)
p2 = Gelfgren::BernsteinPolynomial.new([2.0, 1.0, 1.0], 0.0, 1.0)

puts "p1: #{p1}"
puts "p2: #{p2}"

# Addition
p_sum = p1 + p2
puts "\np1 + p2: #{p_sum}"
puts "At x=0.5: #{p_sum.evaluate(0.5)}"

# Subtraction
p_diff = p1 - p2
puts "\np1 - p2: #{p_diff}"
puts "At x=0.5: #{p_diff.evaluate(0.5)}"

# Multiplication
p_prod = p1 * p2
puts "\np1 * p2: #{p_prod}"
puts "Degree: #{p_prod.degree}"
puts "At x=0.5: #{p_prod.evaluate(0.5)}"

# Scaling
p_scaled = p1.scale(2.0)
puts "\n2 * p1: #{p_scaled}"
puts "At x=0.5: #{p_scaled.evaluate(0.5)}"

# Example 5: Convenience Methods
puts "\n5. Convenience Methods"
puts "-" * 40

# Using module methods
poly = Gelfgren.bernstein([1.0, 2.0, 6.0], a: 0.0, b: 1.0)
puts "Created via Gelfgren.bernstein: #{poly}"

num = Gelfgren.bernstein([0.0, 1.0], a: 0.0, b: 1.0)
den = Gelfgren.bernstein([1.0, 2.0], a: 0.0, b: 1.0)
rat = Gelfgren.rational(num, den)
puts "Created via Gelfgren.rational: #{rat}"

pade = Gelfgren.pade_from_series([1.0, 1.0, 0.5, 1.0 / 6.0], 2, 1)
puts "Created via Gelfgren.pade_from_series: #{pade}"

puts "\nAll examples completed successfully!"
