# frozen_string_literal: true

require_relative "gelfgren/version"

begin
  # Try loading the native extension
  ruby_version = /(\d+\.\d+)/.match(RUBY_VERSION)
  require_relative "gelfgren/#{ruby_version}/gelfgren"
rescue LoadError
  require_relative "gelfgren/gelfgren"
end

module Gelfgren
  class Error < StandardError; end

  # Module methods for convenience
  class << self
    # Get library version
    def version
      VERSION
    end

    # Create a Bernstein polynomial from coefficients
    def bernstein(coeffs, a: 0.0, b: 1.0)
      BernsteinPolynomial.new(coeffs, a, b)
    end

    # Create a rational function from numerator and denominator polynomials
    def rational(numerator, denominator)
      RationalFunction.new(numerator, denominator)
    end

    # Create a Padé approximant from power series
    def pade_from_series(coeffs, n, m, a: -1.0, b: 1.0)
      PadeApproximant.from_power_series(coeffs, n, m, a, b)
    end

    # Create a Padé approximant from derivative values
    def pade_from_derivatives(derivs, n, m, x0: 0.0, a: -1.0, b: 1.0)
      PadeApproximant.from_derivatives(derivs, n, m, x0, a, b)
    end
  end
end
