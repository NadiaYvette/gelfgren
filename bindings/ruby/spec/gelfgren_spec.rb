# frozen_string_literal: true

require "spec_helper"

RSpec.describe Gelfgren do
  it "has a version number" do
    expect(Gelfgren::VERSION).not_to be nil
  end

  describe ".version" do
    it "returns the version string" do
      expect(Gelfgren.version).to eq(Gelfgren::VERSION)
    end
  end

  describe "BernsteinPolynomial" do
    describe ".new" do
      it "creates a polynomial from coefficients" do
        poly = Gelfgren::BernsteinPolynomial.new([1.0, 2.0, 3.0], 0.0, 1.0)
        expect(poly).to be_a(Gelfgren::BernsteinPolynomial)
      end

      it "raises error for empty coefficients" do
        expect {
          Gelfgren::BernsteinPolynomial.new([], 0.0, 1.0)
        }.to raise_error(RuntimeError)
      end

      it "raises error for invalid interval" do
        expect {
          Gelfgren::BernsteinPolynomial.new([1.0, 2.0], 1.0, 0.0)
        }.to raise_error(RuntimeError)
      end
    end

    describe "#degree" do
      it "returns the polynomial degree" do
        poly = Gelfgren::BernsteinPolynomial.new([1.0, 2.0, 3.0], 0.0, 1.0)
        expect(poly.degree).to eq(2)
      end
    end

    describe "#interval" do
      it "returns the interval endpoints" do
        poly = Gelfgren::BernsteinPolynomial.new([1.0, 2.0], 0.0, 1.0)
        expect(poly.interval).to eq([0.0, 1.0])
      end
    end

    describe "#evaluate" do
      it "evaluates at a single point" do
        # Constant polynomial: f(x) = 2
        poly = Gelfgren::BernsteinPolynomial.new([2.0, 2.0], 0.0, 1.0)
        expect(poly.evaluate(0.5)).to be_close_to(2.0)
      end

      it "evaluates at left endpoint" do
        poly = Gelfgren::BernsteinPolynomial.new([1.0, 2.0, 3.0], 0.0, 1.0)
        expect(poly.evaluate(0.0)).to be_close_to(1.0)
      end

      it "evaluates at right endpoint" do
        poly = Gelfgren::BernsteinPolynomial.new([1.0, 2.0, 3.0], 0.0, 1.0)
        expect(poly.evaluate(1.0)).to be_close_to(3.0)
      end

      it "evaluates at multiple points (vectorized)" do
        poly = Gelfgren::BernsteinPolynomial.new([1.0, 2.0, 3.0], 0.0, 1.0)
        results = poly.evaluate([0.0, 0.5, 1.0])
        expect(results).to be_an(Array)
        expect(results.length).to eq(3)
        expect(results[0]).to be_close_to(1.0)
        expect(results[2]).to be_close_to(3.0)
      end
    end

    describe "#derivative" do
      it "computes the derivative" do
        # Linear: f(x) = 2x + 1, f'(x) = 2
        poly = Gelfgren::BernsteinPolynomial.new([1.0, 3.0], 0.0, 1.0)
        dpoly = poly.derivative
        expect(dpoly.degree).to eq(0)
        expect(dpoly.evaluate(0.5)).to be_close_to(2.0, tolerance: 1e-9)
      end
    end

    describe "#integral" do
      it "computes the integral" do
        # Constant: f(x) = 2, ∫f = 2x
        poly = Gelfgren::BernsteinPolynomial.new([2.0, 2.0], 0.0, 1.0)
        integral = poly.integral
        area = integral.evaluate(1.0) - integral.evaluate(0.0)
        expect(area).to be_close_to(2.0, tolerance: 1e-9)
      end
    end

    describe "#+" do
      it "adds two polynomials" do
        p1 = Gelfgren::BernsteinPolynomial.new([1.0, 2.0], 0.0, 1.0)
        p2 = Gelfgren::BernsteinPolynomial.new([2.0, 1.0], 0.0, 1.0)
        p_sum = p1 + p2
        expect(p_sum.evaluate(0.0)).to be_close_to(3.0)
        expect(p_sum.evaluate(1.0)).to be_close_to(3.0)
      end
    end

    describe "#-" do
      it "subtracts two polynomials" do
        p1 = Gelfgren::BernsteinPolynomial.new([3.0, 3.0], 0.0, 1.0)
        p2 = Gelfgren::BernsteinPolynomial.new([1.0, 1.0], 0.0, 1.0)
        p_diff = p1 - p2
        expect(p_diff.evaluate(0.5)).to be_close_to(2.0)
      end
    end

    describe "#*" do
      it "multiplies two polynomials" do
        # (x)(x) = x^2
        p1 = Gelfgren::BernsteinPolynomial.new([0.0, 1.0], 0.0, 1.0)
        p2 = Gelfgren::BernsteinPolynomial.new([0.0, 1.0], 0.0, 1.0)
        p_prod = p1 * p2
        expect(p_prod.degree).to eq(2)
        expect(p_prod.evaluate(0.5)).to be_close_to(0.25, tolerance: 1e-9)
      end
    end

    describe "#scale" do
      it "scales a polynomial by a constant" do
        poly = Gelfgren::BernsteinPolynomial.new([1.0, 2.0], 0.0, 1.0)
        scaled = poly.scale(3.0)
        expect(scaled.evaluate(0.0)).to be_close_to(3.0)
        expect(scaled.evaluate(1.0)).to be_close_to(6.0)
      end
    end

    describe "#elevate_degree" do
      it "elevates the degree without changing the function" do
        poly = Gelfgren::BernsteinPolynomial.new([1.0, 2.0], 0.0, 1.0)
        elevated = poly.elevate_degree(2)
        expect(elevated.degree).to eq(3)
        expect(elevated.evaluate(0.5)).to be_close_to(poly.evaluate(0.5), tolerance: 1e-9)
      end
    end

    describe "#to_s" do
      it "returns a string representation" do
        poly = Gelfgren::BernsteinPolynomial.new([1.0, 2.0], 0.0, 1.0)
        expect(poly.to_s).to include("BernsteinPolynomial")
        expect(poly.to_s).to include("degree=1")
      end
    end
  end

  describe "RationalFunction" do
    describe ".new" do
      it "creates a rational function" do
        num = Gelfgren::BernsteinPolynomial.new([0.0, 1.0], 0.0, 1.0)
        den = Gelfgren::BernsteinPolynomial.new([1.0, 1.0], 0.0, 1.0)
        rat = Gelfgren::RationalFunction.new(num, den)
        expect(rat).to be_a(Gelfgren::RationalFunction)
      end
    end

    describe "#evaluate" do
      it "evaluates the rational function" do
        # x / 1 = x
        num = Gelfgren::BernsteinPolynomial.new([0.0, 1.0], 0.0, 1.0)
        den = Gelfgren::BernsteinPolynomial.new([1.0, 1.0], 0.0, 1.0)
        rat = Gelfgren::RationalFunction.new(num, den)
        expect(rat.evaluate(0.5)).to be_close_to(0.5, tolerance: 1e-9)
      end

      it "evaluates at multiple points" do
        num = Gelfgren::BernsteinPolynomial.new([0.0, 1.0], 0.0, 1.0)
        den = Gelfgren::BernsteinPolynomial.new([1.0, 2.0], 0.0, 1.0)
        rat = Gelfgren::RationalFunction.new(num, den)
        results = rat.evaluate([0.0, 0.5, 1.0])
        expect(results[0]).to be_close_to(0.0, tolerance: 1e-9)
        expect(results[2]).to be_close_to(0.5, tolerance: 1e-9)
      end
    end

    describe "#numerator and #denominator" do
      it "returns the numerator and denominator" do
        num = Gelfgren::BernsteinPolynomial.new([0.0, 1.0], 0.0, 1.0)
        den = Gelfgren::BernsteinPolynomial.new([1.0, 1.0], 0.0, 1.0)
        rat = Gelfgren::RationalFunction.new(num, den)

        expect(rat.numerator.evaluate(0.5)).to be_close_to(num.evaluate(0.5), tolerance: 1e-9)
        expect(rat.denominator.evaluate(0.5)).to be_close_to(den.evaluate(0.5), tolerance: 1e-9)
      end
    end

    describe "#derivative" do
      it "computes the derivative" do
        num = Gelfgren::BernsteinPolynomial.new([0.0, 1.0], 0.0, 1.0)
        den = Gelfgren::BernsteinPolynomial.new([1.0, 1.0], 0.0, 1.0)
        rat = Gelfgren::RationalFunction.new(num, den)
        drat = rat.derivative
        expect(drat).to be_a(Gelfgren::RationalFunction)
      end
    end

    describe "arithmetic operations" do
      let(:r1) do
        num = Gelfgren::BernsteinPolynomial.new([1.0, 1.0], 0.0, 1.0)
        den = Gelfgren::BernsteinPolynomial.new([1.0, 1.0], 0.0, 1.0)
        Gelfgren::RationalFunction.new(num, den)
      end

      let(:r2) do
        num = Gelfgren::BernsteinPolynomial.new([2.0, 2.0], 0.0, 1.0)
        den = Gelfgren::BernsteinPolynomial.new([1.0, 1.0], 0.0, 1.0)
        Gelfgren::RationalFunction.new(num, den)
      end

      it "adds rational functions" do
        r_sum = r1 + r2
        expect(r_sum.evaluate(0.5)).to be_close_to(3.0, tolerance: 1e-9)
      end

      it "subtracts rational functions" do
        r_diff = r2 - r1
        expect(r_diff.evaluate(0.5)).to be_close_to(1.0, tolerance: 1e-9)
      end

      it "multiplies rational functions" do
        r_prod = r1 * r2
        expect(r_prod.evaluate(0.5)).to be_close_to(2.0, tolerance: 1e-9)
      end

      it "divides rational functions" do
        r_quot = r2 / r1
        expect(r_quot.evaluate(0.5)).to be_close_to(2.0, tolerance: 1e-9)
      end
    end

    describe "#reciprocal" do
      it "computes the reciprocal" do
        num = Gelfgren::BernsteinPolynomial.new([2.0, 2.0], 0.0, 1.0)
        den = Gelfgren::BernsteinPolynomial.new([1.0, 1.0], 0.0, 1.0)
        rat = Gelfgren::RationalFunction.new(num, den)
        recip = rat.reciprocal
        expect(recip.evaluate(0.5)).to be_close_to(0.5, tolerance: 1e-9)
      end
    end
  end

  describe "PadeApproximant" do
    describe ".from_power_series" do
      it "creates a Padé approximant from power series" do
        coeffs = [1.0, 1.0, 0.5, 1.0 / 6.0, 1.0 / 24.0]
        pade = Gelfgren::PadeApproximant.from_power_series(coeffs, 2, 2, -1.0, 1.0)
        expect(pade).to be_a(Gelfgren::PadeApproximant)
      end

      it "approximates exp(x) accurately" do
        coeffs = [1.0, 1.0, 0.5, 1.0 / 6.0, 1.0 / 24.0]
        pade = Gelfgren::PadeApproximant.from_power_series(coeffs, 2, 2, -1.0, 1.0)

        x = 0.5
        approx = pade.evaluate(x)
        exact = Math.exp(x)
        expect(approx).to be_close_to(exact, tolerance: 1e-5)
      end
    end

    describe ".from_derivatives" do
      it "creates a Padé approximant from derivatives" do
        derivs = [1.0, 1.0, 1.0, 1.0, 1.0]  # exp(0) and derivatives
        pade = Gelfgren::PadeApproximant.from_derivatives(derivs, 2, 2, 0.0, -1.0, 1.0)
        expect(pade).to be_a(Gelfgren::PadeApproximant)
      end
    end

    describe "#orders" do
      it "returns the numerator and denominator orders" do
        coeffs = [1.0, 1.0, 0.5]
        pade = Gelfgren::PadeApproximant.from_power_series(coeffs, 1, 1, -1.0, 1.0)
        expect(pade.orders).to eq([1, 1])
      end
    end

    describe "#rational" do
      it "returns the underlying rational function" do
        coeffs = [1.0, 1.0, 0.5]
        pade = Gelfgren::PadeApproximant.from_power_series(coeffs, 1, 1, -1.0, 1.0)
        rat = pade.rational
        expect(rat).to be_a(Gelfgren::RationalFunction)
      end
    end

    describe "#evaluate" do
      it "evaluates at multiple points" do
        coeffs = [1.0, 1.0, 0.5, 1.0 / 6.0]
        pade = Gelfgren::PadeApproximant.from_power_series(coeffs, 2, 1, -1.0, 1.0)
        results = pade.evaluate([-0.5, 0.0, 0.5])
        expect(results).to be_an(Array)
        expect(results.length).to eq(3)
        expect(results[1]).to be_close_to(1.0, tolerance: 1e-9)
      end
    end
  end

  describe "Convenience methods" do
    describe ".bernstein" do
      it "creates a Bernstein polynomial" do
        poly = Gelfgren.bernstein([1.0, 2.0, 3.0], a: 0.0, b: 1.0)
        expect(poly).to be_a(Gelfgren::BernsteinPolynomial)
        expect(poly.degree).to eq(2)
      end
    end

    describe ".rational" do
      it "creates a rational function" do
        num = Gelfgren.bernstein([1.0, 2.0], a: 0.0, b: 1.0)
        den = Gelfgren.bernstein([1.0, 1.0], a: 0.0, b: 1.0)
        rat = Gelfgren.rational(num, den)
        expect(rat).to be_a(Gelfgren::RationalFunction)
      end
    end

    describe ".pade_from_series" do
      it "creates a Padé approximant" do
        pade = Gelfgren.pade_from_series([1.0, 1.0, 0.5], 1, 1)
        expect(pade).to be_a(Gelfgren::PadeApproximant)
      end
    end

    describe ".pade_from_derivatives" do
      it "creates a Padé approximant from derivatives" do
        pade = Gelfgren.pade_from_derivatives([1.0, 1.0, 1.0], 1, 1, x0: 0.0)
        expect(pade).to be_a(Gelfgren::PadeApproximant)
      end
    end
  end
end
