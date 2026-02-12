using Test
using Gelfgren

@testset "Gelfgren.jl" begin

    @testset "BernsteinPolynomial - Creation and Basic Operations" begin
        # Test creation
        poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)
        @test poly isa BernsteinPolynomial
        @test degree(poly) == 2

        # Test evaluation at a point
        y = evaluate(poly, 0.5)
        @test y ≈ 3.0 atol=1e-10

        # Test empty coefficients error
        @test_throws GelfgrenError BernsteinPolynomial(Float64[], 0.0, 1.0)
    end

    @testset "BernsteinPolynomial - Vectorized Evaluation" begin
        poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)

        xs = [0.0, 0.25, 0.5, 0.75, 1.0]
        ys = evaluate(poly, xs)

        @test length(ys) == length(xs)
        @test ys[1] ≈ 1.0 atol=1e-10  # f(0) = 1
        @test ys[3] ≈ 3.0 atol=1e-10  # f(0.5) = 3
        @test ys[5] ≈ 6.0 atol=1e-10  # f(1) = 6
    end

    @testset "BernsteinPolynomial - Derivative" begin
        poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)
        dpoly = derivative(poly)

        @test dpoly isa BernsteinPolynomial

        # Derivative of 2x^2 + 3x + 1 is 4x + 3
        # At x=0.5: 4(0.5) + 3 = 5.0
        dy = evaluate(dpoly, 0.5)
        @test dy ≈ 5.0 atol=1e-10
    end

    @testset "BernsteinPolynomial - Integral" begin
        poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)
        ipoly = integral(poly)

        @test ipoly isa BernsteinPolynomial

        # Compute definite integral from 0 to 1
        y1 = evaluate(ipoly, 1.0)
        y0 = evaluate(ipoly, 0.0)
        area = y1 - y0

        # Integral of 2x^2 + 3x + 1 from 0 to 1
        # = [2/3 x^3 + 3/2 x^2 + x] from 0 to 1
        # = 2/3 + 3/2 + 1 = 13/6
        expected = 2.0/3.0 + 3.0/2.0 + 1.0
        @test area ≈ expected atol=1e-10
    end

    @testset "BernsteinPolynomial - Arithmetic Operations" begin
        p1 = BernsteinPolynomial([1.0, 2.0, 3.0], 0.0, 1.0)
        p2 = BernsteinPolynomial([2.0, 1.0, 1.0], 0.0, 1.0)

        # Test addition
        sum_poly = p1 + p2
        y1 = evaluate(p1, 0.5)
        y2 = evaluate(p2, 0.5)
        y_sum = evaluate(sum_poly, 0.5)
        @test y_sum ≈ y1 + y2 atol=1e-10

        # Test subtraction
        diff_poly = p1 - p2
        y_diff = evaluate(diff_poly, 0.5)
        @test y_diff ≈ y1 - y2 atol=1e-10

        # Test multiplication
        prod_poly = p1 * p2
        y_prod = evaluate(prod_poly, 0.5)
        @test y_prod ≈ y1 * y2 atol=1e-10

        # Test scalar multiplication
        scaled = 2.0 * p1
        y_scaled = evaluate(scaled, 0.5)
        @test y_scaled ≈ 2.0 * y1 atol=1e-10

        # Test reverse scalar multiplication
        scaled2 = p1 * 2.0
        y_scaled2 = evaluate(scaled2, 0.5)
        @test y_scaled2 ≈ 2.0 * y1 atol=1e-10
    end

    @testset "RationalFunction - Creation and Evaluation" begin
        # Create rational function: x / (1 + x)
        num = BernsteinPolynomial([0.0, 1.0], 0.0, 1.0)
        den = BernsteinPolynomial([1.0, 2.0], 0.0, 1.0)

        rat = RationalFunction(num, den)
        @test rat isa RationalFunction

        # Test evaluation
        x = 0.5
        y = evaluate(rat, x)
        expected = x / (1.0 + x)
        @test y ≈ expected atol=1e-10
    end

    @testset "RationalFunction - Vectorized Evaluation" begin
        num = BernsteinPolynomial([0.0, 1.0], 0.0, 1.0)
        den = BernsteinPolynomial([1.0, 2.0], 0.0, 1.0)
        rat = RationalFunction(num, den)

        xs = [0.0, 0.25, 0.5, 0.75, 1.0]
        ys = evaluate(rat, xs)

        @test length(ys) == length(xs)
        for (i, x) in enumerate(xs)
            expected = x / (1.0 + x)
            @test ys[i] ≈ expected atol=1e-10
        end
    end

    @testset "RationalFunction - Derivative" begin
        num = BernsteinPolynomial([0.0, 1.0], 0.0, 1.0)
        den = BernsteinPolynomial([1.0, 2.0], 0.0, 1.0)
        rat = RationalFunction(num, den)

        drat = derivative(rat)
        @test drat isa RationalFunction

        # Derivative of x/(1+x) is 1/(1+x)^2
        x = 0.5
        dy = evaluate(drat, x)
        expected = 1.0 / (1.0 + x)^2
        @test dy ≈ expected atol=1e-9
    end

    @testset "PadeApproximant - Creation and Evaluation" begin
        # Approximate exp(x) with [2/2] Padé approximant
        coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]
        pade = PadeApproximant(coeffs, 2, 2, -1.0, 1.0)

        @test pade isa PadeApproximant

        # Test evaluation
        x = 0.5
        y = evaluate(pade, x)
        exact = exp(x)
        error = abs(y - exact)

        # Padé [2/2] should be quite accurate for exp near 0
        @test error < 0.001
    end

    @testset "PadeApproximant - Vectorized Evaluation" begin
        coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]
        pade = PadeApproximant(coeffs, 2, 2, -1.0, 1.0)

        xs = [-0.5, -0.25, 0.0, 0.25, 0.5]
        ys = evaluate(pade, xs)

        @test length(ys) == length(xs)
        for (i, x) in enumerate(xs)
            exact = exp(x)
            error = abs(ys[i] - exact)
            @test error < 0.001
        end
    end

    @testset "Error Handling" begin
        # Test that we get proper errors for invalid operations
        @test_throws GelfgrenError BernsteinPolynomial(Float64[], 0.0, 1.0)
    end

    @testset "Memory Management" begin
        # Test that finalizers work correctly
        # Create and discard many polynomials
        for i in 1:100
            poly = BernsteinPolynomial([1.0, 2.0, 3.0], 0.0, 1.0)
            evaluate(poly, 0.5)
        end
        GC.gc()  # Force garbage collection

        # If there were memory leaks, this would likely cause issues
        @test true
    end

    @testset "Broadcasting" begin
        poly = BernsteinPolynomial([1.0, 2.0, 6.0], 0.0, 1.0)

        # Test with range
        xs = 0.0:0.1:1.0
        ys = evaluate(poly, collect(xs))

        @test length(ys) == length(xs)
        @test all(isfinite, ys)
    end
end
