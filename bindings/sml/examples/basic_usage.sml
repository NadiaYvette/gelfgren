(* Basic usage example for Gelfgren Standard ML bindings *)

fun main () =
    let
        open Gelfgren

        val _ = print "Gelfgren Standard ML Bindings Example\n"
        val _ = print "============================================================\n\n"

        (* Example 1: Bernstein Polynomial *)
        val _ = print "1. Bernstein Polynomial\n"
        val _ = print "------------------------------------------------------------\n"

        val poly = BernsteinPolynomial.create ([1.0, 2.0, 6.0], 0.0, 1.0)
        val deg = BernsteinPolynomial.degree poly
        val _ = print ("Created polynomial with degree: " ^ Int.toString deg ^ "\n")

        val y = BernsteinPolynomial.evaluate (poly, 0.5)
        val _ = print ("f(0.5) = " ^ Real.toString y ^ "\n")

        (* Vectorized evaluation *)
        val _ = print "\nVectorized evaluation:\n"
        val xs = [0.0, 0.25, 0.5, 0.75, 1.0]
        val ys = BernsteinPolynomial.evaluateList (poly, xs)
        val _ = ListPair.app
            (fn (x, y) => print ("  f(" ^ Real.toString x ^ ") = " ^ Real.toString y ^ "\n"))
            (xs, ys)

        (* Derivative *)
        val _ = print "\n"
        val dpoly = BernsteinPolynomial.derivative poly
        val dy = BernsteinPolynomial.evaluate (dpoly, 0.5)
        val _ = print ("Derivative at x=0.5: " ^ Real.toString dy ^ "\n")

        (* Integral *)
        val ipoly = BernsteinPolynomial.integral poly
        val y1 = BernsteinPolynomial.evaluate (ipoly, 1.0)
        val y0 = BernsteinPolynomial.evaluate (ipoly, 0.0)
        val area = y1 - y0
        val _ = print ("Integral from 0 to 1: " ^ Real.toString area ^ "\n")

        (* Clean up *)
        val _ = BernsteinPolynomial.free ipoly
        val _ = BernsteinPolynomial.free dpoly
        val _ = BernsteinPolynomial.free poly

        val _ = print "\n"

        (* Example 2: Rational Function *)
        val _ = print "2. Rational Function\n"
        val _ = print "------------------------------------------------------------\n"

        val num = BernsteinPolynomial.create ([0.0, 1.0], 0.0, 1.0)
        val den = BernsteinPolynomial.create ([1.0, 2.0], 0.0, 1.0)
        val rat = RationalFunction.create (num, den)

        val x = 0.5
        val y_rat = RationalFunction.evaluate (rat, x)
        val expected = x / (1.0 + x)
        val _ = print ("R(0.5) = " ^ Real.toString y_rat ^ "\n")
        val _ = print ("Expected: " ^ Real.toString expected ^ "\n")

        (* Vectorized *)
        val _ = print "\nVectorized evaluation:\n"
        val xs_rat = [0.0, 0.25, 0.5, 0.75, 1.0]
        val ys_rat = RationalFunction.evaluateList (rat, xs_rat)
        val _ = ListPair.app
            (fn (x, y) =>
                let val exp = x / (1.0 + x)
                in print ("  R(" ^ Real.toString x ^ ") = " ^ Real.toString y ^
                          " (expected: " ^ Real.toString exp ^ ")\n")
                end)
            (xs_rat, ys_rat)

        (* Clean up *)
        val _ = RationalFunction.free rat
        val _ = BernsteinPolynomial.free den
        val _ = BernsteinPolynomial.free num

        val _ = print "\n"

        (* Example 3: Padé Approximant *)
        val _ = print "3. Padé Approximant\n"
        val _ = print "------------------------------------------------------------\n"

        val coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]
        val pade = PadeApproximant.fromSeries (coeffs, 2, 2, ~1.0, 1.0)

        val _ = print "Approximation of exp(x):\n"
        val xs_pade = [~0.5, ~0.25, 0.0, 0.25, 0.5]
        val ys_pade = PadeApproximant.evaluateList (pade, xs_pade)
        val _ = ListPair.app
            (fn (x, y) =>
                let
                    val exact = Math.exp x
                    val error = Real.abs (y - exact)
                in
                    print ("  x=" ^ Real.toString x ^ ": Padé=" ^ Real.toString y ^
                           ", exp(x)=" ^ Real.toString exact ^
                           ", error=" ^ Real.toString error ^ "\n")
                end)
            (xs_pade, ys_pade)

        val _ = PadeApproximant.free pade

        val _ = print "\n"
        val _ = print "All examples completed successfully!\n"
    in
        ()
    end

val _ = main ()
