package org.gelfgren;

import java.util.Arrays;

/**
 * Example usage of Gelfgren Java bindings.
 */
public class Example {
    public static void main(String[] args) {
        System.out.println("Gelfgren Java Bindings Example");
        System.out.println("==============================\n");

        demonstratePolynomials();
        demonstrateRationalFunctions();
        demonstrateExceptionHandling();

        System.out.println("\nAll examples completed successfully!");
    }

    private static void demonstratePolynomials() {
        System.out.println("Bernstein Polynomial Operations");
        System.out.println("================================\n");

        // Create polynomial P(x) = 1 + 2x + 3x²
        try (BernsteinPolynomial p = new BernsteinPolynomial(
                new double[]{1.0, 2.0, 3.0}, 0.0, 1.0)) {

            System.out.println("Created: " + p);
            System.out.println("Degree: " + p.getDegree());
            double[] interval = p.getInterval();
            System.out.printf("Interval: [%.1f, %.1f]\n\n", interval[0], interval[1]);

            // Evaluate at several points
            System.out.println("Evaluations:");
            double[] x = {0.0, 0.25, 0.5, 0.75, 1.0};
            for (double xi : x) {
                System.out.printf("  P(%.2f) = %.6f\n", xi, p.evaluate(xi));
            }
            System.out.println();

            // Compute derivative
            System.out.println("Derivative P'(x):");
            try (BernsteinPolynomial pPrime = p.derivative()) {
                for (double xi : x) {
                    System.out.printf("  P'(%.2f) = %.6f\n", xi, pPrime.evaluate(xi));
                }
            }
            System.out.println();

            // Polynomial arithmetic
            System.out.println("Polynomial Arithmetic:");
            try (BernsteinPolynomial q = new BernsteinPolynomial(
                    new double[]{1.0, 0.0}, 0.0, 1.0);
                 BernsteinPolynomial sum = p.add(q);
                 BernsteinPolynomial diff = p.subtract(q);
                 BernsteinPolynomial prod = p.multiply(q)) {

                double xTest = 0.5;
                System.out.printf("  At x = %.1f:\n", xTest);
                System.out.printf("    P(x) + Q(x) = %.6f\n", sum.evaluate(xTest));
                System.out.printf("    P(x) - Q(x) = %.6f\n", diff.evaluate(xTest));
                System.out.printf("    P(x) * Q(x) = %.6f\n", prod.evaluate(xTest));
            }
        }
        System.out.println();
    }

    private static void demonstrateRationalFunctions() {
        System.out.println("Rational Functions");
        System.out.println("==================\n");

        // Create R(x) = (1 + x) / (1 + 2x)
        try (BernsteinPolynomial num = new BernsteinPolynomial(
                new double[]{1.0, 1.0}, 0.0, 1.0);
             BernsteinPolynomial den = new BernsteinPolynomial(
                new double[]{1.0, 2.0}, 0.0, 1.0);
             RationalFunction r = new RationalFunction(num, den)) {

            System.out.println("Created rational function R(x) = (1 + x) / (1 + 2x)\n");

            System.out.println("Evaluations:");
            double[] x = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
            for (double xi : x) {
                System.out.printf("  R(%.1f) = %.6f\n", xi, r.evaluate(xi));
            }
            System.out.println();

            // Derivative
            System.out.println("Derivative R'(x):");
            try (RationalFunction rPrime = r.derivative()) {
                for (double xi : x) {
                    System.out.printf("  R'(%.1f) = %.6f\n", xi, rPrime.evaluate(xi));
                }
            }
        }
        System.out.println();
    }

    private static void demonstrateExceptionHandling() {
        System.out.println("Exception Handling");
        System.out.println("==================\n");

        try {
            // Create rational function with pole
            BernsteinPolynomial num = new BernsteinPolynomial(
                new double[]{1.0, 1.0}, 0.0, 1.0);
            BernsteinPolynomial den = new BernsteinPolynomial(
                new double[]{1.0, -1.0}, 0.0, 1.0);
            RationalFunction r = new RationalFunction(num, den);

            System.out.println("Evaluating rational function with pole...");

            // This should work
            double result = r.evaluate(0.0);
            System.out.printf("  R(0.0) = %.6f (OK)\n", result);

            // This should throw (pole near x = 0.5)
            System.out.print("  R(0.5) = ");
            result = r.evaluate(0.5);
            System.out.printf("%.6f\n", result);

            r.close();
            num.close();
            den.close();

        } catch (GelfgrenException e) {
            System.out.println("Exception caught: " + e.getMessage());
            System.out.println("  (This is expected behavior for poles)");
        }
        System.out.println();
    }
}
