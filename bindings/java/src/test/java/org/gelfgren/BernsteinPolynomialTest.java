package org.gelfgren;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.DisplayName;
import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for BernsteinPolynomial.
 */
public class BernsteinPolynomialTest {

    @Test
    @DisplayName("Create and evaluate polynomial")
    public void testCreateAndEvaluate() {
        try (BernsteinPolynomial p = new BernsteinPolynomial(
                new double[]{1.0, 2.0, 3.0}, 0.0, 1.0)) {

            // Test degree
            assertEquals(2, p.getDegree());

            // Test interval
            double[] interval = p.getInterval();
            assertEquals(0.0, interval[0], 1e-10);
            assertEquals(1.0, interval[1], 1e-10);

            // Test evaluation at endpoints
            assertEquals(1.0, p.evaluate(0.0), 1e-10);
            assertEquals(3.0, p.evaluate(1.0), 1e-10);

            // Test evaluation at midpoint
            double mid = p.evaluate(0.5);
            assertTrue(mid > 1.0 && mid < 3.0);
        }
    }

    @Test
    @DisplayName("Polynomial arithmetic")
    public void testArithmetic() {
        try (BernsteinPolynomial p = new BernsteinPolynomial(
                new double[]{1.0, 2.0}, 0.0, 1.0);
             BernsteinPolynomial q = new BernsteinPolynomial(
                new double[]{3.0, 4.0}, 0.0, 1.0)) {

            // Test addition
            try (BernsteinPolynomial sum = p.add(q)) {
                double expected = p.evaluate(0.5) + q.evaluate(0.5);
                assertEquals(expected, sum.evaluate(0.5), 1e-10);
            }

            // Test subtraction
            try (BernsteinPolynomial diff = p.subtract(q)) {
                double expected = p.evaluate(0.5) - q.evaluate(0.5);
                assertEquals(expected, diff.evaluate(0.5), 1e-10);
            }

            // Test multiplication
            try (BernsteinPolynomial prod = p.multiply(q)) {
                double expected = p.evaluate(0.5) * q.evaluate(0.5);
                assertEquals(expected, prod.evaluate(0.5), 1e-10);
            }
        }
    }

    @Test
    @DisplayName("Derivative computation")
    public void testDerivative() {
        try (BernsteinPolynomial p = new BernsteinPolynomial(
                new double[]{1.0, 2.0, 3.0}, 0.0, 1.0);
             BernsteinPolynomial pPrime = p.derivative()) {

            // Derivative should have degree n-1
            assertEquals(1, pPrime.getDegree());

            // Derivative should be non-negative for this polynomial
            assertTrue(pPrime.evaluate(0.0) >= 0.0);
            assertTrue(pPrime.evaluate(0.5) >= 0.0);
            assertTrue(pPrime.evaluate(1.0) >= 0.0);
        }
    }

    @Test
    @DisplayName("Degree elevation")
    public void testElevate() {
        try (BernsteinPolynomial p = new BernsteinPolynomial(
                new double[]{1.0, 2.0}, 0.0, 1.0);
             BernsteinPolynomial elevated = p.elevate()) {

            // Degree should increase by 1
            assertEquals(p.getDegree() + 1, elevated.getDegree());

            // Values should be preserved
            assertEquals(p.evaluate(0.0), elevated.evaluate(0.0), 1e-10);
            assertEquals(p.evaluate(0.5), elevated.evaluate(0.5), 1e-10);
            assertEquals(p.evaluate(1.0), elevated.evaluate(1.0), 1e-10);
        }
    }

    @Test
    @DisplayName("AutoCloseable resource management")
    public void testAutoClose() {
        BernsteinPolynomial p = new BernsteinPolynomial(
            new double[]{1.0, 2.0}, 0.0, 1.0);
        p.close();

        // Should throw after close
        assertThrows(IllegalStateException.class, () -> p.evaluate(0.5));
    }

    @Test
    @DisplayName("Invalid arguments")
    public void testInvalidArguments() {
        // Null coefficients
        assertThrows(IllegalArgumentException.class,
            () -> new BernsteinPolynomial(null, 0.0, 1.0));

        // Empty coefficients
        assertThrows(IllegalArgumentException.class,
            () -> new BernsteinPolynomial(new double[]{}, 0.0, 1.0));
    }
}
