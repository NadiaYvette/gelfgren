package org.gelfgren;

/**
 * Bernstein polynomial representation.
 *
 * Provides numerically stable polynomial operations using the Bernstein basis.
 * Polynomials are defined on an interval [a, b] with unscaled coefficients.
 *
 * <p>Example usage:
 * <pre>{@code
 * // Create polynomial P(x) = 1 + 2x + 3x²
 * double[] coeffs = {1.0, 2.0, 3.0};
 * BernsteinPolynomial p = new BernsteinPolynomial(coeffs, 0.0, 1.0);
 *
 * // Evaluate at a point
 * double value = p.evaluate(0.5);
 *
 * // Get polynomial properties
 * int degree = p.getDegree();
 * double[] interval = p.getInterval();
 *
 * // Compute derivative
 * BernsteinPolynomial pPrime = p.derivative();
 *
 * // Polynomial arithmetic
 * BernsteinPolynomial q = new BernsteinPolynomial(new double[]{1.0, 1.0}, 0.0, 1.0);
 * BernsteinPolynomial sum = p.add(q);
 *
 * // Clean up (or let GC handle it)
 * p.close();
 * }</pre>
 */
public class BernsteinPolynomial implements AutoCloseable {
    private long nativeHandle;
    private boolean closed = false;

    /**
     * Create a Bernstein polynomial from unscaled coefficients.
     *
     * @param coeffs Unscaled Bernstein coefficients [C₀, C₁, ..., Cₙ]
     * @param a Left endpoint of interval
     * @param b Right endpoint of interval
     * @throws GelfgrenException if creation fails
     */
    public BernsteinPolynomial(double[] coeffs, double a, double b) {
        if (coeffs == null || coeffs.length == 0) {
            throw new IllegalArgumentException("Coefficients cannot be null or empty");
        }
        this.nativeHandle = create(coeffs, coeffs.length - 1, a, b);
        if (this.nativeHandle == 0) {
            throw new GelfgrenException("Failed to create BernsteinPolynomial: "
                + Gelfgren.getLastErrorMessage());
        }
    }

    /**
     * Internal constructor for wrapping existing native handle.
     */
    private BernsteinPolynomial(long nativeHandle) {
        if (nativeHandle == 0) {
            throw new GelfgrenException("Invalid native handle");
        }
        this.nativeHandle = nativeHandle;
    }

    /**
     * Evaluate the polynomial at a point.
     *
     * @param x Evaluation point
     * @return P(x)
     * @throws GelfgrenException if evaluation fails
     */
    public double evaluate(double x) {
        checkClosed();
        double[] result = new double[1];
        int code = evaluateNative(nativeHandle, x, result);
        if (code != 0) {
            throw new GelfgrenException("Evaluation failed: "
                + Gelfgren.getLastErrorMessage());
        }
        return result[0];
    }

    /**
     * Evaluate the polynomial at multiple points.
     *
     * @param x Array of evaluation points
     * @return Array of P(x) values
     */
    public double[] evaluate(double[] x) {
        checkClosed();
        double[] result = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            result[i] = evaluate(x[i]);
        }
        return result;
    }

    /**
     * Get the degree of the polynomial.
     *
     * @return Polynomial degree
     */
    public int getDegree() {
        checkClosed();
        long[] result = new long[1];
        int code = degreeNative(nativeHandle, result);
        if (code != 0) {
            throw new GelfgrenException("Failed to get degree: "
                + Gelfgren.getLastErrorMessage());
        }
        return (int) result[0];
    }

    /**
     * Get the interval [a, b].
     *
     * @return Array [a, b]
     */
    public double[] getInterval() {
        checkClosed();
        double[] result = new double[2];
        int code = intervalNative(nativeHandle, result);
        if (code != 0) {
            throw new GelfgrenException("Failed to get interval: "
                + Gelfgren.getLastErrorMessage());
        }
        return result;
    }

    /**
     * Compute the derivative polynomial.
     *
     * @return Derivative polynomial P'(x)
     */
    public BernsteinPolynomial derivative() {
        checkClosed();
        long handle = derivativeNative(nativeHandle);
        if (handle == 0) {
            throw new GelfgrenException("Failed to compute derivative: "
                + Gelfgren.getLastErrorMessage());
        }
        return new BernsteinPolynomial(handle);
    }

    /**
     * Compute the integral (antiderivative) polynomial.
     *
     * @return Integral polynomial
     */
    public BernsteinPolynomial integral() {
        checkClosed();
        long handle = integralNative(nativeHandle);
        if (handle == 0) {
            throw new GelfgrenException("Failed to compute integral: "
                + Gelfgren.getLastErrorMessage());
        }
        return new BernsteinPolynomial(handle);
    }

    /**
     * Elevate the polynomial degree by 1.
     *
     * @return Elevated polynomial
     */
    public BernsteinPolynomial elevate() {
        checkClosed();
        long handle = elevateNative(nativeHandle);
        if (handle == 0) {
            throw new GelfgrenException("Failed to elevate degree: "
                + Gelfgren.getLastErrorMessage());
        }
        return new BernsteinPolynomial(handle);
    }

    /**
     * Add two polynomials.
     *
     * @param other Polynomial to add
     * @return Sum polynomial P + Q
     */
    public BernsteinPolynomial add(BernsteinPolynomial other) {
        checkClosed();
        other.checkClosed();
        long handle = addNative(nativeHandle, other.nativeHandle);
        if (handle == 0) {
            throw new GelfgrenException("Failed to add polynomials: "
                + Gelfgren.getLastErrorMessage());
        }
        return new BernsteinPolynomial(handle);
    }

    /**
     * Subtract two polynomials.
     *
     * @param other Polynomial to subtract
     * @return Difference polynomial P - Q
     */
    public BernsteinPolynomial subtract(BernsteinPolynomial other) {
        checkClosed();
        other.checkClosed();
        long handle = subtractNative(nativeHandle, other.nativeHandle);
        if (handle == 0) {
            throw new GelfgrenException("Failed to subtract polynomials: "
                + Gelfgren.getLastErrorMessage());
        }
        return new BernsteinPolynomial(handle);
    }

    /**
     * Multiply two polynomials.
     *
     * @param other Polynomial to multiply
     * @return Product polynomial P * Q
     */
    public BernsteinPolynomial multiply(BernsteinPolynomial other) {
        checkClosed();
        other.checkClosed();
        long handle = multiplyNative(nativeHandle, other.nativeHandle);
        if (handle == 0) {
            throw new GelfgrenException("Failed to multiply polynomials: "
                + Gelfgren.getLastErrorMessage());
        }
        return new BernsteinPolynomial(handle);
    }

    /**
     * Check if this polynomial has been closed.
     */
    private void checkClosed() {
        if (closed) {
            throw new IllegalStateException("BernsteinPolynomial has been closed");
        }
    }

    /**
     * Close and free the native resources.
     */
    @Override
    public void close() {
        if (!closed && nativeHandle != 0) {
            free(nativeHandle);
            closed = true;
            nativeHandle = 0;
        }
    }

    @Override
    protected void finalize() throws Throwable {
        try {
            close();
        } finally {
            super.finalize();
        }
    }

    @Override
    public String toString() {
        if (closed) {
            return "BernsteinPolynomial[closed]";
        }
        try {
            double[] interval = getInterval();
            return String.format("BernsteinPolynomial[degree=%d, interval=[%.2f, %.2f]]",
                getDegree(), interval[0], interval[1]);
        } catch (Exception e) {
            return "BernsteinPolynomial[error]";
        }
    }

    // Native method declarations
    private static native long create(double[] coeffs, int degree, double a, double b);
    private static native void free(long handle);
    private static native int evaluateNative(long handle, double x, double[] result);
    private static native int degreeNative(long handle, long[] result);
    private static native int intervalNative(long handle, double[] result);
    private static native long derivativeNative(long handle);
    private static native long integralNative(long handle);
    private static native long elevateNative(long handle);
    private static native long addNative(long handle1, long handle2);
    private static native long subtractNative(long handle1, long handle2);
    private static native long multiplyNative(long handle1, long handle2);
}
