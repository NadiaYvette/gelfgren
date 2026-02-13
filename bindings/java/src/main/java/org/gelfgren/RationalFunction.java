package org.gelfgren;

/**
 * Rational function representation P(x)/Q(x).
 *
 * Represents a rational function as the ratio of two Bernstein polynomials.
 * Automatically detects and reports poles (denominator zeros).
 *
 * <p>Example usage:
 * <pre>{@code
 * BernsteinPolynomial num = new BernsteinPolynomial(new double[]{1.0, 1.0}, 0.0, 1.0);
 * BernsteinPolynomial den = new BernsteinPolynomial(new double[]{1.0, 2.0}, 0.0, 1.0);
 *
 * try (RationalFunction r = new RationalFunction(num, den)) {
 *     double value = r.evaluate(0.5);
 *     RationalFunction rPrime = r.derivative();
 * }
 * }</pre>
 */
public class RationalFunction implements AutoCloseable {
    private long nativeHandle;
    private boolean closed = false;

    /**
     * Create a rational function from numerator and denominator polynomials.
     *
     * @param numerator Numerator polynomial P(x)
     * @param denominator Denominator polynomial Q(x)
     * @throws GelfgrenException if creation fails
     */
    public RationalFunction(BernsteinPolynomial numerator, BernsteinPolynomial denominator) {
        if (numerator == null || denominator == null) {
            throw new IllegalArgumentException("Numerator and denominator cannot be null");
        }
        // Note: We don't close numerator/denominator as they're passed in by reference
        this.nativeHandle = create(numerator.nativeHandle, denominator.nativeHandle);
        if (this.nativeHandle == 0) {
            throw new GelfgrenException("Failed to create RationalFunction: "
                + Gelfgren.getLastErrorMessage());
        }
    }

    /**
     * Internal constructor for wrapping existing native handle.
     */
    private RationalFunction(long nativeHandle) {
        if (nativeHandle == 0) {
            throw new GelfgrenException("Invalid native handle");
        }
        this.nativeHandle = nativeHandle;
    }

    /**
     * Evaluate the rational function at a point.
     *
     * @param x Evaluation point
     * @return R(x) = P(x)/Q(x)
     * @throws GelfgrenException if denominator is zero at x
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
     * Evaluate the rational function at multiple points.
     *
     * @param x Array of evaluation points
     * @return Array of R(x) values
     * @throws GelfgrenException if denominator is zero at any point
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
     * Compute the derivative using the quotient rule.
     *
     * @return Derivative R'(x)
     */
    public RationalFunction derivative() {
        checkClosed();
        long handle = derivativeNative(nativeHandle);
        if (handle == 0) {
            throw new GelfgrenException("Failed to compute derivative: "
                + Gelfgren.getLastErrorMessage());
        }
        return new RationalFunction(handle);
    }

    /**
     * Check if this rational function has been closed.
     */
    private void checkClosed() {
        if (closed) {
            throw new IllegalStateException("RationalFunction has been closed");
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
        return closed ? "RationalFunction[closed]" : "RationalFunction[P/Q]";
    }

    // Native method declarations
    private static native long create(long numeratorHandle, long denominatorHandle);
    private static native void free(long handle);
    private static native int evaluateNative(long handle, double x, double[] result);
    private static native long derivativeNative(long handle);
}
