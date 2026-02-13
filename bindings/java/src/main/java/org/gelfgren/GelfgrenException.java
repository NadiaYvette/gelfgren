package org.gelfgren;

/**
 * Exception thrown when Gelfgren operations fail.
 */
public class GelfgrenException extends RuntimeException {
    /**
     * Create a new exception with a message.
     *
     * @param message Error message
     */
    public GelfgrenException(String message) {
        super(message);
    }

    /**
     * Create a new exception with a message and cause.
     *
     * @param message Error message
     * @param cause Underlying cause
     */
    public GelfgrenException(String message, Throwable cause) {
        super(message, cause);
    }
}
