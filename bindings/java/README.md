# Gelfgren Java Bindings

Java bindings for the Gelfgren piecewise rational interpolation library using JNI.

## Requirements

- Java 11 or higher
- Maven 3.6+
- Rust (for building the native library)

## Building

```bash
# Build native library first
cd ../..
cargo build --release

# Build Java library and run tests
cd bindings/java
mvn clean package

# Run example
mvn exec:java -Dexec.mainClass="org.gelfgren.Example"
```

## Quick Start

```java
import org.gelfgren.*;

// Create a Bernstein polynomial
try (BernsteinPolynomial p = new BernsteinPolynomial(
        new double[]{1.0, 2.0, 3.0}, 0.0, 1.0)) {

    // Evaluate at a point
    double value = p.evaluate(0.5);

    // Compute derivative
    BernsteinPolynomial pPrime = p.derivative();

    // Polynomial arithmetic
    BernsteinPolynomial q = new BernsteinPolynomial(
        new double[]{1.0, 1.0}, 0.0, 1.0);
    BernsteinPolynomial sum = p.add(q);
}
```

## Features

- **BernsteinPolynomial**: Numerically stable polynomial operations
- **RationalFunction**: P(x)/Q(x) with automatic pole detection
- **AutoCloseable**: Proper resource management with try-with-resources
- **Exception Handling**: Java-native exception types

## API Documentation

Generate Javadoc:

```bash
mvn javadoc:javadoc
# Open target/site/apidocs/index.html
```

## Testing

```bash
mvn test
```

## License

MIT OR Apache-2.0
