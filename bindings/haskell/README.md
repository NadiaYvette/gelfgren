# Gelfgren Haskell Bindings

Haskell bindings for the Gelfgren numerical computing library, providing high-performance piecewise rational interpolation methods through FFI.

## Features

- **Type-Safe FFI**: Safe Haskell interface with automatic memory management
- **Foreign Pointers**: Automatic cleanup using `ForeignPtr` with finalizers
- **Exception Handling**: Haskell exceptions for error reporting
- **Type Classes**: Generic `Evaluable` type class for uniform API
- **Bracket Pattern**: Resource management with `withBernstein` and friends
- **Pure Haskell**: No C code in user programs

## Requirements

- GHC 9.2 or later
- Cabal 3.0 or later
- Gelfgren C library (built from source)

## Installation

### From Source

```bash
# Build the Rust library first
cd ../..
cargo build --release

# Build Haskell package
cd bindings/haskell
cabal build

# Run tests
cabal test

# Run example
cabal run gelfgren-example
```

### Using Stack (Alternative)

Create a `stack.yaml`:

```yaml
resolver: lts-21.22  # GHC 9.4.8

packages:
- .

extra-deps: []

extra-lib-dirs:
- ../../target/release

extra-include-dirs:
- ../../include
```

Then:

```bash
stack build
stack test
stack run gelfgren-example
```

## Usage

### Basic Example

```haskell
import Gelfgren

main :: IO ()
main = do
    -- Create a Bernstein polynomial
    poly <- bernsteinPolynomial [1.0, 2.0, 6.0] 0.0 1.0

    -- Evaluate at a point
    y <- evaluate poly 0.5
    print y

    -- Compute derivative
    dpoly <- derivative poly
    dy <- evaluate dpoly 0.5
    print dy
```

### Vectorized Evaluation

```haskell
import Gelfgren

main :: IO ()
main = do
    poly <- bernsteinPolynomial [1.0, 2.0, 6.0] 0.0 1.0

    -- Evaluate at multiple points
    let xs = [0.0, 0.25, 0.5, 0.75, 1.0]
    ys <- evaluateMany poly xs

    -- Print results
    mapM_ print (zip xs ys)
```

### Polynomial Arithmetic

```haskell
import Gelfgren

main :: IO ()
main = do
    p1 <- bernsteinPolynomial [1.0, 2.0, 3.0] 0.0 1.0
    p2 <- bernsteinPolynomial [2.0, 1.0, 1.0] 0.0 1.0

    -- Addition
    pSum <- add p1 p2
    ySum <- evaluate pSum 0.5
    print ySum

    -- Subtraction
    pDiff <- subtract p1 p2
    yDiff <- evaluate pDiff 0.5
    print yDiff

    -- Multiplication
    pProd <- multiply p1 p2
    yProd <- evaluate pProd 0.5
    print yProd

    -- Scaling
    pScaled <- scale 2.0 p1
    yScaled <- evaluate pScaled 0.5
    print yScaled
```

### Rational Functions

```haskell
import Gelfgren

main :: IO ()
main = do
    -- Create numerator and denominator
    num <- bernsteinPolynomial [0.0, 1.0] 0.0 1.0
    den <- bernsteinPolynomial [1.0, 2.0] 0.0 1.0

    -- Create rational function: x / (1 + x)
    rat <- rationalFunction num den

    -- Evaluate
    y <- evaluate rat 0.5
    print y  -- 0.333...
```

### Padé Approximants

```haskell
import Gelfgren

main :: IO ()
main = do
    -- Approximate exp(x) with [2/2] Padé approximant
    let coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0]
    pade <- padeFromSeries coeffs 2 2 (-1.0) 1.0

    -- Evaluate
    y <- evaluate pade 0.5
    let exact = exp 0.5
    print $ "Padé: " ++ show y
    print $ "exp:  " ++ show exact
    print $ "Error: " ++ show (abs (y - exact))
```

### Resource Management with Bracket Pattern

```haskell
import Gelfgren

main :: IO ()
main = do
    -- Automatic cleanup with withBernstein
    result <- withBernstein [1.0, 2.0, 6.0] 0.0 1.0 $ \poly -> do
        y <- evaluate poly 0.5
        dpoly <- derivative poly
        dy <- evaluate dpoly 0.5
        return (y, dy)

    let (y, dy) = result
    print $ "f(0.5) = " ++ show y
    print $ "f'(0.5) = " ++ show dy
```

### Error Handling

```haskell
import Gelfgren
import Control.Exception (catch)

main :: IO ()
main = do
    result <- (Just <$> bernsteinPolynomial [1.0, 2.0] 0.0 1.0)
        `catch` \(e :: GelfgrenError) -> do
            print $ "Error: " ++ show e
            return Nothing

    case result of
        Just poly -> do
            y <- evaluate poly 0.5
            print y
        Nothing -> putStrLn "Failed to create polynomial"
```

## API Reference

### Types

#### `BernsteinPolynomial`

Represents a Bernstein polynomial on an interval [a, b].

#### `RationalFunction`

Represents a rational function P(x)/Q(x).

#### `PadeApproximant`

Represents a Padé approximant [n/m].

#### `GelfgrenError`

Exception type for errors:
- `NullPointerError String`
- `InvalidIntervalError String`
- `EmptyCoefficientsError String`
- `SingularMatrixError String`
- `DivisionByZeroError String`
- `PoleError String`
- `UnknownError String`

### Functions

#### Creating Polynomials

```haskell
bernsteinPolynomial :: [Double] -> Double -> Double -> IO BernsteinPolynomial
```

Create a Bernstein polynomial from coefficients on interval [a, b].

```haskell
rationalFunction :: BernsteinPolynomial -> BernsteinPolynomial -> IO RationalFunction
```

Create a rational function from numerator and denominator.

```haskell
padeFromSeries :: [Double] -> Int -> Int -> Double -> Double -> IO PadeApproximant
```

Create a Padé approximant from power series coefficients.

#### Evaluation

```haskell
evaluate :: Evaluable a => a -> Double -> IO Double
```

Evaluate at a single point. Works for `BernsteinPolynomial`, `RationalFunction`, and `PadeApproximant`.

```haskell
evaluateMany :: Evaluable a => a -> [Double] -> IO [Double]
```

Evaluate at multiple points.

#### Calculus

```haskell
derivative :: BernsteinPolynomial -> IO BernsteinPolynomial
```

Compute the derivative.

```haskell
integral :: BernsteinPolynomial -> IO BernsteinPolynomial
```

Compute the antiderivative.

```haskell
degree :: BernsteinPolynomial -> IO Int
```

Get the polynomial degree.

#### Arithmetic

```haskell
add :: BernsteinPolynomial -> BernsteinPolynomial -> IO BernsteinPolynomial
```

Add two polynomials.

```haskell
subtract :: BernsteinPolynomial -> BernsteinPolynomial -> IO BernsteinPolynomial
```

Subtract two polynomials.

```haskell
multiply :: BernsteinPolynomial -> BernsteinPolynomial -> IO BernsteinPolynomial
```

Multiply two polynomials.

```haskell
scale :: Double -> BernsteinPolynomial -> IO BernsteinPolynomial
```

Multiply a polynomial by a scalar.

#### Resource Management

```haskell
withBernstein :: [Double] -> Double -> Double -> (BernsteinPolynomial -> IO a) -> IO a
```

Bracket-style helper for automatic cleanup.

```haskell
withRational :: BernsteinPolynomial -> BernsteinPolynomial -> (RationalFunction -> IO a) -> IO a
```

Bracket-style helper for rational functions.

```haskell
withPade :: [Double] -> Int -> Int -> Double -> Double -> (PadeApproximant -> IO a) -> IO a
```

Bracket-style helper for Padé approximants.

## Memory Management

The library uses `ForeignPtr` with automatic finalizers, so memory is automatically cleaned up by the garbage collector. However, for long-running programs, you may want to use the bracket-style helpers (`withBernstein`, etc.) to ensure timely cleanup.

## Type Classes

### `Evaluable`

Generic type class for objects that can be evaluated:

```haskell
class Evaluable a where
    evaluateImpl :: a -> Double -> IO Double
```

Instances: `BernsteinPolynomial`, `RationalFunction`, `PadeApproximant`

## Performance

The Haskell bindings have minimal overhead:
- FFI calls are fast (marked `unsafe` where appropriate)
- `ForeignPtr` adds negligible overhead
- No data copying for evaluation

For best performance:
- Use `evaluateMany` for multiple points
- Reuse polynomial objects
- Compile with optimizations (`-O2`)

## GHC Compatibility

Tested with:
- GHC 9.2.8
- GHC 9.4.8
- GHC 9.6.3

## Building Documentation

Generate Haddock documentation:

```bash
cabal haddock
```

View documentation:

```bash
open dist-newstyle/build/*/ghc-*/gelfgren-*/doc/html/gelfgren/index.html
```

## Testing

Run the test suite:

```bash
cabal test --test-show-details=direct
```

With coverage:

```bash
cabal test --enable-coverage
```

## Examples

Complete working examples are in the `app/` directory:
- `Main.hs` - Comprehensive demonstration of all features

Run with:

```bash
cabal run gelfgren-example
```

## Troubleshooting

### Library Not Found

If you get "cannot load library" errors:

```bash
export LD_LIBRARY_PATH=/path/to/gelfgren/target/release:$LD_LIBRARY_PATH
# On macOS:
export DYLD_LIBRARY_PATH=/path/to/gelfgren/target/release:$DYLD_LIBRARY_PATH
```

Or configure in `cabal.project.local`:

```cabal
package gelfgren
  extra-lib-dirs: /absolute/path/to/gelfgren/target/release
  extra-include-dirs: /absolute/path/to/gelfgren/include
```

### Linking Issues

If Cabal can't find the library, create `cabal.project.local`:

```cabal
package gelfgren
  extra-lib-dirs: ../../target/release
  extra-include-dirs: ../../include
```

### GHC Version Issues

The package requires GHC 9.2 or later. Check your version:

```bash
ghc --version
```

Use Stack or ghcup to install the correct version.

## Contributing

Contributions are welcome! Please see the main [CONTRIBUTING.md](../../CONTRIBUTING.md) for guidelines.

## License

Licensed under either of MIT or Apache-2.0 at your option.

## See Also

- [Main Gelfgren Documentation](https://yourusername.github.io/gelfgren)
- [Rust API Documentation](https://docs.rs/gelfgren-core)
- [C API Header](../../include/gelfgren.h)
- [Mathematical Background](../../docs/mathematics.md)
- [Haddock Documentation](dist-newstyle/build/*/ghc-*/gelfgren-*/doc/html/gelfgren/index.html)
