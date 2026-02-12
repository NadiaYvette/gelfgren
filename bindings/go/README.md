# Gelfgren Go Bindings

Go bindings for the Gelfgren piecewise rational interpolation library using cgo.

## Features

- **Idiomatic Go API** with proper error handling
- **Automatic Memory Management** via finalizers
- **Type Safety** with Go's strong typing
- **Vectorized Operations** for arrays
- **Zero-Copy FFI** through cgo
- **Comprehensive Tests** and benchmarks

## Installation

```bash
# Build Rust library first
cd ../..
cargo build --release

# Use Go module
cd bindings/go
export CGO_LDFLAGS="-L../../target/release"
export LD_LIBRARY_PATH="../../target/release:$LD_LIBRARY_PATH"
go build ./gelfgren
go test ./gelfgren
```

## Quick Start

```go
package main

import (
    "fmt"
    "github.com/yourusername/gelfgren-go/gelfgren"
)

func main() {
    // Create polynomial
    poly, err := gelfgren.NewBernsteinPolynomial([]float64{1.0, 2.0, 6.0}, 0.0, 1.0)
    if err != nil {
        panic(err)
    }
    defer poly.Free() // Automatic cleanup

    // Evaluate
    y, _ := poly.Evaluate(0.5)
    fmt.Printf("f(0.5) = %f\n", y)

    // Vectorized
    ys, _ := poly.EvaluateSlice([]float64{0.0, 0.5, 1.0})
}
```

## API Reference

### BernsteinPolynomial

```go
poly, err := gelfgren.NewBernsteinPolynomial(coeffs []float64, a, b float64)
y, err := poly.Evaluate(x float64)
ys, err := poly.EvaluateSlice(xs []float64)
dpoly, err := poly.Derivative()
ipoly, err := poly.Integral()
deg, err := poly.Degree()
sum, err := poly.Add(other *BernsteinPolynomial)
diff, err := poly.Subtract(other *BernsteinPolynomial)
prod, err := poly.Multiply(other *BernsteinPolynomial)
scaled, err := poly.Scale(scalar float64)
poly.Free() // Manual cleanup if needed
```

### RationalFunction

```go
rat, err := gelfgren.NewRationalFunction(num, den *BernsteinPolynomial)
y, err := rat.Evaluate(x float64)
ys, err := rat.EvaluateSlice(xs []float64)
drat, err := rat.Derivative()
rat.Free()
```

### PadeApproximant

```go
pade, err := gelfgren.NewPadeApproximant(coeffs []float64, n, m int, a, b float64)
y, err := pade.Evaluate(x float64)
ys, err := pade.EvaluateSlice(xs []float64)
pade.Free()
```

## Testing

```bash
go test ./gelfgren -v
go test ./gelfgren -bench=.
```

## Examples

```bash
cd examples/basic
go run main.go
```

## License

MIT or Apache-2.0
