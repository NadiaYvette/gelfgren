# Gelfgren Standard ML Bindings

Standard ML bindings for Gelfgren using MLton FFI.

## Requirements

- MLton compiler
- Gelfgren C library

## Building

```bash
# Build Rust library
cd ../..
cargo build --release

# Compile with MLton
cd bindings/sml
mlton -link-opt -L../../target/release -link-opt -lgelfgren examples/basic_usage.mlb

# Run
export LD_LIBRARY_PATH=../../target/release:$LD_LIBRARY_PATH
./basic_usage
```

## Usage

```sml
open Gelfgren

val poly = BernsteinPolynomial.create ([1.0, 2.0, 6.0], 0.0, 1.0)
val y = BernsteinPolynomial.evaluate (poly, 0.5)
val _ = BernsteinPolynomial.free poly
```

## License

MIT or Apache-2.0
