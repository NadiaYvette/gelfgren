# Gelfgren OCaml Bindings

OCaml bindings for the Gelfgren numerical computing library, providing high-performance piecewise rational interpolation methods through type-safe Ctypes FFI.

## Features

- **Type-Safe FFI**: Ctypes for safe and idiomatic C interop
- **Automatic Resource Management**: `with_*` functions ensure cleanup
- **Functional API**: Immutable operations with operator overloading
- **Exception Handling**: OCaml exceptions for error reporting
- **Array Operations**: Vectorized evaluation support
- **Zero-Copy**: Efficient data transfer via Ctypes

## Requirements

- OCaml 4.14.0 or later
- Dune 3.0 or later
- Ctypes 0.20.0 or later
- Gelfgren C library (built from source)

## Installation

### Using opam (when published)

```bash
opam install gelfgren
```

### Building from Source

```bash
# Build the Rust library first
cd ../..
cargo build --release

# Build OCaml bindings
cd bindings/ocaml
dune build

# Run tests
dune test

# Run example
dune exec examples/basic_usage.exe
```

## Usage

### Basic Example

```ocaml
open Gelfgren

let () =
  (* Create a Bernstein polynomial with automatic cleanup *)
  BernsteinPolynomial.with_poly [|1.0; 2.0; 6.0|] 0.0 1.0 (fun poly ->
      (* Evaluate at a point *)
      let y = BernsteinPolynomial.evaluate poly 0.5 in
      Printf.printf "f(0.5) = %.6f\n" y;

      (* Compute derivative *)
      let dpoly = BernsteinPolynomial.derivative poly in
      let dy = BernsteinPolynomial.evaluate dpoly 0.5 in
      Printf.printf "f'(0.5) = %.6f\n" dy;
      BernsteinPolynomial.free dpoly
    )
```

### Polynomial Arithmetic with Operators

```ocaml
open Gelfgren
open BernsteinPolynomial

let p1 = create [|1.0; 2.0; 3.0|] 0.0 1.0 in
let p2 = create [|2.0; 1.0; 1.0|] 0.0 1.0 in

(* Use infix operators *)
let sum = p1 + p2 in
let diff = p1 - p2 in
let prod = p1 * p2 in
let scaled = 2.0 *. p1 in

(* Evaluate *)
let y = evaluate sum 0.5 in

(* Clean up *)
free p1;
free p2;
free sum;
free diff;
free prod;
free scaled
```

### Vectorized Evaluation

```ocaml
(* Evaluate at multiple points *)
let xs = [|0.0; 0.25; 0.5; 0.75; 1.0|] in
let ys = BernsteinPolynomial.evaluate_array poly xs in

Array.iter2
  (fun x y -> Printf.printf "f(%.2f) = %.6f\n" x y)
  xs ys
```

### Rational Functions

```ocaml
(* Create rational function: x / (1 + x) *)
let num = BernsteinPolynomial.create [|0.0; 1.0|] 0.0 1.0 in
let den = BernsteinPolynomial.create [|1.0; 2.0|] 0.0 1.0 in

RationalFunction.with_rational num den (fun rat ->
    let y = RationalFunction.evaluate rat 0.5 in
    Printf.printf "R(0.5) = %.6f\n" y
  );

BernsteinPolynomial.free num;
BernsteinPolynomial.free den
```

### Padé Approximants

```ocaml
(* Approximate exp(x) with [2/2] Padé approximant *)
let coeffs = [|1.0; 1.0; 0.5; 1.0/.6.0; 1.0/.24.0|] in

PadeApproximant.with_pade coeffs 2 2 (-1.0) 1.0 (fun pade ->
    let y = PadeApproximant.evaluate pade 0.5 in
    let exact = exp 0.5 in
    let error = abs_float (y -. exact) in
    Printf.printf "Padé: %.6f, exp: %.6f, error: %.8f\n" y exact error
  )
```

### Error Handling

```ocaml
try
  let poly = BernsteinPolynomial.create [||] 0.0 1.0 in
  BernsteinPolynomial.free poly
with Gelfgren_error msg ->
  Printf.eprintf "Error: %s\n" msg
```

## API Reference

### Module: BernsteinPolynomial

#### Types

```ocaml
type t  (* Bernstein polynomial on interval [a, b] *)
```

#### Functions

```ocaml
val create : float array -> float -> float -> t
(* Create polynomial from coefficients on [a, b] *)

val free : t -> unit
(* Free polynomial resources *)

val with_poly : float array -> float -> float -> (t -> 'a) -> 'a
(* Create, use, and automatically free polynomial *)

val evaluate : t -> float -> float
(* Evaluate at a point *)

val evaluate_array : t -> float array -> float array
(* Evaluate at multiple points *)

val derivative : t -> t
(* Compute derivative (returns new polynomial) *)

val integral : t -> t
(* Compute antiderivative (returns new polynomial) *)

val degree : t -> int
(* Get polynomial degree *)

val add : t -> t -> t
val ( + ) : t -> t -> t
(* Add two polynomials *)

val subtract : t -> t -> t
val ( - ) : t -> t -> t
(* Subtract polynomials *)

val multiply : t -> t -> t
val ( * ) : t -> t -> t
(* Multiply polynomials *)

val scale : float -> t -> t
val ( *. ) : float -> t -> t
(* Scale polynomial by scalar *)
```

### Module: RationalFunction

```ocaml
type t  (* Rational function P(x)/Q(x) *)

val create : BernsteinPolynomial.t -> BernsteinPolynomial.t -> t
val free : t -> unit
val with_rational : BernsteinPolynomial.t -> BernsteinPolynomial.t -> (t -> 'a) -> 'a
val evaluate : t -> float -> float
val evaluate_array : t -> float array -> float array
val derivative : t -> t
```

### Module: PadeApproximant

```ocaml
type t  (* Padé approximant [n/m] *)

val from_series : float array -> int -> int -> float -> float -> t
(* from_series coeffs n m a b - create [n/m] approximant from series *)

val free : t -> unit
val with_pade : float array -> int -> int -> float -> float -> (t -> 'a) -> 'a
val evaluate : t -> float -> float
val evaluate_array : t -> float array -> float array
```

### Exception

```ocaml
exception Gelfgren_error of string
(* Raised for all Gelfgren operation errors *)
```

## Memory Management

The OCaml bindings provide two patterns for memory management:

### Manual Management

```ocaml
let poly = BernsteinPolynomial.create [|1.0; 2.0|] 0.0 1.0 in
(* use poly *)
BernsteinPolynomial.free poly  (* Must call explicitly *)
```

### Automatic Management (Recommended)

```ocaml
BernsteinPolynomial.with_poly [|1.0; 2.0|] 0.0 1.0 (fun poly ->
    (* use poly - automatically freed when function returns *)
  )
```

**Important**: Derived polynomials (from `derivative`, `add`, etc.) are new objects that must be freed separately.

## Performance

The OCaml bindings offer excellent performance:

- **Zero-copy FFI**: Ctypes.CArray for efficient array passing
- **No boxing overhead**: Direct float operations
- **Lazy evaluation**: Ctypes foreign calls are JIT-compiled
- **Native speed**: FFI calls compile to direct C function calls

For best performance:
- Use `evaluate_array` for vectorized operations
- Compile with `dune build --profile release`
- Consider using `flambda` optimizations

## Building Documentation

Generate API documentation:

```bash
dune build @doc
```

View in browser:

```bash
open _build/default/_doc/_html/gelfgren/index.html
```

## Examples

Complete examples are in the `examples/` directory:

- `basic_usage.ml` - Comprehensive demonstration

Run with:

```bash
dune exec examples/basic_usage.exe
```

## Testing

Run the test suite:

```bash
dune test
```

The tests use OUnit2 and cover:
- Polynomial creation and evaluation
- Arithmetic operations
- Calculus operations
- Rational functions
- Padé approximants
- Error handling

## Integration with OCaml Ecosystem

### Using with Core/Base

```ocaml
open Core
open Gelfgren

let evaluate_at_points poly points =
  List.map points ~f:(fun x ->
      (x, BernsteinPolynomial.evaluate poly x))
```

### Using with Owl (Scientific Computing)

```ocaml
open Owl
open Gelfgren

let evaluate_mat poly mat =
  Mat.map (fun x -> BernsteinPolynomial.evaluate poly x) mat
```

## Troubleshooting

### Library Not Found

If you get "cannot find -lgelfgren" errors:

```bash
# Add to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/path/to/gelfgren/target/release:$LD_LIBRARY_PATH

# Or set in dune file
(env
 (dev
  (env-vars
   (LD_LIBRARY_PATH /path/to/gelfgren/target/release))))
```

### Ctypes Version Issues

Ensure you have a compatible Ctypes version:

```bash
opam install ctypes.0.20.0 ctypes-foreign.0.18.0
```

## Contributing

Contributions are welcome! See the main [CONTRIBUTING.md](../../CONTRIBUTING.md).

## License

Licensed under either of MIT or Apache-2.0 at your option.

## See Also

- [Main Gelfgren Documentation](../../README.md)
- [OCaml Ctypes Documentation](https://github.com/ocamllabs/ocaml-ctypes)
- [Dune Documentation](https://dune.readthedocs.io/)
- [C API Header](../../include/gelfgren.h)
