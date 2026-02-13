# Gelfgren Mercury Bindings

Mercury bindings for the Gelfgren numerical computing library, providing high-performance piecewise rational interpolation methods through foreign_proc pragmas.

## Features

- **Type-Safe FFI**: Mercury's foreign_proc for safe C interop
- **Opaque Types**: Proper abstraction of C pointers
- **Maybe Error Types**: Idiomatic Mercury error handling
- **Pure Predicates**: All operations marked pure where applicable
- **Memory Management**: Explicit cleanup with free predicates
- **List-Based API**: Natural Mercury list handling

## Requirements

- Mercury compiler 22.01 or later
- Gelfgren C library (built from source)
- C compiler (GCC, Clang, or MSVC)

## Installation

### Building

```bash
# Build the Rust library first
cd ../..
cargo build --release

# Build Mercury module
cd bindings/mercury
mmc --make libgelfgren

# Build examples
mmc --make basic_usage \
    --mld ../../target/release \
    -L ../../target/release \
    -l gelfgren

# Run example (set library path)
export LD_LIBRARY_PATH=../../target/release:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=../../target/release:$DYLD_LIBRARY_PATH
./basic_usage
```

### Using Mmakefile

```bash
# Build library
mmake library

# Build examples
mmake examples

# Run example
mmake run-example

# Clean
mmake clean
```

## Usage

### Basic Example

```mercury
:- module example.
:- interface.
:- import_module io.

:- pred main(io::di, io::uo) is det.

:- implementation.
:- import_module gelfgren.
:- import_module maybe.

main(!IO) :-
    % Create a Bernstein polynomial
    Coeffs = [1.0, 2.0, 6.0],
    create_bernstein(Coeffs, 0.0, 1.0, Result, !IO),

    (
        Result = ok(Poly),

        % Evaluate at a point
        evaluate_bernstein(Poly, 0.5, EvalResult, !IO),
        (
            EvalResult = ok(Y),
            io.format("f(0.5) = %.6f\n", [f(Y)], !IO)
        ;
            EvalResult = error(Err),
            io.write_string(error_to_string(Err), !IO)
        ),

        % Clean up
        free_bernstein(Poly, !IO)
    ;
        Result = error(Err),
        io.write_string(error_to_string(Err), !IO)
    ).
```

### Error Handling

All operations return `gelfgren_result(T)` which is defined as:

```mercury
:- type gelfgren_result(T) == maybe_error(T, gelfgren_error).
```

Use pattern matching to handle results:

```mercury
create_bernstein(Coeffs, A, B, Result, !IO),
(
    Result = ok(Poly),
    % Use Poly
    free_bernstein(Poly, !IO)
;
    Result = error(Err),
    io.write_string(error_to_string(Err), !IO)
).
```

### Vectorized Evaluation

```mercury
% Evaluate at multiple points
XValues = [0.0, 0.25, 0.5, 0.75, 1.0],
list.foldl(evaluate_and_print(Poly), XValues, !IO).

:- pred evaluate_and_print(bernstein_polynomial::in, float::in,
    io::di, io::uo) is det.

evaluate_and_print(Poly, X, !IO) :-
    evaluate_bernstein(Poly, X, Result, !IO),
    (
        Result = ok(Y),
        io.format("f(%.2f) = %.6f\n", [f(X), f(Y)], !IO)
    ;
        Result = error(Err),
        io.write_string(error_to_string(Err), !IO)
    ).
```

### Polynomial Arithmetic

```mercury
% Create two polynomials
create_bernstein([1.0, 2.0, 3.0], 0.0, 1.0, P1Result, !IO),
create_bernstein([2.0, 1.0, 1.0], 0.0, 1.0, P2Result, !IO),

(
    P1Result = ok(P1),
    P2Result = ok(P2)
->
    % Addition
    add_bernstein(P1, P2, SumResult, !IO),
    (
        SumResult = ok(Sum),
        evaluate_bernstein(Sum, 0.5, EvalResult, !IO),
        % ... use result ...
        free_bernstein(Sum, !IO)
    ;
        SumResult = error(Err),
        io.write_string(error_to_string(Err), !IO)
    ),

    % Scaling
    scale_bernstein(2.0, P1, ScaledResult, !IO),
    (
        ScaledResult = ok(Scaled),
        % ... use scaled polynomial ...
        free_bernstein(Scaled, !IO)
    ;
        ScaledResult = error(Err),
        io.write_string(error_to_string(Err), !IO)
    ),

    free_bernstein(P1, !IO),
    free_bernstein(P2, !IO)
;
    io.write_string("Error creating polynomials\n", !IO)
).
```

### Rational Functions

```mercury
% Create numerator and denominator
create_bernstein([0.0, 1.0], 0.0, 1.0, NumResult, !IO),
create_bernstein([1.0, 2.0], 0.0, 1.0, DenResult, !IO),

(
    NumResult = ok(Num),
    DenResult = ok(Den)
->
    % Create rational function: x / (1 + x)
    create_rational(Num, Den, RatResult, !IO),
    (
        RatResult = ok(Rat),

        % Evaluate
        evaluate_rational(Rat, 0.5, EvalResult, !IO),
        (
            EvalResult = ok(Y),
            io.format("R(0.5) = %.6f\n", [f(Y)], !IO)
        ;
            EvalResult = error(Err),
            io.write_string(error_to_string(Err), !IO)
        ),

        free_rational(Rat, !IO)
    ;
        RatResult = error(Err),
        io.write_string(error_to_string(Err), !IO)
    ),

    free_bernstein(Num, !IO),
    free_bernstein(Den, !IO)
;
    io.write_string("Error creating polynomials\n", !IO)
).
```

### Padé Approximants

```mercury
% Approximate exp(x) with [2/2] Padé approximant
Coeffs = [1.0, 1.0, 0.5, 1.0/6.0, 1.0/24.0],
create_pade_from_series(Coeffs, 2, 2, -1.0, 1.0, PadeResult, !IO),

(
    PadeResult = ok(Pade),

    % Evaluate
    evaluate_pade(Pade, 0.5, EvalResult, !IO),
    (
        EvalResult = ok(Y),
        Exact = math.exp(0.5),
        Error = float.abs(Y - Exact),
        io.format("Padé: %.6f, exp: %.6f, error: %.8f\n",
            [f(Y), f(Exact), f(Error)], !IO)
    ;
        EvalResult = error(Err),
        io.write_string(error_to_string(Err), !IO)
    ),

    free_pade(Pade, !IO)
;
    PadeResult = error(Err),
    io.write_string(error_to_string(Err), !IO)
).
```

## API Reference

### Types

#### `bernstein_polynomial`

Opaque type representing a Bernstein polynomial on interval [a, b].

#### `rational_function`

Opaque type representing a rational function P(x)/Q(x).

#### `pade_approximant`

Opaque type representing a Padé approximant [n/m].

#### `gelfgren_error`

Discriminated union for errors:
- `null_pointer_error(string)`
- `invalid_interval_error(string)`
- `empty_coefficients_error(string)`
- `singular_matrix_error(string)`
- `division_by_zero_error(string)`
- `pole_error(string)`
- `unknown_error(string)`

#### `gelfgren_result(T)`

Type alias for `maybe_error(T, gelfgren_error)`.

### Predicates

#### Bernstein Polynomial

```mercury
:- pred create_bernstein(list(float)::in, float::in, float::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.
```

Create a Bernstein polynomial from coefficients on interval [a, b].

```mercury
:- pred free_bernstein(bernstein_polynomial::in, io::di, io::uo) is det.
```

Free a Bernstein polynomial (required to avoid memory leaks).

```mercury
:- pred evaluate_bernstein(bernstein_polynomial::in, float::in,
    gelfgren_result(float)::out, io::di, io::uo) is det.
```

Evaluate the polynomial at a point.

```mercury
:- pred derivative_bernstein(bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.
```

Compute the derivative.

```mercury
:- pred integral_bernstein(bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.
```

Compute the antiderivative.

```mercury
:- pred degree_bernstein(bernstein_polynomial::in,
    gelfgren_result(int)::out, io::di, io::uo) is det.
```

Get the polynomial degree.

```mercury
:- pred add_bernstein(bernstein_polynomial::in, bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.
```

Add two polynomials.

```mercury
:- pred subtract_bernstein(bernstein_polynomial::in, bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.
```

Subtract two polynomials.

```mercury
:- pred multiply_bernstein(bernstein_polynomial::in, bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.
```

Multiply two polynomials.

```mercury
:- pred scale_bernstein(float::in, bernstein_polynomial::in,
    gelfgren_result(bernstein_polynomial)::out, io::di, io::uo) is det.
```

Multiply a polynomial by a scalar.

#### Rational Function

```mercury
:- pred create_rational(bernstein_polynomial::in, bernstein_polynomial::in,
    gelfgren_result(rational_function)::out, io::di, io::uo) is det.
```

Create a rational function from numerator and denominator.

```mercury
:- pred free_rational(rational_function::in, io::di, io::uo) is det.
```

Free a rational function.

```mercury
:- pred evaluate_rational(rational_function::in, float::in,
    gelfgren_result(float)::out, io::di, io::uo) is det.
```

Evaluate the rational function at a point.

```mercury
:- pred derivative_rational(rational_function::in,
    gelfgren_result(rational_function)::out, io::di, io::uo) is det.
```

Compute the derivative.

#### Padé Approximant

```mercury
:- pred create_pade_from_series(list(float)::in, int::in, int::in,
    float::in, float::in, gelfgren_result(pade_approximant)::out,
    io::di, io::uo) is det.
```

Create a Padé approximant from power series coefficients.

```mercury
:- pred free_pade(pade_approximant::in, io::di, io::uo) is det.
```

Free a Padé approximant.

```mercury
:- pred evaluate_pade(pade_approximant::in, float::in,
    gelfgren_result(float)::out, io::di, io::uo) is det.
```

Evaluate the Padé approximant at a point.

#### Utilities

```mercury
:- func error_to_string(gelfgren_error) = string.
```

Convert error to string for display.

## Memory Management

**Important**: Unlike some Mercury programs, these bindings require explicit cleanup. Always call `free_*` predicates when done with objects to avoid memory leaks:

```mercury
create_bernstein(Coeffs, A, B, Result, !IO),
(
    Result = ok(Poly),
    % ... use Poly ...
    free_bernstein(Poly, !IO)  % Required!
;
    Result = error(_),
    % Nothing to free
    true
).
```

## Performance

The Mercury bindings have minimal overhead:
- Direct C calls via foreign_proc
- No data copying for opaque types
- Efficient list-to-array conversion
- Thread-safe where marked

For best performance:
- Reuse polynomial objects when possible
- Compile with optimizations (`-O4` or higher)
- Use appropriate grade for your needs

## Mercury Grades

The bindings work with most Mercury grades. Recommended grades:
- `asm_fast.gc` - Fast with garbage collection (default)
- `hlc.gc` - High-level C backend with GC
- `hl` - High-level without GC (manual memory management)

## Building Documentation

Generate Mercury documentation:

```bash
mmc --make libgelfgren.mh
mercury_config --print-doc-generator
# Follow Mercury documentation generation instructions
```

## Examples

Complete working examples are in the `examples/` directory:
- `basic_usage.m` - Comprehensive demonstration of all features

Run with:

```bash
mmake run-example
```

## Troubleshooting

### Library Not Found

If you get "cannot open shared object file" errors:

```bash
export LD_LIBRARY_PATH=/path/to/gelfgren/target/release:$LD_LIBRARY_PATH
# On macOS:
export DYLD_LIBRARY_PATH=/path/to/gelfgren/target/release:$DYLD_LIBRARY_PATH
```

Or add to Mmakefile:
```make
RPATH_OPT = -Wl,-rpath,/absolute/path/to/gelfgren/target/release
MLLIBS += $(RPATH_OPT)
```

### Compilation Errors

Ensure Mercury compiler can find the C library:

```bash
mmc --make your_module \
    -I ../../include \
    -L ../../target/release \
    -l gelfgren
```

### Mercury Version Issues

Check your Mercury version:

```bash
mmc --version
```

The bindings require Mercury 22.01 or later.

## Contributing

Contributions are welcome! Please see the main [CONTRIBUTING.md](../../CONTRIBUTING.md) for guidelines.

## License

Licensed under either of MIT or Apache-2.0 at your option.

## See Also

- [Main Gelfgren Documentation](https://yourusername.github.io/gelfgren)
- [Mercury Language](https://mercurylang.org/)
- [Rust API Documentation](https://docs.rs/gelfgren-core)
- [C API Header](../../include/gelfgren.h)
- [Mathematical Background](../../docs/mathematics.md)
