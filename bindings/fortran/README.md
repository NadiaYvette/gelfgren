# Gelfgren Fortran Bindings

Fortran bindings for the Gelfgren numerical computing library, providing high-performance piecewise rational interpolation methods through ISO_C_BINDING.

## Features

- **ISO_C_BINDING Interface**: Standard Fortran 2003 C interoperability
- **Object-Oriented API**: Type-bound procedures for intuitive usage
- **Error Handling**: Optional error reporting with detailed messages
- **Array Support**: Vectorized evaluation of polynomials and functions
- **Memory Safety**: Automatic cleanup with type-bound destructors

## Requirements

- Fortran compiler with ISO_C_BINDING support (gfortran 8+, ifort 19+)
- Gelfgren C library (built from source)
- Fortran Package Manager (fpm) - optional, recommended

## Installation

### Using Fortran Package Manager (fpm)

```bash
# Install fpm if not already installed
cargo install fpm

# Build the library
cd bindings/fortran
fpm build --profile release

# Run tests
fpm test

# Run examples
fpm run basic_usage
```

### Using Make

```bash
# Build the Rust library first
cd ../..
cargo build --release

# Build Fortran module and example
cd bindings/fortran
make

# Run example
make run-example

# Run tests
make run-test
```

### Manual Compilation

```bash
# Build Rust library
cd ../..
cargo build --release

# Compile Fortran module
cd bindings/fortran
gfortran -c src/gelfgren_mod.f90 -o build/gelfgren_mod.o -Jbuild

# Compile and link example
gfortran -Ibuild example/basic_usage.f90 build/gelfgren_mod.o \
    -L../../target/release -lgelfgren -o build/basic_usage

# Run (may need to set LD_LIBRARY_PATH)
export LD_LIBRARY_PATH=../../target/release:$LD_LIBRARY_PATH
./build/basic_usage
```

## Usage

### Basic Example

```fortran
program example
    use gelfgren_mod
    implicit none

    type(bernstein_polynomial) :: poly
    type(gelfgren_error) :: error
    real(8) :: coeffs(3), y

    ! Create polynomial
    coeffs = [1.0d0, 2.0d0, 6.0d0]
    call poly%create(coeffs, 0.0d0, 1.0d0, error)

    ! Check for errors
    if (error%code /= GELFGREN_SUCCESS) then
        print *, "Error:", trim(error%message)
        stop 1
    end if

    ! Evaluate
    y = poly%evaluate(0.5d0, error)
    print *, "f(0.5) =", y

    ! Clean up
    call poly%free()
end program example
```

### Vectorized Evaluation

```fortran
program vectorized
    use gelfgren_mod
    implicit none

    type(bernstein_polynomial) :: poly
    type(gelfgren_error) :: error
    real(8) :: coeffs(3)
    real(8), allocatable :: x(:), y(:)
    integer :: i

    coeffs = [1.0d0, 2.0d0, 6.0d0]
    call poly%create(coeffs, 0.0d0, 1.0d0, error)

    ! Evaluate at multiple points
    allocate(x(100))
    x = [(i/99.0d0, i = 0, 99)]
    y = poly%evaluate(x, error)

    ! Use results
    do i = 1, 10
        print *, "f(", x(i), ") =", y(i)
    end do

    call poly%free()
    deallocate(x, y)
end program vectorized
```

### Derivative and Integral

```fortran
program calculus
    use gelfgren_mod
    implicit none

    type(bernstein_polynomial) :: poly, dpoly, ipoly
    type(gelfgren_error) :: error
    real(8) :: coeffs(3), dy, area

    coeffs = [1.0d0, 2.0d0, 6.0d0]
    call poly%create(coeffs, 0.0d0, 1.0d0, error)

    ! Derivative
    dpoly = poly%derivative()
    dy = dpoly%evaluate(0.5d0, error)
    print *, "f'(0.5) =", dy

    ! Integral
    ipoly = poly%integral()
    area = ipoly%evaluate(1.0d0, error) - ipoly%evaluate(0.0d0, error)
    print *, "Area =", area

    ! Clean up
    call poly%free()
    call dpoly%free()
    call ipoly%free()
end program calculus
```

### Polynomial Arithmetic

```fortran
program arithmetic
    use gelfgren_mod
    implicit none

    type(bernstein_polynomial) :: p1, p2, result_poly
    type(gelfgren_error) :: error
    real(8) :: coeffs(3), y

    ! Create polynomials
    coeffs = [1.0d0, 2.0d0, 3.0d0]
    call p1%create(coeffs, 0.0d0, 1.0d0, error)

    coeffs = [2.0d0, 1.0d0, 1.0d0]
    call p2%create(coeffs, 0.0d0, 1.0d0, error)

    ! Addition
    result_poly = p1%add(p2)
    y = result_poly%evaluate(0.5d0, error)
    print *, "p1 + p2 at x=0.5:", y
    call result_poly%free()

    ! Multiplication
    result_poly = p1%multiply(p2)
    y = result_poly%evaluate(0.5d0, error)
    print *, "p1 * p2 at x=0.5:", y
    call result_poly%free()

    ! Scaling
    result_poly = p1%scale(2.0d0)
    y = result_poly%evaluate(0.5d0, error)
    print *, "2 * p1 at x=0.5:", y
    call result_poly%free()

    call p1%free()
    call p2%free()
end program arithmetic
```

### Rational Functions

```fortran
program rational_example
    use gelfgren_mod
    implicit none

    type(bernstein_polynomial) :: num, den
    type(rational_function) :: rat
    type(gelfgren_error) :: error
    real(8) :: coeffs(2), y

    ! Create numerator: x
    coeffs = [0.0d0, 1.0d0]
    call num%create(coeffs, 0.0d0, 1.0d0, error)

    ! Create denominator: 1 + x
    coeffs = [1.0d0, 2.0d0]
    call den%create(coeffs, 0.0d0, 1.0d0, error)

    ! Create rational: x / (1 + x)
    call rat%create(num, den, error)

    ! Evaluate
    y = rat%evaluate(0.5d0, error)
    print *, "R(0.5) =", y

    ! Clean up
    call num%free()
    call den%free()
    call rat%free()
end program rational_example
```

### Padé Approximants

```fortran
program pade_example
    use gelfgren_mod
    implicit none

    type(pade_approximant) :: pade
    type(gelfgren_error) :: error
    real(8) :: coeffs(5), x, y

    ! Approximate exp(x) with [2/2] Padé
    coeffs = [1.0d0, 1.0d0, 0.5d0, 1.0d0/6.0d0, 1.0d0/24.0d0]
    call pade%create_from_series(coeffs, 2, 2, -1.0d0, 1.0d0, error)

    ! Evaluate
    x = 0.5d0
    y = pade%evaluate(x, error)
    print *, "Padé(", x, ") =", y
    print *, "exp(", x, ") =", exp(x)
    print *, "Error =", abs(y - exp(x))

    call pade%free()
end program pade_example
```

## API Reference

### Types

#### `bernstein_polynomial`

Represents a Bernstein polynomial on an interval [a, b].

**Methods:**
- `create(coeffs, a, b, error)` - Create from coefficients
- `free()` - Free memory
- `evaluate(x, error)` - Evaluate at point(s), returns scalar or array
- `derivative()` - Compute derivative (returns new polynomial)
- `integral()` - Compute antiderivative (returns new polynomial)
- `degree(error)` - Get polynomial degree
- `add(other)` - Add polynomials
- `subtract(other)` - Subtract polynomials
- `multiply(other)` - Multiply polynomials
- `scale(scalar)` - Multiply by scalar

#### `rational_function`

Represents a rational function P(x)/Q(x).

**Methods:**
- `create(numerator, denominator, error)` - Create from polynomials
- `free()` - Free memory
- `evaluate(x, error)` - Evaluate at point(s)
- `derivative()` - Compute derivative (returns new rational)

#### `pade_approximant`

Represents a Padé approximant [n/m].

**Methods:**
- `create_from_series(coeffs, n, m, a, b, error)` - From power series
- `free()` - Free memory
- `evaluate(x, error)` - Evaluate at point(s)

#### `gelfgren_error`

Error reporting type.

**Fields:**
- `code` (integer) - Error code (GELFGREN_SUCCESS = 0)
- `message` (character(256)) - Error description

**Error Codes:**
- `GELFGREN_SUCCESS` = 0
- `GELFGREN_NULL_POINTER` = -1
- `GELFGREN_INVALID_INTERVAL` = -2
- `GELFGREN_EMPTY_COEFFICIENTS` = -3
- `GELFGREN_SINGULAR_MATRIX` = -4
- `GELFGREN_DIVISION_BY_ZERO` = -5
- `GELFGREN_POLE` = -6

## Error Handling

All methods accept an optional `error` argument of type `gelfgren_error`:

```fortran
type(gelfgren_error) :: error

y = poly%evaluate(0.5d0, error)

if (error%code /= GELFGREN_SUCCESS) then
    print *, "Error:", error%code, trim(error%message)
    stop 1
end if
```

If `error` is omitted, errors are not reported (use with caution).

## Memory Management

All types require explicit cleanup with `call obj%free()`:

```fortran
type(bernstein_polynomial) :: poly

call poly%create([1.0d0, 2.0d0], 0.0d0, 1.0d0)
! ... use poly ...
call poly%free()  ! Important: free memory
```

Failure to call `free()` will result in memory leaks.

## Compiler Compatibility

Tested with:
- **GFortran** 8.0+ (recommended)
- **Intel Fortran** 19.0+
- **NAG Fortran** 7.0+

All compilers must support ISO_C_BINDING (Fortran 2003 standard).

## Performance

The Fortran bindings have minimal overhead:
- Single function call: ~5 ns FFI overhead
- Vectorized operations: overhead amortized across array elements
- No data copying for C-compatible arrays

For best performance:
- Use vectorized evaluation for large arrays
- Reuse polynomial objects when possible
- Compile with optimization flags (`-O2` or `-O3`)

## Examples

Complete working examples are in the `example/` directory:
- `basic_usage.f90` - Comprehensive demonstration of all features

Run with:
```bash
fpm run basic_usage
# or
make run-example
```

## Testing

Run the test suite:

```bash
fpm test
# or
make run-test
```

Tests cover:
- Polynomial creation and evaluation
- Derivatives and integrals
- Arithmetic operations
- Rational functions
- Padé approximants
- Error handling

## Troubleshooting

### Library Not Found

If you get "cannot open shared object file" errors:

```bash
export LD_LIBRARY_PATH=/path/to/gelfgren/target/release:$LD_LIBRARY_PATH
# On macOS:
export DYLD_LIBRARY_PATH=/path/to/gelfgren/target/release:$DYLD_LIBRARY_PATH
```

### Module File Issues

Ensure module files are in the include path:

```bash
gfortran -I/path/to/modules your_program.f90
```

### Linking Issues

Explicitly link the Gelfgren library:

```bash
gfortran your_program.f90 -L/path/to/lib -lgelfgren
```

## Contributing

Contributions are welcome! Please see the main [CONTRIBUTING.md](../../CONTRIBUTING.md) for guidelines.

## License

Licensed under either of MIT or Apache-2.0 at your option.

## See Also

- [Main Gelfgren Documentation](https://yourusername.github.io/gelfgren)
- [Rust API Documentation](https://docs.rs/gelfgren-core)
- [C API Header](../../include/gelfgren.h)
- [Mathematical Background](../../docs/mathematics.md)
