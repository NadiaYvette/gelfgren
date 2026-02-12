# Gelfgren Guile Scheme Bindings

Guile Scheme bindings for Gelfgren using FFI.

## Requirements

- GNU Guile 3.0+
- Gelfgren C library

## Usage

```scheme
(use-modules (gelfgren))

(let ((poly (make-bernstein-polynomial '(1.0 2.0 6.0) 0.0 1.0)))
  (let ((y (evaluate-bernstein poly 0.5)))
    (display y))
  (free-bernstein-polynomial poly))
```

## Running Examples

```bash
export LD_LIBRARY_PATH=../../target/release:$LD_LIBRARY_PATH
guile -L src examples/basic-usage.scm
```

## License

MIT or Apache-2.0
