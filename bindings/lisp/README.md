# Gelfgren Common Lisp Bindings

Common Lisp bindings for Gelfgren using CFFI.

## Requirements

- Common Lisp implementation (SBCL, CCL, etc.)
- CFFI library
- Gelfgren C library

## Installation

```lisp
(ql:quickload :gelfgren)
```

## Usage

```lisp
(use-package :gelfgren)

;; Create polynomial
(with-bernstein-polynomial (poly '(1.0d0 2.0d0 6.0d0) 0.0d0 1.0d0)
  (let ((y (evaluate-bernstein poly 0.5d0)))
    (format t "f(0.5) = ~F~%" y)))
```

## API

- `make-bernstein-polynomial`, `evaluate-bernstein`, `derivative-bernstein`
- `make-rational-function`, `evaluate-rational`
- `make-pade-approximant`, `evaluate-pade`
- `with-*` macros for automatic cleanup

## License

MIT or Apache-2.0
