#!/usr/bin/env guile
!#

(add-to-load-path "../src")
(use-modules (gelfgren))

(define (main)
  (display "Gelfgren Guile Scheme Bindings Example\n")
  (display (make-string 60 #\=))
  (newline)
  (newline)

  ;; Example 1: Bernstein Polynomial
  (display "1. Bernstein Polynomial\n")
  (display (make-string 60 #\-))
  (newline)

  (let ((poly (make-bernstein-polynomial '(1.0 2.0 6.0) 0.0 1.0)))
    (let ((deg (degree-bernstein poly)))
      (format #t "Created polynomial with degree: ~a\n" deg))

    (let ((y (evaluate-bernstein poly 0.5)))
      (format #t "f(0.5) = ~a\n" y))

    (newline)
    (display "Vectorized evaluation:\n")
    (for-each
     (lambda (x)
       (format #t "  f(~a) = ~a\n" x (evaluate-bernstein poly x)))
     '(0.0 0.25 0.5 0.75 1.0))

    (free-bernstein-polynomial poly))

  (newline)

  ;; Example 2: Rational Function
  (display "2. Rational Function\n")
  (display (make-string 60 #\-))
  (newline)

  (let* ((num (make-bernstein-polynomial '(0.0 1.0) 0.0 1.0))
         (den (make-bernstein-polynomial '(1.0 2.0) 0.0 1.0))
         (rat (make-rational-function num den)))
    (let* ((x 0.5)
           (y (evaluate-rational rat x))
           (expected (/ x (+ 1.0 x))))
      (format #t "R(0.5) = ~a\n" y)
      (format #t "Expected: ~a\n" expected))

    (free-rational-function rat)
    (free-bernstein-polynomial den)
    (free-bernstein-polynomial num))

  (newline)
  (display "All examples completed successfully!\n"))

(main)
