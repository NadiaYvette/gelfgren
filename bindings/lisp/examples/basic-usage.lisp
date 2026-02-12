;;;; Basic usage example for Gelfgren Common Lisp bindings

(ql:quickload :gelfgren)

(defun main ()
  (format t "Gelfgren Common Lisp Bindings Example~%")
  (format t "~60~~%~%")

  ;; Example 1: Bernstein Polynomial
  (format t "1. Bernstein Polynomial~%")
  (format t "~60~-~%")

  (with-bernstein-polynomial (poly '(1.0d0 2.0d0 6.0d0) 0.0d0 1.0d0)
    (let ((deg (degree-bernstein poly)))
      (format t "Created polynomial with degree: ~D~%" deg))

    (let ((y (evaluate-bernstein poly 0.5d0)))
      (format t "f(0.5) = ~F~%" y))

    ;; Vectorized evaluation
    (format t "~%Vectorized evaluation:~%")
    (dolist (x '(0.0d0 0.25d0 0.5d0 0.75d0 1.0d0))
      (format t "  f(~F) = ~F~%" x (evaluate-bernstein poly x)))

    ;; Derivative
    (format t "~%")
    (let ((dpoly (derivative-bernstein poly)))
      (format t "Derivative at x=0.5: ~F~%" (evaluate-bernstein dpoly 0.5d0))
      (free-bernstein-polynomial dpoly))

    ;; Integral
    (let ((ipoly (integral-bernstein poly)))
      (let ((area (- (evaluate-bernstein ipoly 1.0d0)
                     (evaluate-bernstein ipoly 0.0d0))))
        (format t "Integral from 0 to 1: ~F~%" area))
      (free-bernstein-polynomial ipoly)))

  (format t "~%")

  ;; Example 2: Rational Function
  (format t "2. Rational Function~%")
  (format t "~60~-~%")

  (with-bernstein-polynomial (num '(0.0d0 1.0d0) 0.0d0 1.0d0)
    (with-bernstein-polynomial (den '(1.0d0 2.0d0) 0.0d0 1.0d0)
      (with-rational-function (rat num den)
        (let* ((x 0.5d0)
               (y (evaluate-rational rat x))
               (expected (/ x (+ 1.0d0 x))))
          (format t "R(0.5) = ~F~%" y)
          (format t "Expected: ~F~%" expected))

        ;; Vectorized
        (format t "~%Vectorized evaluation:~%")
        (dolist (x '(0.0d0 0.25d0 0.5d0 0.75d0 1.0d0))
          (let ((y (evaluate-rational rat x))
                (expected (/ x (+ 1.0d0 x))))
            (format t "  R(~F) = ~F (expected: ~F)~%" x y expected))))))

  (format t "~%")

  ;; Example 3: Padé Approximant
  (format t "3. Padé Approximant~%")
  (format t "~60~-~%")

  (with-pade-approximant (pade '(1.0d0 1.0d0 0.5d0 0.16666666666666666d0 0.041666666666666664d0)
                                2 2 -1.0d0 1.0d0)
    (format t "Approximation of exp(x):~%")
    (dolist (x '(-0.5d0 -0.25d0 0.0d0 0.25d0 0.5d0))
      (let* ((y (evaluate-pade pade x))
             (exact (exp x))
             (error (abs (- y exact))))
        (format t "  x=~F: Padé=~F, exp(x)=~F, error=~F~%" x y exact error))))

  (format t "~%")
  (format t "All examples completed successfully!~%"))

(main)
