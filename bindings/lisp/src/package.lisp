;;;; Package definition for Gelfgren

(defpackage #:gelfgren
  (:use #:cl #:cffi)
  (:export
   ;; Error handling
   #:gelfgren-error
   #:gelfgren-error-message

   ;; Bernstein Polynomial
   #:bernstein-polynomial
   #:make-bernstein-polynomial
   #:free-bernstein-polynomial
   #:with-bernstein-polynomial
   #:evaluate-bernstein
   #:derivative-bernstein
   #:integral-bernstein
   #:degree-bernstein
   #:add-bernstein
   #:scale-bernstein

   ;; Rational Function
   #:rational-function
   #:make-rational-function
   #:free-rational-function
   #:with-rational-function
   #:evaluate-rational
   #:derivative-rational

   ;; Padé Approximant
   #:pade-approximant
   #:make-pade-approximant
   #:free-pade-approximant
   #:with-pade-approximant
   #:evaluate-pade))
