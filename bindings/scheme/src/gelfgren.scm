;;; Guile Scheme bindings for Gelfgren

(define-module (gelfgren)
  #:use-module (system foreign)
  #:use-module (rnrs bytevectors)
  #:export (make-bernstein-polynomial
            free-bernstein-polynomial
            evaluate-bernstein
            derivative-bernstein
            integral-bernstein
            degree-bernstein
            add-bernstein
            scale-bernstein
            make-rational-function
            free-rational-function
            evaluate-rational
            make-pade-approximant
            free-pade-approximant
            evaluate-pade))

;;; Load library
(define libgelfgren
  (dynamic-link "libgelfgren"))

(define %last-error
  (pointer->procedure '* (dynamic-func "gelfgren_last_error_message" libgelfgren) '()))

;;; Error handling
(define (check-error code operation)
  (unless (= code 0)
    (error (format #f "~a: ~a" operation (pointer->string (%last-error))))))

(define (check-null ptr operation)
  (when (null-pointer? ptr)
    (error (format #f "~a: ~a" operation (pointer->string (%last-error)))))
  ptr)

;;; Bernstein Polynomial FFI
(define %bernstein-create
  (pointer->procedure '* (dynamic-func "gelfgren_bernstein_create" libgelfgren)
                      (list '* size_t double double)))

(define %bernstein-free
  (pointer->procedure void (dynamic-func "gelfgren_bernstein_free" libgelfgren) (list '*)))

(define %bernstein-evaluate
  (pointer->procedure int (dynamic-func "gelfgren_bernstein_evaluate" libgelfgren)
                      (list '* double '*)))

(define %bernstein-derivative
  (pointer->procedure '* (dynamic-func "gelfgren_bernstein_derivative" libgelfgren) (list '*)))

(define %bernstein-integral
  (pointer->procedure '* (dynamic-func "gelfgren_bernstein_integral" libgelfgren) (list '*)))

(define %bernstein-degree
  (pointer->procedure int (dynamic-func "gelfgren_bernstein_degree" libgelfgren) (list '* '*)))

(define %bernstein-add
  (pointer->procedure '* (dynamic-func "gelfgren_bernstein_add" libgelfgren) (list '* '*)))

(define %bernstein-scale
  (pointer->procedure '* (dynamic-func "gelfgren_bernstein_scale" libgelfgren) (list double '*)))

;;; High-level API
(define (make-bernstein-polynomial coeffs a b)
  "Create Bernstein polynomial from coefficient list on interval [a, b]"
  (let* ((degree (- (length coeffs) 1))
         (bv (make-f64vector (length coeffs))))
    (let loop ((i 0) (cs coeffs))
      (unless (null? cs)
        (f64vector-set! bv i (exact->inexact (car cs)))
        (loop (+ i 1) (cdr cs))))
    (let ((ptr (%bernstein-create (bytevector->pointer bv) degree a b)))
      (check-null ptr "create-bernstein"))))

(define (free-bernstein-polynomial poly)
  (%bernstein-free poly))

(define (evaluate-bernstein poly x)
  (let ((result (make-f64vector 1)))
    (let ((code (%bernstein-evaluate poly x (bytevector->pointer result))))
      (check-error code "evaluate")
      (f64vector-ref result 0))))

(define (derivative-bernstein poly)
  (check-null (%bernstein-derivative poly) "derivative"))

(define (integral-bernstein poly)
  (check-null (%bernstein-integral poly) "integral"))

(define (degree-bernstein poly)
  (let ((deg (make-s32vector 1)))
    (let ((code (%bernstein-degree poly (bytevector->pointer deg))))
      (check-error code "degree")
      (s32vector-ref deg 0))))

(define (add-bernstein p1 p2)
  (check-null (%bernstein-add p1 p2) "add"))

(define (scale-bernstein scalar poly)
  (check-null (%bernstein-scale scalar poly) "scale"))

;;; Rational Function FFI
(define %rational-create
  (pointer->procedure '* (dynamic-func "gelfgren_rational_create" libgelfgren) (list '* '*)))

(define %rational-free
  (pointer->procedure void (dynamic-func "gelfgren_rational_free" libgelfgren) (list '*)))

(define %rational-evaluate
  (pointer->procedure int (dynamic-func "gelfgren_rational_evaluate" libgelfgren) (list '* double '*)))

(define (make-rational-function num den)
  (check-null (%rational-create num den) "create-rational"))

(define (free-rational-function rat)
  (%rational-free rat))

(define (evaluate-rational rat x)
  (let ((result (make-f64vector 1)))
    (let ((code (%rational-evaluate rat x (bytevector->pointer result))))
      (check-error code "evaluate-rational")
      (f64vector-ref result 0))))

;;; Padé Approximant FFI
(define %pade-from-series
  (pointer->procedure '* (dynamic-func "gelfgren_pade_from_series" libgelfgren)
                      (list '* size_t int int double double)))

(define %pade-free
  (pointer->procedure void (dynamic-func "gelfgren_pade_free" libgelfgren) (list '*)))

(define %pade-evaluate
  (pointer->procedure int (dynamic-func "gelfgren_pade_evaluate" libgelfgren) (list '* double '*)))

(define (make-pade-approximant coeffs n m a b)
  (let* ((len (length coeffs))
         (bv (make-f64vector len)))
    (let loop ((i 0) (cs coeffs))
      (unless (null? cs)
        (f64vector-set! bv i (exact->inexact (car cs)))
        (loop (+ i 1) (cdr cs))))
    (check-null (%pade-from-series (bytevector->pointer bv) len n m a b) "create-pade")))

(define (free-pade-approximant pade)
  (%pade-free pade))

(define (evaluate-pade pade x)
  (let ((result (make-f64vector 1)))
    (let ((code (%pade-evaluate pade x (bytevector->pointer result))))
      (check-error code "evaluate-pade")
      (f64vector-ref result 0))))
