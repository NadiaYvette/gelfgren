;;;; High-level Gelfgren API

(in-package #:gelfgren)

;;; Error handling
(define-condition gelfgren-error (error)
  ((message :initarg :message :reader gelfgren-error-message))
  (:report (lambda (condition stream)
             (format stream "Gelfgren error: ~A" (gelfgren-error-message condition)))))

(defun check-error (code operation)
  (unless (zerop code)
    (error 'gelfgren-error :message (format nil "~A: ~A" operation (%gelfgren-last-error-message)))))

(defun check-null (ptr operation)
  (when (null-pointer-p ptr)
    (error 'gelfgren-error :message (format nil "~A: ~A" operation (%gelfgren-last-error-message))))
  ptr)

;;; Bernstein Polynomial
(defclass bernstein-polynomial ()
  ((pointer :initarg :pointer :accessor bernstein-pointer)
   (a :initarg :a :reader bernstein-a)
   (b :initarg :b :reader bernstein-b)))

(defmethod print-object ((obj bernstein-polynomial) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "[~A, ~A]" (bernstein-a obj) (bernstein-b obj))))

(defun make-bernstein-polynomial (coeffs a b)
  "Create a Bernstein polynomial from coefficients on interval [A, B]."
  (let* ((degree (1- (length coeffs)))
         (coeffs-array (foreign-alloc :double :initial-contents coeffs))
         (ptr (%bernstein-create coeffs-array degree (coerce a 'double-float) (coerce b 'double-float))))
    (foreign-free coeffs-array)
    (check-null ptr "create-bernstein")
    (let ((poly (make-instance 'bernstein-polynomial :pointer ptr :a a :b b)))
      (trivial-garbage:finalize poly (lambda () (%bernstein-free ptr)))
      poly)))

(defun free-bernstein-polynomial (poly)
  "Explicitly free a Bernstein polynomial."
  (trivial-garbage:cancel-finalization poly)
  (%bernstein-free (bernstein-pointer poly))
  (setf (bernstein-pointer poly) (null-pointer)))

(defmacro with-bernstein-polynomial ((var coeffs a b) &body body)
  "Create polynomial, execute body, ensure cleanup."
  `(let ((,var (make-bernstein-polynomial ,coeffs ,a ,b)))
     (unwind-protect (progn ,@body)
       (free-bernstein-polynomial ,var))))

(defun evaluate-bernstein (poly x)
  "Evaluate Bernstein polynomial at point X."
  (with-foreign-object (result :double)
    (let ((code (%bernstein-evaluate (bernstein-pointer poly) (coerce x 'double-float) result)))
      (check-error code "evaluate")
      (mem-ref result :double))))

(defun derivative-bernstein (poly)
  "Compute derivative of Bernstein polynomial."
  (let ((ptr (%bernstein-derivative (bernstein-pointer poly))))
    (check-null ptr "derivative")
    (let ((dpoly (make-instance 'bernstein-polynomial
                                :pointer ptr
                                :a (bernstein-a poly)
                                :b (bernstein-b poly))))
      (trivial-garbage:finalize dpoly (lambda () (%bernstein-free ptr)))
      dpoly)))

(defun integral-bernstein (poly)
  "Compute antiderivative of Bernstein polynomial."
  (let ((ptr (%bernstein-integral (bernstein-pointer poly))))
    (check-null ptr "integral")
    (let ((ipoly (make-instance 'bernstein-polynomial
                                :pointer ptr
                                :a (bernstein-a poly)
                                :b (bernstein-b poly))))
      (trivial-garbage:finalize ipoly (lambda () (%bernstein-free ptr)))
      ipoly)))

(defun degree-bernstein (poly)
  "Get degree of Bernstein polynomial."
  (with-foreign-object (deg :int)
    (let ((code (%bernstein-degree (bernstein-pointer poly) deg)))
      (check-error code "degree")
      (mem-ref deg :int))))

(defun add-bernstein (p1 p2)
  "Add two Bernstein polynomials."
  (let ((ptr (%bernstein-add (bernstein-pointer p1) (bernstein-pointer p2))))
    (check-null ptr "add")
    (let ((sum (make-instance 'bernstein-polynomial
                              :pointer ptr
                              :a (bernstein-a p1)
                              :b (bernstein-b p1))))
      (trivial-garbage:finalize sum (lambda () (%bernstein-free ptr)))
      sum)))

(defun scale-bernstein (scalar poly)
  "Scale Bernstein polynomial by scalar."
  (let ((ptr (%bernstein-scale (coerce scalar 'double-float) (bernstein-pointer poly))))
    (check-null ptr "scale")
    (let ((scaled (make-instance 'bernstein-polynomial
                                 :pointer ptr
                                 :a (bernstein-a poly)
                                 :b (bernstein-b poly))))
      (trivial-garbage:finalize scaled (lambda () (%bernstein-free ptr)))
      scaled)))

;;; Rational Function
(defclass rational-function ()
  ((pointer :initarg :pointer :accessor rational-pointer)))

(defmethod print-object ((obj rational-function) stream)
  (print-unreadable-object (obj stream :type t)))

(defun make-rational-function (num den)
  "Create rational function from numerator and denominator."
  (let ((ptr (%rational-create (bernstein-pointer num) (bernstein-pointer den))))
    (check-null ptr "create-rational")
    (let ((rat (make-instance 'rational-function :pointer ptr)))
      (trivial-garbage:finalize rat (lambda () (%rational-free ptr)))
      rat)))

(defun free-rational-function (rat)
  "Explicitly free rational function."
  (trivial-garbage:cancel-finalization rat)
  (%rational-free (rational-pointer rat))
  (setf (rational-pointer rat) (null-pointer)))

(defmacro with-rational-function ((var num den) &body body)
  `(let ((,var (make-rational-function ,num ,den)))
     (unwind-protect (progn ,@body)
       (free-rational-function ,var))))

(defun evaluate-rational (rat x)
  "Evaluate rational function at point X."
  (with-foreign-object (result :double)
    (let ((code (%rational-evaluate (rational-pointer rat) (coerce x 'double-float) result)))
      (check-error code "evaluate-rational")
      (mem-ref result :double))))

;;; Padé Approximant
(defclass pade-approximant ()
  ((pointer :initarg :pointer :accessor pade-pointer)))

(defmethod print-object ((obj pade-approximant) stream)
  (print-unreadable-object (obj stream :type t)))

(defun make-pade-approximant (coeffs n m a b)
  "Create Padé approximant [N/M] from power series coefficients."
  (let* ((len (length coeffs))
         (coeffs-array (foreign-alloc :double :initial-contents coeffs))
         (ptr (%pade-from-series coeffs-array len n m (coerce a 'double-float) (coerce b 'double-float))))
    (foreign-free coeffs-array)
    (check-null ptr "create-pade")
    (let ((pade (make-instance 'pade-approximant :pointer ptr)))
      (trivial-garbage:finalize pade (lambda () (%pade-free ptr)))
      pade)))

(defun free-pade-approximant (pade)
  "Explicitly free Padé approximant."
  (trivial-garbage:cancel-finalization pade)
  (%pade-free (pade-pointer pade))
  (setf (pade-pointer pade) (null-pointer)))

(defmacro with-pade-approximant ((var coeffs n m a b) &body body)
  `(let ((,var (make-pade-approximant ,coeffs ,n ,m ,a ,b)))
     (unwind-protect (progn ,@body)
       (free-pade-approximant ,var))))

(defun evaluate-pade (pade x)
  "Evaluate Padé approximant at point X."
  (with-foreign-object (result :double)
    (let ((code (%pade-evaluate (pade-pointer pade) (coerce x 'double-float) result)))
      (check-error code "evaluate-pade")
      (mem-ref result :double))))
