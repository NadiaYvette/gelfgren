;;;; CFFI bindings for Gelfgren C library

(in-package #:gelfgren)

;;; Load library
(define-foreign-library libgelfgren
  (:unix (:or "libgelfgren.so" "libgelfgren.dylib"))
  (t (:default "libgelfgren")))

(use-foreign-library libgelfgren)

;;; Error handling
(defcfun ("gelfgren_last_error_message" %gelfgren-last-error-message) :string)

;;; Bernstein Polynomial FFI
(defcfun ("gelfgren_bernstein_create" %bernstein-create) :pointer
  (coeffs :pointer)
  (degree :size)
  (a :double)
  (b :double))

(defcfun ("gelfgren_bernstein_free" %bernstein-free) :void
  (poly :pointer))

(defcfun ("gelfgren_bernstein_evaluate" %bernstein-evaluate) :int
  (poly :pointer)
  (x :double)
  (result :pointer))

(defcfun ("gelfgren_bernstein_derivative" %bernstein-derivative) :pointer
  (poly :pointer))

(defcfun ("gelfgren_bernstein_integral" %bernstein-integral) :pointer
  (poly :pointer))

(defcfun ("gelfgren_bernstein_degree" %bernstein-degree) :int
  (poly :pointer)
  (degree :pointer))

(defcfun ("gelfgren_bernstein_add" %bernstein-add) :pointer
  (p1 :pointer)
  (p2 :pointer))

(defcfun ("gelfgren_bernstein_scale" %bernstein-scale) :pointer
  (scalar :double)
  (poly :pointer))

;;; Rational Function FFI
(defcfun ("gelfgren_rational_create" %rational-create) :pointer
  (num :pointer)
  (den :pointer))

(defcfun ("gelfgren_rational_free" %rational-free) :void
  (rat :pointer))

(defcfun ("gelfgren_rational_evaluate" %rational-evaluate) :int
  (rat :pointer)
  (x :double)
  (result :pointer))

(defcfun ("gelfgren_rational_derivative" %rational-derivative) :pointer
  (rat :pointer))

;;; Padé Approximant FFI
(defcfun ("gelfgren_pade_from_series" %pade-from-series) :pointer
  (coeffs :pointer)
  (len :size)
  (n :int)
  (m :int)
  (a :double)
  (b :double))

(defcfun ("gelfgren_pade_free" %pade-free) :void
  (pade :pointer))

(defcfun ("gelfgren_pade_evaluate" %pade-evaluate) :int
  (pade :pointer)
  (x :double)
  (result :pointer))
