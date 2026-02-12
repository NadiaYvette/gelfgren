;;;; ASDF system definition for Gelfgren

(asdf:defsystem #:gelfgren
  :description "Common Lisp bindings for Gelfgren piecewise rational interpolation library"
  :author "Nadia Chambers <nadia.chambers@iohk.io>"
  :license "MIT OR Apache-2.0"
  :version "0.1.0"
  :serial t
  :depends-on (#:cffi)
  :components ((:module "src"
                :components
                ((:file "package")
                 (:file "bindings")
                 (:file "gelfgren"))))
  :in-order-to ((test-op (test-op #:gelfgren/test))))

(asdf:defsystem #:gelfgren/test
  :depends-on (#:gelfgren #:fiveam)
  :components ((:module "test"
                :components
                ((:file "test-gelfgren"))))
  :perform (test-op (o c) (symbol-call :fiveam :run! :gelfgren-suite)))
