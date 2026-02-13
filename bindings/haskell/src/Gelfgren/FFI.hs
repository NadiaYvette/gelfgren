{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE EmptyDataDecls #-}

-- | Low-level FFI bindings to the Gelfgren C library
module Gelfgren.FFI
    ( -- * Opaque pointer types
      GelfgrenBernstein
    , GelfgrenRational
    , GelfgrenPade
    , GelfgrenErrorCode(..)

      -- * Bernstein polynomial functions
    , c_bernstein_create
    , c_bernstein_free
    , c_bernstein_evaluate
    , c_bernstein_derivative
    , c_bernstein_integral
    , c_bernstein_degree
    , c_bernstein_add
    , c_bernstein_subtract
    , c_bernstein_multiply
    , c_bernstein_scale

      -- * Rational function functions
    , c_rational_create
    , c_rational_free
    , c_rational_evaluate
    , c_rational_derivative

      -- * Padé approximant functions
    , c_pade_from_power_series
    , c_pade_free
    , c_pade_evaluate

      -- * Error handling
    , c_last_error_message
    ) where

import Foreign.C.Types
import Foreign.C.String
import Foreign.Ptr

-- | Opaque pointer to Bernstein polynomial
data GelfgrenBernstein

-- | Opaque pointer to rational function
data GelfgrenRational

-- | Opaque pointer to Padé approximant
data GelfgrenPade

-- | Error codes returned by C functions
newtype GelfgrenErrorCode = GelfgrenErrorCode CInt
    deriving (Eq, Show)

pattern GELFGREN_SUCCESS :: GelfgrenErrorCode
pattern GELFGREN_SUCCESS = GelfgrenErrorCode 0

pattern GELFGREN_NULL_POINTER :: GelfgrenErrorCode
pattern GELFGREN_NULL_POINTER = GelfgrenErrorCode (-1)

pattern GELFGREN_INVALID_INTERVAL :: GelfgrenErrorCode
pattern GELFGREN_INVALID_INTERVAL = GelfgrenErrorCode (-2)

pattern GELFGREN_EMPTY_COEFFICIENTS :: GelfgrenErrorCode
pattern GELFGREN_EMPTY_COEFFICIENTS = GelfgrenErrorCode (-3)

pattern GELFGREN_SINGULAR_MATRIX :: GelfgrenErrorCode
pattern GELFGREN_SINGULAR_MATRIX = GelfgrenErrorCode (-4)

pattern GELFGREN_DIVISION_BY_ZERO :: GelfgrenErrorCode
pattern GELFGREN_DIVISION_BY_ZERO = GelfgrenErrorCode (-5)

pattern GELFGREN_POLE :: GelfgrenErrorCode
pattern GELFGREN_POLE = GelfgrenErrorCode (-6)

-- Bernstein polynomial functions

foreign import ccall unsafe "gelfgren_bernstein_create"
    c_bernstein_create :: Ptr CDouble  -- ^ coefficients
                       -> CSize        -- ^ degree
                       -> CDouble      -- ^ a
                       -> CDouble      -- ^ b
                       -> IO (Ptr GelfgrenBernstein)

foreign import ccall unsafe "gelfgren_bernstein_free"
    c_bernstein_free :: Ptr GelfgrenBernstein -> IO ()

foreign import ccall unsafe "gelfgren_bernstein_evaluate"
    c_bernstein_evaluate :: Ptr GelfgrenBernstein  -- ^ polynomial
                         -> CDouble                -- ^ x
                         -> Ptr CDouble            -- ^ result
                         -> IO GelfgrenErrorCode

foreign import ccall unsafe "gelfgren_bernstein_derivative"
    c_bernstein_derivative :: Ptr GelfgrenBernstein -> IO (Ptr GelfgrenBernstein)

foreign import ccall unsafe "gelfgren_bernstein_integral"
    c_bernstein_integral :: Ptr GelfgrenBernstein -> IO (Ptr GelfgrenBernstein)

foreign import ccall unsafe "gelfgren_bernstein_degree"
    c_bernstein_degree :: Ptr GelfgrenBernstein  -- ^ polynomial
                       -> Ptr CSize              -- ^ degree (output)
                       -> IO GelfgrenErrorCode

foreign import ccall unsafe "gelfgren_bernstein_add"
    c_bernstein_add :: Ptr GelfgrenBernstein  -- ^ poly1
                    -> Ptr GelfgrenBernstein  -- ^ poly2
                    -> IO (Ptr GelfgrenBernstein)

foreign import ccall unsafe "gelfgren_bernstein_subtract"
    c_bernstein_subtract :: Ptr GelfgrenBernstein  -- ^ poly1
                         -> Ptr GelfgrenBernstein  -- ^ poly2
                         -> IO (Ptr GelfgrenBernstein)

foreign import ccall unsafe "gelfgren_bernstein_multiply"
    c_bernstein_multiply :: Ptr GelfgrenBernstein  -- ^ poly1
                         -> Ptr GelfgrenBernstein  -- ^ poly2
                         -> IO (Ptr GelfgrenBernstein)

foreign import ccall unsafe "gelfgren_bernstein_scale"
    c_bernstein_scale :: Ptr GelfgrenBernstein  -- ^ polynomial
                      -> CDouble                -- ^ scalar
                      -> IO (Ptr GelfgrenBernstein)

-- Rational function functions

foreign import ccall unsafe "gelfgren_rational_create"
    c_rational_create :: Ptr GelfgrenBernstein  -- ^ numerator
                      -> Ptr GelfgrenBernstein  -- ^ denominator
                      -> IO (Ptr GelfgrenRational)

foreign import ccall unsafe "gelfgren_rational_free"
    c_rational_free :: Ptr GelfgrenRational -> IO ()

foreign import ccall unsafe "gelfgren_rational_evaluate"
    c_rational_evaluate :: Ptr GelfgrenRational  -- ^ rational
                        -> CDouble               -- ^ x
                        -> Ptr CDouble           -- ^ result
                        -> IO GelfgrenErrorCode

foreign import ccall unsafe "gelfgren_rational_derivative"
    c_rational_derivative :: Ptr GelfgrenRational -> IO (Ptr GelfgrenRational)

-- Padé approximant functions

foreign import ccall unsafe "gelfgren_pade_from_power_series"
    c_pade_from_power_series :: Ptr CDouble  -- ^ coefficients
                             -> CSize        -- ^ n_coeffs
                             -> CSize        -- ^ n
                             -> CSize        -- ^ m
                             -> CDouble      -- ^ a
                             -> CDouble      -- ^ b
                             -> IO (Ptr GelfgrenPade)

foreign import ccall unsafe "gelfgren_pade_free"
    c_pade_free :: Ptr GelfgrenPade -> IO ()

foreign import ccall unsafe "gelfgren_pade_evaluate"
    c_pade_evaluate :: Ptr GelfgrenPade  -- ^ pade
                    -> CDouble           -- ^ x
                    -> Ptr CDouble       -- ^ result
                    -> IO GelfgrenErrorCode

-- Error handling

foreign import ccall unsafe "gelfgren_last_error_message"
    c_last_error_message :: IO CString
