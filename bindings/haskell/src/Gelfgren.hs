{-# LANGUAGE GeneralizedNewtypeDeriving #-}

-- | High-level Haskell interface to Gelfgren numerical computing library
--
-- This module provides safe, idiomatic Haskell bindings to the Gelfgren library,
-- which implements piecewise rational interpolation methods.
--
-- Example usage:
--
-- @
-- import Gelfgren
--
-- main :: IO ()
-- main = do
--     -- Create a Bernstein polynomial
--     poly <- bernsteinPolynomial [1.0, 2.0, 6.0] 0.0 1.0
--
--     -- Evaluate at a point
--     y <- evaluate poly 0.5
--     print y
--
--     -- Compute derivative
--     dpoly <- derivative poly
--     dy <- evaluate dpoly 0.5
--     print dy
-- @
module Gelfgren
    ( -- * Types
      BernsteinPolynomial
    , RationalFunction
    , PadeApproximant
    , GelfgrenError(..)

      -- * Creating polynomials
    , bernsteinPolynomial
    , rationalFunction
    , padeFromSeries
    , padeFromDerivatives

      -- * Evaluation
    , evaluate
    , evaluateMany

      -- * Calculus
    , derivative
    , integral
    , degree

      -- * Arithmetic
    , add
    , subtract
    , multiply
    , scale

      -- * Utilities
    , withBernstein
    , withRational
    , withPade
    ) where

import Prelude hiding (subtract)
import Control.Exception (Exception, throwIO, bracket)
import Control.Monad (when)
import Data.Typeable (Typeable)
import Foreign.C.Types
import Foreign.C.String
import Foreign.Ptr
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc
import Foreign.Marshal.Array
import Foreign.Storable

import qualified Gelfgren.FFI as FFI

-- | Bernstein polynomial on an interval [a, b]
newtype BernsteinPolynomial = BernsteinPolynomial (ForeignPtr FFI.GelfgrenBernstein)
    deriving (Eq)

-- | Rational function P(x)/Q(x)
newtype RationalFunction = RationalFunction (ForeignPtr FFI.GelfgrenRational)
    deriving (Eq)

-- | Padé approximant [n/m]
newtype PadeApproximant = PadeApproximant (ForeignPtr FFI.GelfgrenPade)
    deriving (Eq)

-- | Errors that can occur during Gelfgren operations
data GelfgrenError
    = NullPointerError String
    | InvalidIntervalError String
    | EmptyCoefficientsError String
    | SingularMatrixError String
    | DivisionByZeroError String
    | PoleError String
    | UnknownError String
    deriving (Show, Eq, Typeable)

instance Exception GelfgrenError

-- | Get last error message from C library
getLastErrorMessage :: IO String
getLastErrorMessage = do
    cstr <- FFI.c_last_error_message
    if cstr == nullPtr
        then return "Unknown error"
        else peekCString cstr

-- | Convert error code to Haskell error
errorCodeToError :: FFI.GelfgrenErrorCode -> IO GelfgrenError
errorCodeToError code = do
    msg <- getLastErrorMessage
    return $ case code of
        FFI.GELFGREN_NULL_POINTER -> NullPointerError msg
        FFI.GELFGREN_INVALID_INTERVAL -> InvalidIntervalError msg
        FFI.GELFGREN_EMPTY_COEFFICIENTS -> EmptyCoefficientsError msg
        FFI.GELFGREN_SINGULAR_MATRIX -> SingularMatrixError msg
        FFI.GELFGREN_DIVISION_BY_ZERO -> DivisionByZeroError msg
        FFI.GELFGREN_POLE -> PoleError msg
        _ -> UnknownError msg

-- | Check error code and throw exception if needed
checkError :: FFI.GelfgrenErrorCode -> IO ()
checkError code =
    when (code /= FFI.GELFGREN_SUCCESS) $ do
        err <- errorCodeToError code
        throwIO err

-- | Create a Bernstein polynomial from coefficients
--
-- > poly <- bernsteinPolynomial [1.0, 2.0, 6.0] 0.0 1.0
bernsteinPolynomial :: [Double]  -- ^ Coefficients (length n+1 for degree n)
                    -> Double    -- ^ Left endpoint a
                    -> Double    -- ^ Right endpoint b
                    -> IO BernsteinPolynomial
bernsteinPolynomial coeffs a b = do
    let n = length coeffs
    when (n == 0) $ throwIO $ EmptyCoefficientsError "Coefficient list is empty"

    withArray (map realToFrac coeffs) $ \coeffsPtr -> do
        ptr <- FFI.c_bernstein_create coeffsPtr (fromIntegral (n - 1)) (realToFrac a) (realToFrac b)

        when (ptr == nullPtr) $ do
            msg <- getLastErrorMessage
            throwIO $ NullPointerError msg

        fptr <- newForeignPtr_ ptr
        addForeignPtrFinalizer finalizerBernstein fptr
        return $ BernsteinPolynomial fptr

foreign import ccall unsafe "&gelfgren_bernstein_free"
    finalizerBernstein :: FunPtr (Ptr FFI.GelfgrenBernstein -> IO ())

-- | Create a rational function from numerator and denominator polynomials
--
-- > rat <- rationalFunction num den
rationalFunction :: BernsteinPolynomial  -- ^ Numerator
                 -> BernsteinPolynomial  -- ^ Denominator
                 -> IO RationalFunction
rationalFunction (BernsteinPolynomial numFptr) (BernsteinPolynomial denFptr) =
    withForeignPtr numFptr $ \numPtr ->
    withForeignPtr denFptr $ \denPtr -> do
        ptr <- FFI.c_rational_create numPtr denPtr

        when (ptr == nullPtr) $ do
            msg <- getLastErrorMessage
            throwIO $ NullPointerError msg

        fptr <- newForeignPtr_ ptr
        addForeignPtrFinalizer finalizerRational fptr
        return $ RationalFunction fptr

foreign import ccall unsafe "&gelfgren_rational_free"
    finalizerRational :: FunPtr (Ptr FFI.GelfgrenRational -> IO ())

-- | Create a Padé approximant from power series coefficients
--
-- > pade <- padeFromSeries [1.0, 1.0, 0.5, 1.0/6.0] 2 1 (-1.0) 1.0
padeFromSeries :: [Double]  -- ^ Power series coefficients
               -> Int       -- ^ Numerator degree n
               -> Int       -- ^ Denominator degree m
               -> Double    -- ^ Left endpoint a
               -> Double    -- ^ Right endpoint b
               -> IO PadeApproximant
padeFromSeries coeffs n m a b = do
    let nCoeffs = length coeffs
    when (nCoeffs == 0) $ throwIO $ EmptyCoefficientsError "Coefficient list is empty"

    withArray (map realToFrac coeffs) $ \coeffsPtr -> do
        ptr <- FFI.c_pade_from_power_series
                coeffsPtr
                (fromIntegral nCoeffs)
                (fromIntegral n)
                (fromIntegral m)
                (realToFrac a)
                (realToFrac b)

        when (ptr == nullPtr) $ do
            msg <- getLastErrorMessage
            throwIO $ NullPointerError msg

        fptr <- newForeignPtr_ ptr
        addForeignPtrFinalizer finalizerPade fptr
        return $ PadeApproximant fptr

foreign import ccall unsafe "&gelfgren_pade_free"
    finalizerPade :: FunPtr (Ptr FFI.GelfgrenPade -> IO ())

-- | Placeholder for Padé from derivatives (not yet implemented in C API)
padeFromDerivatives :: [Double] -> Int -> Int -> Double -> Double -> Double -> IO PadeApproximant
padeFromDerivatives _ _ _ _ _ _ =
    throwIO $ UnknownError "padeFromDerivatives not yet implemented"

-- | Type class for evaluable objects
class Evaluable a where
    evaluateImpl :: a -> Double -> IO Double

instance Evaluable BernsteinPolynomial where
    evaluateImpl (BernsteinPolynomial fptr) x =
        withForeignPtr fptr $ \ptr ->
        alloca $ \resultPtr -> do
            code <- FFI.c_bernstein_evaluate ptr (realToFrac x) resultPtr
            checkError code
            result <- peek resultPtr
            return $ realToFrac result

instance Evaluable RationalFunction where
    evaluateImpl (RationalFunction fptr) x =
        withForeignPtr fptr $ \ptr ->
        alloca $ \resultPtr -> do
            code <- FFI.c_rational_evaluate ptr (realToFrac x) resultPtr
            checkError code
            result <- peek resultPtr
            return $ realToFrac result

instance Evaluable PadeApproximant where
    evaluateImpl (PadeApproximant fptr) x =
        withForeignPtr fptr $ \ptr ->
        alloca $ \resultPtr -> do
            code <- FFI.c_pade_evaluate ptr (realToFrac x) resultPtr
            checkError code
            result <- peek resultPtr
            return $ realToFrac result

-- | Evaluate a polynomial, rational function, or Padé approximant at a point
--
-- > y <- evaluate poly 0.5
evaluate :: Evaluable a => a -> Double -> IO Double
evaluate = evaluateImpl

-- | Evaluate at multiple points
--
-- > ys <- evaluateMany poly [0.0, 0.25, 0.5, 0.75, 1.0]
evaluateMany :: Evaluable a => a -> [Double] -> IO [Double]
evaluateMany obj xs = mapM (evaluate obj) xs

-- | Compute the derivative of a Bernstein polynomial
--
-- > dpoly <- derivative poly
derivative :: BernsteinPolynomial -> IO BernsteinPolynomial
derivative (BernsteinPolynomial fptr) =
    withForeignPtr fptr $ \ptr -> do
        dptr <- FFI.c_bernstein_derivative ptr
        when (dptr == nullPtr) $ do
            msg <- getLastErrorMessage
            throwIO $ NullPointerError msg

        dfptr <- newForeignPtr_ dptr
        addForeignPtrFinalizer finalizerBernstein dfptr
        return $ BernsteinPolynomial dfptr

-- | Compute the integral of a Bernstein polynomial
--
-- > ipoly <- integral poly
integral :: BernsteinPolynomial -> IO BernsteinPolynomial
integral (BernsteinPolynomial fptr) =
    withForeignPtr fptr $ \ptr -> do
        iptr <- FFI.c_bernstein_integral ptr
        when (iptr == nullPtr) $ do
            msg <- getLastErrorMessage
            throwIO $ NullPointerError msg

        ifptr <- newForeignPtr_ iptr
        addForeignPtrFinalizer finalizerBernstein ifptr
        return $ BernsteinPolynomial ifptr

-- | Get the degree of a Bernstein polynomial
--
-- > d <- degree poly
degree :: BernsteinPolynomial -> IO Int
degree (BernsteinPolynomial fptr) =
    withForeignPtr fptr $ \ptr ->
    alloca $ \degPtr -> do
        code <- FFI.c_bernstein_degree ptr degPtr
        checkError code
        deg <- peek degPtr
        return $ fromIntegral deg

-- | Add two Bernstein polynomials
--
-- > sum_poly <- add poly1 poly2
add :: BernsteinPolynomial -> BernsteinPolynomial -> IO BernsteinPolynomial
add (BernsteinPolynomial fptr1) (BernsteinPolynomial fptr2) =
    withForeignPtr fptr1 $ \ptr1 ->
    withForeignPtr fptr2 $ \ptr2 -> do
        sptr <- FFI.c_bernstein_add ptr1 ptr2
        when (sptr == nullPtr) $ do
            msg <- getLastErrorMessage
            throwIO $ NullPointerError msg

        sfptr <- newForeignPtr_ sptr
        addForeignPtrFinalizer finalizerBernstein sfptr
        return $ BernsteinPolynomial sfptr

-- | Subtract two Bernstein polynomials
--
-- > diff_poly <- subtract poly1 poly2
subtract :: BernsteinPolynomial -> BernsteinPolynomial -> IO BernsteinPolynomial
subtract (BernsteinPolynomial fptr1) (BernsteinPolynomial fptr2) =
    withForeignPtr fptr1 $ \ptr1 ->
    withForeignPtr fptr2 $ \ptr2 -> do
        dptr <- FFI.c_bernstein_subtract ptr1 ptr2
        when (dptr == nullPtr) $ do
            msg <- getLastErrorMessage
            throwIO $ NullPointerError msg

        dfptr <- newForeignPtr_ dptr
        addForeignPtrFinalizer finalizerBernstein dfptr
        return $ BernsteinPolynomial dfptr

-- | Multiply two Bernstein polynomials
--
-- > prod_poly <- multiply poly1 poly2
multiply :: BernsteinPolynomial -> BernsteinPolynomial -> IO BernsteinPolynomial
multiply (BernsteinPolynomial fptr1) (BernsteinPolynomial fptr2) =
    withForeignPtr fptr1 $ \ptr1 ->
    withForeignPtr fptr2 $ \ptr2 -> do
        pptr <- FFI.c_bernstein_multiply ptr1 ptr2
        when (pptr == nullPtr) $ do
            msg <- getLastErrorMessage
            throwIO $ NullPointerError msg

        pfptr <- newForeignPtr_ pptr
        addForeignPtrFinalizer finalizerBernstein pfptr
        return $ BernsteinPolynomial pfptr

-- | Scale a Bernstein polynomial by a constant
--
-- > scaled_poly <- scale 2.0 poly
scale :: Double -> BernsteinPolynomial -> IO BernsteinPolynomial
scale scalar (BernsteinPolynomial fptr) =
    withForeignPtr fptr $ \ptr -> do
        sptr <- FFI.c_bernstein_scale ptr (realToFrac scalar)
        when (sptr == nullPtr) $ do
            msg <- getLastErrorMessage
            throwIO $ NullPointerError msg

        sfptr <- newForeignPtr_ sptr
        addForeignPtrFinalizer finalizerBernstein sfptr
        return $ BernsteinPolynomial sfptr

-- | Bracket-style helper for Bernstein polynomials
--
-- > withBernstein [1, 2, 3] 0 1 $ \poly -> evaluate poly 0.5
withBernstein :: [Double] -> Double -> Double -> (BernsteinPolynomial -> IO a) -> IO a
withBernstein coeffs a b = bracket (bernsteinPolynomial coeffs a b) freeBernstein
  where
    freeBernstein (BernsteinPolynomial fptr) = finalizeForeignPtr fptr

-- | Bracket-style helper for rational functions
withRational :: BernsteinPolynomial -> BernsteinPolynomial -> (RationalFunction -> IO a) -> IO a
withRational num den = bracket (rationalFunction num den) freeRational
  where
    freeRational (RationalFunction fptr) = finalizeForeignPtr fptr

-- | Bracket-style helper for Padé approximants
withPade :: [Double] -> Int -> Int -> Double -> Double -> (PadeApproximant -> IO a) -> IO a
withPade coeffs n m a b = bracket (padeFromSeries coeffs n m a b) freePade
  where
    freePade (PadeApproximant fptr) = finalizeForeignPtr fptr
