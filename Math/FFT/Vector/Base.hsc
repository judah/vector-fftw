-- | A basic interface between Vectors and the fftw library.
module Math.FFT.Vector.Base(
            -- * Planners
            Planner(..),
            planOfType,
            PlanType(..),
            plan,
            run,
            -- * Plans
            Plan(),
            planInputSize,
            planOutputSize,
            execute,
            -- * Unsafe C stuff
            CFlags,
            CPlan,
            -- * Normalization helpers
            Scalable(..),
            modifyInput,
            modifyOutput,
            constMultOutput,
            multC,
            unsafeModify,
            ) where

import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as MS
import qualified Data.Vector.Generic as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as UM
import Control.Monad.Primitive (RealWorld)
import Control.Monad(forM_)
import Foreign (Storable, Ptr, unsafePerformIO, FunPtr,
                ForeignPtr, withForeignPtr, newForeignPtr)
import Foreign.C (CInt, CUInt)
import Data.Bits ( (.&.) )
import Data.Complex(Complex(..))
import Foreign.Storable.Complex()



#include <fftw3.h>

---------------------
-- Creating FFTW plans

-- First, the planner flags:
data PlanType = Estimate | Measure | Patient | Exhaustive
data Preservation = PreserveInput | DestroyInput

type CFlags = CUInt

-- | Marshal the planner flags for use by fftw.
planInitFlags :: PlanType -> Preservation -> CFlags
planInitFlags pt pr = planTypeInt .&. preservationInt
  where
    planTypeInt = case pt of
                    Estimate -> #const FFTW_ESTIMATE
                    Measure -> #const FFTW_MEASURE
                    Patient -> #const FFTW_PATIENT
                    Exhaustive -> #const FFTW_EXHAUSTIVE
    preservationInt = case pr of
                    PreserveInput -> #const FFTW_PRESERVE_INPUT
                    DestroyInput -> #const FFTW_DESTROY_INPUT

newtype CPlan = CPlan {unCPlan :: ForeignPtr CPlan}

withPlan :: CPlan -> (Ptr CPlan -> IO a) -> IO a
withPlan = withForeignPtr . unCPlan

foreign import ccall unsafe fftw_execute :: Ptr CPlan -> IO ()
foreign import ccall "&" fftw_destroy_plan :: FunPtr (Ptr CPlan -> IO ())

newPlan :: Ptr CPlan -> IO CPlan
newPlan = fmap CPlan . newForeignPtr fftw_destroy_plan

----------------------------------------
-- vector-fftw plans

data Plan a b = Plan {
                    planInput :: !(VS.MVector RealWorld a),
                    planOutput :: !(VS.MVector RealWorld b),
                    planExecute :: IO ()
                }

planInputSize :: Storable a => Plan a b -> Int
planInputSize = MS.length . planInput

planOutputSize :: Storable b => Plan a b -> Int
planOutputSize = MS.length . planOutput

execute :: (Storable a, Storable b, U.Unbox a, U.Unbox b)
            => Plan a b -> U.Vector a -> U.Vector b
execute Plan{..} v
    | n /= V.length v
        = error $ "execute: size mismatch; expected " ++ show n
                    ++ ", got " ++ show (V.length v)
    | otherwise = unsafePerformIO $ do
                        forM_ [0..n-1] $ \k -> MS.unsafeWrite planInput k
                                                $ V.unsafeIndex v k
                        planExecute
                        v' <- UM.unsafeNew m
                        forM_ [0..m-1] $ \k -> MS.unsafeRead planOutput k
                                                >>= UM.unsafeWrite v' k
                        U.unsafeFreeze v'
  where
    n = MS.length planInput
    m = MS.length planOutput

-----------------------
-- Planners: methods of plan creation.

data Planner a b = Planner {
                        inputSize :: Int -> Int,
                        outputSize :: Int -> Int,
                        creationSizeFromInput :: Int -> Int,
                        makePlan :: CInt -> Ptr a -> Ptr b -> CFlags -> IO (Ptr CPlan),
                        normalization :: Int -> Plan a b -> Plan a b
                    }


planOfType :: (Storable a, Storable b) => PlanType
                                -> Planner a b -> Int -> Plan a b
planOfType ptype Planner{..} n
  | m_in <= 0 || m_out <= 0 = error "Can't (yet) plan for empty arrays!"
  | otherwise  = unsafePerformIO $ do
    planInput <- MS.unsafeNew m_in
    planOutput <- MS.unsafeNew m_out
    MS.unsafeWith planInput $ \inP -> MS.unsafeWith planOutput $ \outP -> do
    pPlan <- makePlan (toEnum n) inP outP $ planInitFlags ptype DestroyInput
    cPlan <- newPlan pPlan
    -- Use unsafeWith here to ensure that the Storable MVectors' ForeignPtrs
    -- aren't released too soon:
    let planExecute = MS.unsafeWith planInput $ \_ ->
                        MS.unsafeWith planOutput $ \_ ->
                          withPlan cPlan fftw_execute
    return $ normalization n $ Plan {..}
  where
    m_in = inputSize n
    m_out = outputSize n

plan :: (Storable a, Storable b) => Planner a b -> Int -> Plan a b
plan = planOfType Estimate

run :: (Storable a, Storable b, U.Unbox a, U.Unbox b)
            => Planner a b -> U.Vector a -> U.Vector b
run p v = execute
            (planOfType Estimate p $ creationSizeFromInput p $ V.length v)
            v

---------------------------
-- For scaling input/output:

class Scalable a where
    scaleByD :: Double -> a -> a
    {-# INLINE scaleByD #-}

instance Scalable Double where
    scaleByD = (*)
    {-# INLINE scaleByD #-}

instance Scalable (Complex Double) where
    scaleByD s (x:+y) = s*x :+ s*y
    {-# INLINE scaleByD #-}


{-# INLINE modifyInput #-}
modifyInput :: (MS.MVector RealWorld a -> IO ()) -> Plan a b -> Plan a b
modifyInput f p@Plan{..} = p {planExecute = f planInput >> planExecute}

{-# INLINE modifyOutput #-}
modifyOutput :: (MS.MVector RealWorld b -> IO ()) -> Plan a b -> Plan a b
modifyOutput f p@Plan{..} = p {planExecute = planExecute >> f planOutput}

{-# INLINE constMultOutput #-}
constMultOutput :: (Storable b, Scalable b) => Double -> Plan a b -> Plan a b
constMultOutput !s = modifyOutput (multC s)

{-# INLINE multC #-}
multC :: (Storable a, Scalable a) => Double -> MS.MVector RealWorld a -> IO ()
multC !s v = forM_ [0..n-1] $ \k -> unsafeModify v k (scaleByD s)
  where !n = MS.length v

-- | Helper function; seems like it should be in the vector package...
{-# INLINE unsafeModify #-}
unsafeModify :: (Storable a)
                => MS.MVector RealWorld a -> Int -> (a -> a) -> IO ()
unsafeModify v k f = MS.unsafeRead v k >>= MS.unsafeWrite v k . f
