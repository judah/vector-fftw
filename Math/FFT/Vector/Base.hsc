-- | A basic interface between Vectors and the fftw library.
module Math.FFT.Vector.Base where

import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Generic as V
import Control.Monad.Primitive
import Data.Complex as DC
import Foreign
import Foreign.C
import Data.Bits


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

data Plan sh a b = Plan {
                    planInput :: VS.MVector RealWorld a,
                    planOutput :: VS.MVector RealWorld b,
                    planExecute :: IO ()
                }

{-
-- TODO: Allow arbitrary Shape sh
execute :: (Storable a, Elt a, Storable b, Elt b)
            => Plan (Z:.Int) a b -> Array (Z:.Int) a -> Array (Z:.Int) b
execute p a
    | extent a /= planInputSize p = error "execute: shape mismatch: expected "
                                        ++ show (planInputSize p)
                                        ++ ", got " ++ show (extent a)
    | otherwise = unsafePerformIO
                        $ withForeignPtr (planInputArray p) $ \p_in ->
                        $ withForeignPtr (planOutputArray p) $ \p_out -> do
                            forM_ [0..n-1] $ \k -> pokeElemOff p_in k $ a :! (Z:.k)
                            planExecute p
                            -- TODO: Hacky McHackerstein!
                            return $! force $ fromFunction (Z:.n)
                                                (\_:.k -> unsafeInlineIO
                                                            $ peekElemOff p_out k)
  where
    _:.n = extent a

-----------------------
-- Planners: methods of plan creation.

data Planner sh a b = Planner {
                        plannerSizes :: sh -> (sh,sh) -- (input,output)
                        creationSizeFromInput :: sh -> sh,
                        makePlan :: sh -> Ptr a -> Ptr b -> CFlags -> IO (Ptr CPlan),
                        normalizeInput :: sh -> Plan sh a b -> Plan sh a b
                    }


planOfType :: (Storable a, Storable b) => PlanType
                                -> Planner Int a b -> Int -> Plan Int a b
planOfType ptype Planner{..} n
  | inputSize n <= 0 || outputSize n <= 0 = error "Can't (yet) plan for empty arrays!"
  | otherwise  = unsafePerformIO $ do
    planInput@(MArray _ inFP) <- newFFTWArr_ $ inputSize n
    planOutput@(MArray _ outFP) <- newFFTWArr_ $ outputSize n
    withForeignPtr inFP $ \inP -> withForeignPtr outFP $ \outP -> do
    pPlan <- makePlan (toEnum n) inP outP $ flagsInt ptype DestroyInput
    cPlan <- newPlan pPlan
    let planExecute = withPlan cPlan fftw_execute
    return $ normalization n $ Plan {..}

plan :: (Storable a, Storable b) => Planner a b -> Int -> Plan a b
plan = planOfType Estimate
-}
