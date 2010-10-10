-- | Unnormalized versions of the transforms.  This is faster but may cause math problems.
module Math.FFT.Vector.Unnormalized(
                    -- * Creating and executing 'Plan's
                    run,
                    plan,
                    execute,
                    -- * Complex-to-complex transforms
                    dft,
                    idft,
                    -- * Real-to-complex transforms
                    dftC2R,
                    dftR2C,
                    -- * Discrete cosine transforms
                    dct1,
                    dct2,
                    dct3,
                    dct4,
                    -- * Discrete sine transforms
                    dst1,
                    dst2,
                    dst3,
                    dst4,
                    ) where

import Math.FFT.Vector.Base
import Foreign
import Foreign.C
import Data.Complex

#include <fftw3.h>

-- | Whether the complex fft is forwards or backwards.
type CDirection = CInt

-- | The type of the cosine or sine transform.
type CKind = (#type fftw_r2r_kind)

foreign import ccall unsafe fftw_plan_dft_1d
    :: CInt -> Ptr (Complex Double) -> Ptr (Complex Double) -> CDirection
        -> CFlags -> IO (Ptr CPlan)

foreign import ccall unsafe fftw_plan_dft_r2c_1d
    :: CInt -> Ptr Double -> Ptr (Complex Double) -> CFlags -> IO (Ptr CPlan)

foreign import ccall unsafe fftw_plan_dft_c2r_1d
    :: CInt -> Ptr (Complex Double) -> Ptr Double -> CFlags -> IO (Ptr CPlan)

foreign import ccall unsafe fftw_plan_r2r_1d
    :: CInt -> Ptr Double -> Ptr Double -> CKind -> CFlags -> IO (Ptr CPlan)

dft1D :: CDirection -> Planner (Complex Double) (Complex Double)
dft1D d = Planner {
            inputSize = id,
            outputSize = id,
            creationSizeFromInput = id,
            makePlan = \n a b -> fftw_plan_dft_1d n a b d,
            normalization = const id
            }

dft, idft :: Planner (Complex Double) (Complex Double)
dft = dft1D (#const FFTW_FORWARD)
idft = dft1D (#const FFTW_BACKWARD)

dftR2C :: Planner Double (Complex Double)
dftR2C = Planner {
            inputSize = id,
            outputSize = \n -> n `div` 2 + 1,
            creationSizeFromInput = id,
            makePlan = fftw_plan_dft_r2c_1d,
            normalization = const id
        }

dftC2R :: Planner (Complex Double) Double
dftC2R = Planner {
            inputSize = \n -> n `div` 2 + 1,
            outputSize = id,
            creationSizeFromInput = \n -> 2 * (n-1),
            makePlan = fftw_plan_dft_c2r_1d,
            normalization = const id
        }

r2rPlanner :: CKind -> Planner Double Double
r2rPlanner kind = Planner {
                    inputSize = id,
                    outputSize = id,
                    creationSizeFromInput = id,
                    makePlan = \n a b -> fftw_plan_r2r_1d n a b kind,
                    normalization = const id
                }

dct1, dct2, dct3, dct4 :: Planner Double Double
dct1 = r2rPlanner (#const  FFTW_REDFT00)
dct2 = r2rPlanner (#const  FFTW_REDFT10)
dct3 = r2rPlanner (#const  FFTW_REDFT01)
dct4 = r2rPlanner (#const  FFTW_REDFT11)

dst1, dst2, dst3, dst4 :: Planner Double Double
dst1 = r2rPlanner (#const  FFTW_RODFT00)
dst2 = r2rPlanner (#const  FFTW_RODFT10)
dst3 = r2rPlanner (#const  FFTW_RODFT01)
dst4 = r2rPlanner (#const FFTW_RODFT11)
