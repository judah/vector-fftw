{- |
This module provides normalized versions of the transforms in @fftw@.

All of the transforms are normalized so that

 - Each transform is orthogonal, i.e., preserves the sum-of-squares norm of its input.

 - Each backwards transform is the inverse of the corresponding forwards transform.

(Both conditions only hold approximately, due to floating point precision.)

-}
module Math.FFT.Vector.Orthogonal(
                -- * Creating and executing 'Plan's
                run,
                plan,
                execute,
                -- * Complex-to-complex transforms
                dft,
                idft,
                -- * Real-to-complex transforms
                dftR2C,
                dftC2R,
                -- * Discrete cosine transforms
                -- $dct_inv
                dct2,
                dct3,
                dct4,
                ) where

import Math.FFT.Vector.Base
import qualified Math.FFT.Vector.Unnormalized as U
import Data.Complex
import qualified Data.Vector.Storable.Mutable as MS
import Control.Monad.Primitive(RealWorld)

dft, idft :: Planner (Complex Double) (Complex Double)
dft = U.dft {normalization = \n -> constMultOutput $ 1 / sqrt (toEnum n)}
idft = U.idft {normalization = \n -> constMultOutput $ 1 / sqrt (toEnum n)}

dftR2C :: Planner Double (Complex Double)
dftR2C = U.dftR2C {normalization = \n -> modifyOutput $
                    complexR2CScaling (sqrt 2) n
        }

dftC2R :: Planner (Complex Double) Double
dftC2R = U.dftC2R {normalization = \n -> modifyInput $
                    complexR2CScaling (sqrt 0.5) n
        }

complexR2CScaling :: Double -> Int -> MS.MVector RealWorld (Complex Double) -> IO ()
complexR2CScaling !t !n !a = do
    let !s1 = sqrt (1/toEnum n)
    let !s2 = t * s1
    let len = MS.length a
    -- Justification for the use of unsafeModify:
    -- The output size is 2n+1; so if n>0 then the output size is >=1;
    -- and if n even then the output size is >=3.
    unsafeModify a 0 $ scaleByD s1
    if odd n
        then multC s2 (MS.unsafeSlice 1 (len-1) a)
        else do
            unsafeModify a (len-1) $ scaleByD s1
            multC s2 (MS.unsafeSlice 1 (len-2) a)

{- $dct_inv
Note that 'dct4' is its own inverse, whereas 'dct2' and 'dct3' are inverses of each other.
-}

dct4 :: Planner Double Double
dct4 = U.dct4 {normalization = \n -> constMultOutput $ 1 / sqrt (2 * toEnum n)}

dct2 :: Planner Double Double
dct2 = U.dct2 {normalization = \n -> modifyOutput $ \a -> do
    let n' = toEnum n
    let !s1 = sqrt $ 1 / (4*n')
    let !s2 = sqrt $ 1 / (2*n')
    unsafeModify a 0 (*s1)
    multC s2 (MS.unsafeSlice 1 (MS.length a-1) a)
    }

dct3 :: Planner Double Double
dct3 = U.dct3 {normalization = \n -> modifyInput $ \a -> do
    let n' = toEnum n
    let !s1 = sqrt $ 1 / n'
    let !s2 = sqrt $ 1 / (2*n')
    unsafeModify a 0 (*s1)
    multC s2 (MS.unsafeSlice 1 (MS.length a-1) a)
    }
