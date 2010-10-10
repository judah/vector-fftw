module Math.FFT.Vector.Invertible(
                    -- * Creating and executing 'Plan's
                    run,
                    plan,
                    execute,
                    -- * Complex-to-complex transforms
                    U.dft,
                    idft,
                    -- * Real-to-complex transforms
                    U.dftR2C,
                    dftC2R,
                    -- * Discrete cosine transforms
                    U.dct1,
                    idct1,
                    U.dct2,
                    idct2,
                    U.dct3,
                    idct3,
                    U.dct4,
                    idct4,
                    -- * Discrete sine transforms
                    U.dst1,
                    idst1,
                    U.dst2,
                    idst2,
                    U.dst3,
                    idst3,
                    U.dst4,
                    idst4,
                    ) where

import Math.FFT.Vector.Base
import qualified Math.FFT.Vector.Unnormalized as U
import Data.Complex

idft :: Planner (Complex Double) (Complex Double)
idft = U.idft {normalization = \n -> constMultOutput $ 1 / toEnum n}

dftC2R :: Planner (Complex Double) Double
dftC2R = U.dftC2R {normalization = \n -> constMultOutput $ 1 / toEnum n}

-- OK, the inverse of each unnormalized operation.

idct1 :: Planner Double Double
idct1 = U.dct1 {normalization = \n -> constMultOutput $ 1 / toEnum (2 * (n-1))}

idct2 :: Planner Double Double
idct2 = U.dct3 {normalization = \n -> constMultOutput $ 1 / toEnum (2 * n)}

idct3 :: Planner Double Double
idct3 = U.dct2 {normalization = \n -> constMultOutput $ 1 / toEnum (2 * n)}

idct4 :: Planner Double Double
idct4 = U.dct4 {normalization = \n -> constMultOutput $ 1 / toEnum (2 * n)}

idst1 :: Planner Double Double
idst1 = U.dst1 {normalization = \n -> constMultOutput $ 1 / toEnum (2 * (n+1))}

idst2 :: Planner Double Double
idst2 = U.dst3 {normalization = \n -> constMultOutput $ 1 / toEnum (2 * n)}

idst3 :: Planner Double Double
idst3 = U.dst2 {normalization = \n -> constMultOutput $ 1 / toEnum (2 * n)}

idst4 :: Planner Double Double
idst4 = U.dst4 {normalization = \n -> constMultOutput $ 1 / toEnum (2 * n)}

