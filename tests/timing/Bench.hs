module Main (main) where

import Data.Vector.Unboxed as V
import Math.FFT.Vector.Unnormalized as FFT
import Criterion.Main
import Data.Complex

main = do
    let n = 1024
    let v = V.replicate n 17 :: Vector (Complex Double)
    let !p = plan dft n
    defaultMain
        [ bench "fftRun" $ whnf (FFT.run dft) v
        , bench "fftPlan" $ whnf (execute p) v
        ]
