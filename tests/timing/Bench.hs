module Main (main) where

import Data.Vector.Unboxed as V
import Math.FFT.Vector.Unnormalized as U
import Math.FFT.Vector.Invertible as I
import Criterion.Main
import Data.Complex

main = do
    let mkSize planner k = do
            let !v = V.replicate k 17 :: Vector (Complex Double)
            let !p = plan planner k
            return $ bench ("fftPlan" Prelude.++ show k) $ whnf (\w -> execute p w) v
    benches <- Prelude.mapM (mkSize U.idft) [64,128,256,512,1024]
    benches' <- Prelude.mapM (mkSize I.idft) [64,128,256,512,1024]
    defaultMain $ benches Prelude.++ benches'
