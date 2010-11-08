module Main (main) where

import Data.Vector.Storable as V
import Math.FFT.Vector.Unnormalized as U
import Math.FFT.Vector.Invertible as I
import Criterion.Main
import Data.Complex

main = do
    let mkSize planner k = do
            let !v = V.replicate k 17 :: Vector (Complex Double)
            let !p = plan planner k
            return $ bench ("fftPlanA" Prelude.++ show k) $ whnf (\w -> execute p w) v
    let mkSize1 k = do
            let !v = V.replicate k 17 :: Vector (Complex Double)
            let !p = plan dft k
            return $ bench ("fftPlanB" Prelude.++ show k) $ whnf (\w -> execute p w) v
    let mkSize2 k = do
            let !v = V.replicate k 17 :: Vector (Complex Double)
            return $ bench ("fftPlanC" Prelude.++ show k) $ whnf (\w -> I.run dft w) v
    
    let sizes = [64,256,2048,8192]
    benches <- Prelude.mapM (mkSize U.idft) sizes
    benches1 <- Prelude.mapM mkSize1 sizes
    benches2 <- Prelude.mapM mkSize2 sizes
    defaultMain $ benches Prelude.++ benches1 Prelude.++ benches2
