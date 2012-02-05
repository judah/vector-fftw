{-# LANGUAGE BangPatterns #-}
module Main (main) where

import Data.Vector.Storable as V
import Numeric.FFT.Vector.Unnormalized as U
import Numeric.FFT.Vector.Invertible as I
import Numeric.FFT.Vector.Unitary as O
import Criterion.Main
import Data.Complex

main = do
    let mkSize planner k = do
            let !v = V.replicate k 17 :: Vector (Complex Double)
            let !p = plan planner k
            return $ bench (show k) $ whnf (\w -> execute p w) v
    let mkSize1 k = do
            let !v = V.replicate k 17 :: Vector (Complex Double)
            let !p = plan U.dft k
            return $ bench (show k) $ whnf (\w -> execute p w) v
    let mkSize2 k = do
            let !v = V.replicate k 17 :: Vector (Complex Double)
            let !p = plan O.dft k
            return $ bench (show k) $ whnf (\w -> execute p w) v
    
    let mkSize3 k = do
            let !v = V.replicate k 17 :: Vector (Complex Double)
            return $ bench (show k) $ whnf (\w -> I.run O.dft w) v
    
    let sizes = [64,256,2048,8192]
    let sizedGroup name f = fmap (bgroup name) $ Prelude.mapM f sizes
    benches <- sizedGroup "U.idft_1" (mkSize U.idft)
    benches1 <- sizedGroup "U.dft" mkSize1
    benches2 <- sizedGroup "O.dft" mkSize2
    benches3 <- sizedGroup "run_O.dft" mkSize3
    defaultMain $ [benches, benches1, benches2, benches3]
