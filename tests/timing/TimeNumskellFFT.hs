module Main(main) where

import qualified Data.Vector.Unboxed as V
import Data.Complex
import Math.FFT.Vector.Unnormalized

import Control.Exception
import Control.Monad(forM_)

main = do
    let n = 1024
    let p = plan dft n
    V.forM_ (V.enumFromN 0 10000) $ \k -> do
        let v = V.replicate n (k:+0)
        evaluate $ execute p v
        -- evaluate $ run dft v
