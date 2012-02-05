{-# LANGUAGE BangPatterns #-}
module Main where

import Data.Vector.Storable as V
import Data.Vector.Storable.Mutable as M
import Numeric.FFT.Vector.Unnormalized as U
import Numeric.FFT.Vector.Plan
import Data.Numskell.Vector as N
import Criterion.Main
import Data.Complex

import System.Environment

main = do
    [n] <- fmap (fmap Prelude.read) getArgs
    let numIters = 1000 * 10
    vIn <- M.unsafeNew n
    vOut <- M.unsafeNew n
    N.write vIn =: pure n 17
    let !p = plan dft n
    N.sequence_ $ pure numIters $ executeM p vIn vOut
