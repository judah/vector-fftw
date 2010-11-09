{-# LANGUAGE ScopedTypeVariables #-}
-- Run: 
-- ./Bench -u all.csv
-- runghc PlotResults.hs all.csv
module Main where

import Text.CSV
import System.Process
import qualified Data.Map as Map
import System.IO
import System.FilePath
import Data.List
import System.Directory

import System.Environment
import Control.Exception


main = do
    [file] <- getArgs
    pd <- readRunInfos file
    renderPlotData pd

-- each entry in this map is one line in the graph.
type PlotData = Map.Map String [RunInfo]

data RunInfo = RunInfo {
                runSize :: Int,
                runMean :: Double 
            }
    deriving Show

runGnuplot :: [String] -> IO ()
runGnuplot ss = do
    (Just hIn,Nothing,Nothing,ph) <- createProcess (proc "gnuplot" ["-persist"])
                                        {std_in = CreatePipe}
    hPutStrLn hIn $ unlines $ ["set terminal x11"] ++ ss ++ ["quit"]
    waitForProcess ph
    return ()

-- read from a criterion -u dump.       
readRunInfos :: FilePath -> IO PlotData
readRunInfos path = do
    -- ignore the header
    Right (_:csv) <- parseCSVFromFile path
    return $ foldl' (\m (n,r) -> Map.insertWith (++) n [r] m) Map.empty
            $ fmap parseRun
            $ filter (any (not . null)) csv
  where
    parseRun csv = let
        sizeStr = takeFileName $ csv !! 0
        testName = takeFileName $ takeDirectory $ csv !! 0
        size = tryRead "size" sizeStr
        mean = tryRead "mean" $ csv !! 1
        in (testName, RunInfo size mean)

tempDataFile :: [RunInfo] -> IO FilePath
tempDataFile runs = do
    dir <- getTemporaryDirectory
    (path,h) <- openTempFile dir "plots.dat"
    mapM_ (hPutStrLn h)
        $ map (\r -> show (runSize r) ++ " " ++ show (runMean r))
            runs
    hClose h
    return path
                    

renderPlotData :: PlotData -> IO ()
renderPlotData pd = do
    let pdList = Map.toList pd
    fs <- mapM tempDataFile $ fmap snd pdList
    let cmds = [ "set logscale x"
               , "set logscale y"
               , "set xtics"
               , "plot " ++ intercalate ", " 
                            [show f 
                                    ++ " using 1:2:xticlabels(1)"
                                    ++ " title " ++ show n
                                    ++ " with lines "
                                | (n,f) <- zip (fmap fst pdList) fs]

               ]
    runGnuplot cmds
    
    mapM_ removeFile fs

tryRead :: Read a => String -> String -> a
tryRead descr s = mapException (\(e::SomeException) -> userError ("descr: " ++ show s ++ " " ++ show e))
                    $ read s
