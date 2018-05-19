{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}
module Taiji.Core (builder) where

import           Bio.Data.Experiment
import           Bio.Data.Experiment.Parser   (readHiC, readHiCTSV)
import           Control.Lens
import           Control.Monad.IO.Class       (liftIO)
import           Control.Monad.Reader         (asks)
import qualified Data.Map.Strict              as M
import           Data.Maybe                   (fromJust)
import           Data.Monoid                  ((<>))
import           Scientific.Workflow

import           Taiji.Core.Network
import           Taiji.Core.RegulatoryElement
import           Taiji.Types                  (_taiji_input)

builder :: Builder ()
builder = do
    nodePS 1 "Find_Active_Promoter" 'findActivePromoters $ do
        note .= "Identify active promoters. Promoters are defined by " <>
            "-5000 ~ +1000 regions around annotated transcription start sites. " <>
            "If a promoter is overlapped with ATAC-seq peaks, we assume it is active."

    nodeS "HiC_Read_Input" [| \_ -> do
        input <- asks _taiji_input
        liftIO $ do
            hic <- if ".tsv" == reverse (take 4 $ reverse input)
                then readHiCTSV input "HiC"
                else readHiC input "HiC"
            return $ getHiCLoops hic
        |] $ do
            submitToRemote .= Just False
            note .= "Read HiC loops from input file."

    node' "Find_TF_Target_Prep" [| \(activePro, tfbs, hic) ->
        let mkDict xs = M.fromList $ map (\x ->
                (x^.groupName, head $ x^..replicates.folded.files)) xs
            lookup' x xs = M.lookup (x^.groupName) xs
            tfbs' = mkDict tfbs
            hic' = mkDict hic
        in flip map activePro $ \e ->
            ( e & replicates.mapped.files %~ (\f -> (f, fromJust $ lookup' e tfbs'))
            , M.lookup (e^.groupName) hic' )
        |] $ do
            note .= "Prepare for parallel execution."
            submitToRemote .= Just False

    nodePS 1 "Find_TF_Target" [| \(x1,x2) -> findTargets x1 x2 |] $ do
        remoteParam .= "--mem=30000 -p gpu"
        note .= "Assign TFs to their target genes. We borrow the concept of " <>
                "gene regulatory domain from GREAT. Gene regulatory domain " <>
                "definition: Active promoters are used to define the basal " <>
                "regulatory domains of genes. The gene regulatory domain is " <>
                "extended in both directions to the nearest gene's basal domain " <>
                "but no more than the maximum extension in one direction." <>
                "TF binding sites located in gene regulatory domains are then " <>
                "assigned to corresponding genes."
    path ["Find_TF_Target_Prep", "Find_TF_Target"]

    node' "Create_Linkage_Prep" [| \(assignment, atac_peaks, chip_peaks) ->
        let getFile x = (x^.groupName._Just, x^.replicates._2.files)
            atacFileMap = fmap Left $ M.fromList $ map getFile atac_peaks
            chipFileMap = M.fromList $ map getFile chip_peaks
        in flip map assignment $ \e ->
            let grp = e^.groupName._Just
                pro = M.findWithDefault undefined grp atacFileMap
                enh = M.findWithDefault pro grp chipFileMap
            in (e, pro, enh)
        |] $ do
            note .= "Prepare for parallel execution."
            submitToRemote .= Just False
    nodePS 1 "Create_Linkage" 'createLinkage $ return ()
    path ["Create_Linkage_Prep", "Create_Linkage"]

    node' "Compute_Ranks_Prep" [| \(links, expr) -> zip links $ repeat expr |] $ do
        note .= "Prepare for parallel execution."
        submitToRemote .= Just False
    nodePS 1 "Compute_Ranks" 'computeRanks $ do
        note .= "Perform personalized Pagerank."
        remoteParam .= "--mem=20000 -p gpu"
    nodeS "Output_Ranks" 'outputRanks $ do
        remoteParam .= "--mem=20000 -p gpu"
    path ["Compute_Ranks_Prep", "Compute_Ranks", "Output_Ranks"]
