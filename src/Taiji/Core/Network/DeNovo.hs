-- | Infer network de novo from data

{-# LANGUAGE DataKinds         #-}
{-# LANGUAGE FlexibleContexts  #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE RecordWildCards   #-}

module Taiji.Core.Network.DeNovo
    ( createLinkage
    , mkNetwork
    , readNodesAndEdges
    ) where

import Control.Arrow ((&&&))
import           Bio.Utils.Misc                    (readDouble)
import           Bio.Data.Experiment
import           Bio.Pipeline.Utils                (getPath, asDir)
import Control.Monad.State.Strict
import           Bio.Data.Bed
import           Conduit
import qualified Data.IntervalMap.Strict as IM
import           Control.Lens                      hiding (pre, to)
import qualified Data.Set as S
import           Control.Monad.Reader              (asks)
import qualified Data.ByteString.Char8             as B
import           Data.CaseInsensitive              (mk)
import Data.Ord (comparing)
import           Data.List                         (foldl', maximumBy)
import qualified Data.HashMap.Strict                   as M
import           Data.Maybe                        (fromJust, mapMaybe)
import qualified Data.Text                         as T
import IGraph
import           Scientific.Workflow               hiding (_data)
import System.IO

import Taiji.Core.Network.Utils
import Taiji.Core.RegulatoryElement (findTargets)
import           Taiji.Core.Config                 ()
import           Taiji.Types
import           Taiji.Constants (edge_weight_cutoff)

createLinkage :: ( ATACSeq S ( File tag1 'Bed         -- ^ Active promoters
                             , File tag2 'Bed         -- ^ TFBS
                             )
                 , Either (File t1 'NarrowPeak) (File t2 'BroadPeak)  -- ^ promoter activity
                 , Either (File t1 'NarrowPeak) (File t2 'BroadPeak)  -- ^ enhancer activity
                 , Maybe (File '[ChromosomeLoop] 'Bed)  -- ^ HiC loops
                 , Maybe (File '[] 'Tsv)          -- ^ Expression
                 )
              -> WorkflowConfig TaijiConfig
                    (ATACSeq S (File '[] 'Other, File '[] 'Other))
createLinkage (atac, pro, enh, hic, expr) = do
    dir <- asks
        ((<> "/Network/" <> asDir (T.unpack grp)) . _taiji_output_dir)
        >>= getPath
    let netEdges = dir ++ "/edges_combined.csv"
        netNodes = dir ++ "/nodes.csv"
        bindingEdges = dir ++ "/edges_binding.csv"
    liftIO $ do
        expr' <- case expr of
            Nothing -> return M.empty
            Just e -> readExpression 1 (B.pack $ T.unpack grp ) $ e^.location
        activityPro <- fmap (bedToTree max . map (\x -> (x, fromJust $ x^.npPvalue))) $
            readBed' fl_pro
        activityEnh <- fmap (bedToTree max . map (\x -> (x, fromJust $ x^.npPvalue))) $
            readBed' fl_enh
        withFile netEdges WriteMode $ \h1 -> withFile bindingEdges WriteMode $ \h2 -> do
            B.hPutStrLn h1 ":START_ID,:END_ID,weight,:TYPE"
            B.hPutStrLn h2 $ ":START_ID,:END_ID,chr,start:int,end:int," <>
                "annotation,affinity,:TYPE"
            let conduit = findTargets active_pro tfbs hic .|
                    createLinkage_ activityPro activityEnh expr' .|
                    mapM_C (liftIO . outputEdge h1 h2)
            s <- execStateT (runConduit conduit) S.empty
            let nodeHeader = "geneName:ID,expression,expressionZScore"
            B.writeFile netNodes $ B.unlines $ (nodeHeader:) $
                map nodeToLine $ S.toList s
    return $ atac & replicates.mapped.files .~
        ( emptyFile & location .~ netNodes
        , emptyFile & location .~ netEdges )
  where
    outputEdge h1 h2 e = B.hPutStrLn hdl $ edgeToLine e
      where
        hdl = case _edge_type e of
            Combined _ -> h1
            Binding{..} -> h2
    fl_pro = either (^.location) (^.location) pro
    fl_enh = either (^.location) (^.location) enh
    (active_pro, tfbs) = atac^.replicates._2.files
    grp = atac^.groupName._Just
{-# INLINE createLinkage #-}

createLinkage_ :: BEDTree Double   -- ^ Promoter activities
               -> BEDTree Double   -- ^ Enhancer activities
               -> M.HashMap GeneName (Double, Double)   -- ^ Gene expression
               -> ConduitT (GeneName, ([BED], [BED]))
                           NetEdge
                           (StateT (S.Set NetNode) IO) ()
createLinkage_ act_pro act_enh expr = concatMapMC $ \(geneName, (ps, es)) -> do
    let tfEnhancer = M.toList $ fmap getBestMotif $ M.fromListWith (++) $
            mapMaybe (getWeight act_enh) es
        edgeEnhancer = flip concatMap tfEnhancer $ \(tfName, sites) ->
            flip map sites $ \st -> NetEdge
                { _edge_from = tfName
                , _edge_to = geneName
                , _edge_type = Binding
                    { _edge_binding_locus = convert st
                    , _edge_binding_annotation = "enhancer"
                    , _edge_binding_affinity = fromJust $ st^.score }
                }
        tfPromoter = M.toList $ fmap getBestMotif $ M.fromListWith (++) $
            mapMaybe (getWeight act_pro) ps
        edgePromoter = flip concatMap tfPromoter $ \(tfName, sites) ->
            flip map sites $ \st -> NetEdge
                { _edge_from = tfName
                , _edge_to = geneName
                , _edge_type = Binding
                    { _edge_binding_locus = convert st
                    , _edge_binding_annotation = "promoter"
                    , _edge_binding_affinity = fromJust $ st^.score }
                }
        tfs = M.toList $ fmap (lp 2 . map (fromJust . (^.score))) $
            M.fromListWith (++) $ tfEnhancer ++ tfPromoter
        (geneExpr, scaledGeneExpr) = M.lookupDefault (0.1, 0) geneName expr
        geneNode = NetNode { _node_name = geneName
                           , _node_expression = Just geneExpr
                           , _node_scaled_expression = Just scaledGeneExpr }
    modify' $ S.insert geneNode
    edgeCombined <- forM tfs $ \(tfName, w) -> do
        let (tfExpr, scaledTfExpr) = M.lookupDefault (0.1, 0) tfName expr
            tfNode = NetNode { _node_name = tfName
                             , _node_expression = Just tfExpr
                             , _node_scaled_expression = Just scaledTfExpr }
        modify' $ S.insert tfNode
        return $ NetEdge { _edge_from = tfName
                         , _edge_to = geneName
                         , _edge_type = Combined (w * sqrt tfExpr) }
    return $ edgePromoter ++ edgeEnhancer ++ edgeCombined
  where
    getBestMotif xs = runIdentity $ runConduit $
        mergeBedWith (maximumBy (comparing (^.score))) xs .| sinkList
    getWeight act bed = case IM.elems (intersecting act bed) of
        [] -> Nothing
        xs -> let p = transform_peak_height $ maximum xs
                  w = sqrt $ transform_site_pvalue (fromJust $ bed^.score) * p
              in if w < edge_weight_cutoff
                then Nothing else Just (getTFName bed, [bed & score .~ Just w])
    getTFName x = mk $ head $ B.split '+' $ x^.name._Just
    transform_peak_height x = 1 / (1 + exp (-(x - 5)))
    transform_site_pvalue x' = 1 / (1 + exp (-(x - 5)))
      where
        x = negate $ logBase 10 $ max 1e-20 x'
{-# INLINE createLinkage_ #-}

lp :: Int -> [Double] -> Double
lp p = (**(1/fromIntegral p)) . foldl' (+) 0 . map (**fromIntegral p)
{-# INLINE lp #-}

repressors :: S.Set GeneName
repressors = S.fromList ["AES", "BAZ1A", "Bhlhe40", "BRCA2", "BRIP1", "BTG2"
    , "CALR", "CARHSP1", "CBX2", "CSDA", "DEK", "DENND4A", "FBXW7", "FHL2"
    , "HELLS", "HMGN2", "HOPX", "KLF11", "Klf8", "MCM5", "MLF1IP", "MLLT1"
    , "Prdm1", "RBL2", "SAP30", "SART1", "SETD8", "SLA2", "SMAD7", "SNAI3"
    , "SUV39H1", "TSC22D3", "TSHZ1", "UHRF1", "ZFP598", "Atf7", "Bcl6", "Bcl6b"
    , "Bhlhe40", "Cic", "Crebzf", "Egr3", "Erf", "Etv3", "Hbp1", "Hesx1"
    , "Hic1", "Hivep1", "Homez", "Id1", "Jdp2", "Kdm2b", "Klf8", "Mbd1"
    , "Mbd2", "Mecp2", "Msc", "Mypop", "Nfil3", "Nr2f6", "Pou6f1", "Prdm1"
    , "Sp2", "Srebf2", "Tgif1", "Tgif2", "Vax2", "Zbtb1", "Zbtb4", "Zfp128"
    , "Zfp161", "Zfp187", "Zfp202", "Zfp263", "Zfp281", "Zscan10"]

both' :: S.Set GeneName
both' = S.fromList ["Ar", "Arid5a", "Arid5b", "Ascl2", "Atf2", "Atf3", "Bach1"
    , "Bach2", "Cdc5l", "Cebpa", "Cebpb", "Cebpd", "Cebpg", "Clock"
    , "Creb1", "Crem", "Ctcf", "Cux1", "E2f1", "E2f2", "E2f3", "E2f4"
    , "E2f5", "E2f6", "E2f7", "E2f8", "Egr1", "Elf2", "Elf3", "Elk3"
    , "Elk4", "Eomes", "Esr1", "Esrra", "Etv6", "Fli1", "Foxk1", "Foxo1"
    , "Foxo3", "Foxo4", "Foxp1", "Foxp3", "Gata1", "Gata3", "Glis2", "Hes1"
    , "Hes5", "Hhex", "Hmbox1", "Hmga2", "Hoxa7", "Hoxb4", "Hsf1", "Hsf4"
    , "Ikzf1", "Irf1", "Irf2", "Irf3", "Jun", "Junb", "Klf12", "Klf4", "Lcor"
    , "Lef1", "Mafk", "Max", "Mef2a", "Mef2c", "Mitf", "Mlx", "Mnt", "Myb"
    , "Myc", "Mzf1", "Nanog", "Nfatc1", "Nfatc2", "Nfatc3", "Nfix", "Nfkb1"
    , "Nr1h2", "Nr2c1", "Nr3c1", "Nr4a2", "Nr4a3", "Ovol1", "Patz1", "Pax5"
    , "Phf21a", "Plag1", "Pou2f1", "Pou5f1", "Ppard", "Pparg", "Prdm11"
    , "Rara", "Rarb", "Rarg", "Rel", "Rela", "Rfx3", "Rfx5", "Rora", "Rorc"
    , "Rreb1", "Runx1", "Runx2", "Rxra", "Sfpi1", "Smad3", "Smarcc2"
    , "Snai1", "Sox5", "Sox6", "Sp1", "Srebf1", "Stat3", "Stat6", "Tbx21"
    , "Tcf3", "Tcf4", "Tcf7", "Tcf7l1", "Tcf7l2", "Terf2", "Trp53", "Vdr"
    , "Xbp1", "Yy1", "Zbtb7b", "Zeb1", "Zfp143", "Zfp238", "Zfp691", "ATAD2"
    , "ATXN1", "BRCA1", "DNMT3B", "E2f2", "E2f3", "Fli1", "FMR1", "FOXM1"
    , "Gata3", "HLX", "HMGB2", "HMGN5", "HNRNPD", "Hsf4", "Irf1", "KDM4A"
    , "Klf12", "MAF", "PHF10", "PHF19", "Runx2", "Runx3", "Rxra", "Tbx21"
    , "TCF19", "Tcf7l2", "TLE1", "YY2", "ZBED4"]

-- | Build the network from files containing the information of nodes and edges.
mkNetwork :: FilePath  -- ^ nodes
          -> FilePath  -- ^ edges
          -> IO (Graph 'D NetNode Double)
mkNetwork nodeFl edgeFl = do
    nodeMap <- M.fromList . map ((_node_name &&& id) . nodeFromLine) .
        tail . B.lines <$> B.readFile nodeFl
    fmap (filterGraph (repressors <> both')) $ runResourceT $ fromLabeledEdges' edgeFl (toEdge nodeMap)
  where
    toEdge nodeMap fl = sourceFileBS fl .| linesUnboundedAsciiC .|
        (dropC 1 >> mapC f) 
      where
        f l = ( ( M.lookupDefault undefined (mk f2) nodeMap
                , M.lookupDefault undefined (mk f1) nodeMap )
              , readDouble f3 )
          where
            [f1,f2,f3,_] = B.split ',' l
{-# INLINE mkNetwork #-}

filterGraph :: S.Set GeneName -> Graph 'D NetNode Double -> Graph 'D NetNode Double
filterGraph s gr = efilter f gr
  where
    f ((fr, to), _) = not $ _node_name (nodeLab gr to) `S.member` s &&
        fromJust (_node_scaled_expression $ nodeLab gr fr) > 0

-- | Read network files as nodes and edges
readNodesAndEdges :: FilePath   -- ^ nodes
                  -> FilePath   -- ^ edges
                  -> IO ([NetNode], [((GeneName, GeneName), Double)])
readNodesAndEdges nodeFl edgeFl = do
    nds <- map nodeFromLine . tail . B.lines <$> B.readFile nodeFl
    es <- map f . tail . B.lines <$> B.readFile edgeFl
    return (nds, es)
  where
    f l = ( ( mk f2, mk f1), readDouble f3 )
        where
        [f1,f2,f3,_] = B.split ',' l