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

group_1_genes :: S.Set GeneName
group_1_genes = S.fromList ["Arid3b","CENPK","CHD7","ELL2","FKBP2","FKBP5","H2AFZ","HAT1","HMGB3","JMJD5","LMO4","MED30","Mybl2","NAP1L1","Nfic","Nr2c2","PAXIP1","RFC3","SNRPD1","SUB1","TAF4A","TOE1","Trp73","UPF1","ATAD2","ATXN1","BRCA1","DNMT3B","E2f2","E2f3","Fli1","FMR1","FOXM1","Gata3","HLX","HMGB2","HMGN5","HNRNPD","Hsf4","Irf1","KDM4A","Klf12","MAF","PHF10","PHF19","Runx2","Runx3","Rxra","Tbx21","TCF19","Tcf7l2","TLE1","YY2","ZBED4","AES","BAZ1A","Bhlhe40","BRCA2","BRIP1","BTG2","CALR","CARHSP1","CBX2","CSDA","DEK","DENND4A","FBXW7","FHL2","HELLS","HMGN2","HOPX","KLF11","Klf8","MCM5","MLF1IP","MLLT1","Prdm1","RBL2","SAP30","SART1","SETD8","SLA2","SMAD7","SNAI3","SUV39H1","TSC22D3","TSHZ1","UHRF1","ZFP598"]

group_4_genes :: S.Set GeneName
group_4_genes = S.fromList ["Pcyox1","Tspan13","Prmt6","Fbxo31","Prickle1","Morn1","Adamts3","Tsen2","Slc52a2","Peg13","Stx11","2010002M12Rik","Afp","Prrt2","2310001H17Rik","Galk1","Ninj1","Gprin3","Rab31","Zfp839","App","Gm166","Chchd6","Pear1","1700052K11Rik","C130039O16Rik","Eif4enif1","Aars2","Rasgef1a","Dstyk","Ipmk","Jund","1110019D14Rik","Ipo4","Grhl1","2900005J15Rik","Inpp4b","Dusp6","Gm10825","Nf2","Ift88","Eno2","Adssl1","Inpp5f","Dusp4","Ears2","Dlg4","Wdr41","1700049G17Rik","Xcl1","Enpp2","Nudt8","Mfsd4","Ces2d-ps","Trim62","Trip6","Hck","Alpl","2010016I18Rik","Pecr","Scmh1","Ift74","Dnase1","1190005F20Rik","Pfkm","2410022L05Rik","Apex1","Cnn3","Inppl1","Spint1","Rltpr","Oas1g","Agl","Gm17455","Tnfrsf4","A230050P20Rik","Gna12","Slc30a3","Fbf1","Scamp1","Thap2","Atp1a3","Adck4","Ralb","Ptger2","Irf6","Fsd2","Acadsb","A430107P09Rik","Zfp937","Parp12","Trat1","Morn4","Tgif1","Dhx34","Akap8l","H6pd","Hdac4","Gdap10","Nr1d2","Acsf2","Pdk1","Lclat1","Dgki","Serinc5","Ryk","Insr","Taf1c","Plcb4","Grcc10","Mtmr7","Bcl2l14","Adamts6","Akr7a5","Chst3","Cd101","Bambi-ps1","Lgmn","2900008C10Rik","Polr1b","Gabrr2","Fto","Rundc3a","Per2","Fam179b","Vangl2","Egr2","Fasn","Hyal2","Ptpn13","Ccdc30","BC017612","Tmem86a","Tmem170b","Pbx4","Nle1","Cxxc5","Mapre3","Ankrd37","Pkp2","Dtd1","Rnf128","Atad3a","Maml3","Acyp1","Dgkz","Espn","Lpcat1","B230217O12Rik","Wfs1","Rdh12","Prss2","Prss41","Gnl3","Mx1","Hmgn3","Nfix","Csgalnact1","Arid5b","Padi4","Rasgrp4","A430071A18Rik","Ntrk3","Rhob","D030028A08Rik","Il6st","Oas2","Atp6v1g3","Wdr26","Zmynd8","Zfp870","Icam5","Pstpip2","Zfp53","Slc22a15","Gabrd","1500011K16Rik","BC057022","Fbxl20","Dctd","Irf7","Lysmd2","Zfp119b","Npnt","Frat1","Mbnl2","4930452B06Rik","1700021C14Rik","Dtx4","Ccdc164","Fam83g","Bbs12","Cxx1b","Dusp14","Pex11c","Gpm6b","Fam132a","Pctp","Zhx2","Ext1","Fam120c","Slc5a3","Nab1","Pdxk","Ccdc122","Ifrd2","Pprc1","Naprt1","Atp9a","Spsb1","Zwint","Vamp5","C80913","Trit1","Eml5","Top1mt","Atcay","Ephb6","5330417C22Rik","Bbs1","Tppp","Btla","Bcor","Tank","Fbxo10","Opn3","Tbc1d19","Tox","2010011I20Rik","Cbx7","Pde8b","Zfp239","Rtp4","Etnk1","Fam71b","Mpp1","Lyrm4","Dennd2d","Wdr34","Fam108c","Zfp526","Tmem158","Arrb1","Mov10","Slc26a2","Galnt11","Cyp4v3","Fam110a","3110082I17Rik","Camkk1","Pdcd1","Rab11fip1","Pitpnm1","B3gnt9-ps","Zfp128","Pou2f2","Msh6","Rhbdd3","Eme2","Zfp169","Galr3","4933403G14Rik","2310035K24Rik","Kifap3","Ttll12","Limk1","Tmem106c","Gramd4","Zfp619","Mdn1","Cxx1c","2410002F23Rik","Dync2li1","Pikfyve","Khk","Zfp821","Akap6","Ankrd55","Ctxn1","Aicda","Idua","Tasp1","Orai1","Spock2","Rnf19a","Pip4k2a","Wdr59","Ifi44","Gadd45b","Herpud2","Dnajb2","Tctex1d2","Ctss","Tiam1","2810405K02Rik","Ift80","Tspan3","9430091E24Rik","Zc4h2","Hivep2","Tmem131","Maz","N4bp2","2010109I03Rik","4930578N16Rik","Zfp763","Echdc3","Serac1","Cul7","Pcdhgb7","Oasl1","Efnb1","Clec11a","4632415K11Rik","E130308A19Rik","4933403F05Rik","Ttc27","Gpr133","Dst","Dock6","Tmem67","Gdpd1","Nbeal1","Stat5a","Oas1a","Myo10","Gm16894","Srxn1","Sh3d19","Icam4","Dnase2a","Acsl3","Msi2","1700019D03Rik","Batf","Isg20","Ttll5","Kcnab1","Arxes2","Padi3","Fscn1","Rpgr","Trim59","Acacb","Dleu2","Gbp9","Gm19705","Lef1","Nmnat3","Fyb","Tuft1","Mrrf","Ncor2","Plcxd1","Adck3","Pptc7","Thnsl2","Lipt2","Pcdhgb4","St3gal2","Shpk","Rab32","Wfikkn2","Repin1","Adk","Hexa","Tcp11l2","Als2cl","Iqcb1","4632427E13Rik","Ralgps1","Tmc3","Tcf7","A530088E08Rik","F2rl1","Coro2b","Cyp39a1","S1pr2","Spaca1","Cd200","Slc25a42","Ramp1","Supt3h","Rnf167","Usp6nl","Gas8","Fhit","Lnx2","Ryr3","Gstt3","Tdrd7","Il16","Klhdc2","Trmt61a","Zfp566","Hhat","Mblac1","Ap1m2","Sgsh","Crtc3","Nadsyn1","Gcg","Stard5","Ece1","Siah3","Jmy","Hap1","Tubb2a","Cyp4f13","Slc39a11","Klf3","Rassf6","Fars2","Slc25a26","Cebpd","Abhd8","Ppcdc","Tanc1","Myo1h","Rnaset2b","Hsdl1","Nrp","Capsl","Neu1","Dnase1l3","Fuz","Cass4","Slc19a1","Sv2a","Thbs3","A430088P11Rik","Fntb","Elmo3","Tjp2","Mrpl52","Btg1","A930024E05Rik","C920006O11Rik","Stxbp1","Ndrg3","Mpzl1","Grwd1","Plk1s1","Rel","Tert","Nenf","Rab26","Tmc4","Thumpd2","Nckipsd","Endog","Sgsm2","Stx2","Chd3","Dag1","Tpd52","Nefh","Wdr19","Neil1","Cela1","Sh3rf1","Ldhb","Ttc30b","Tlr7","Spata7","Epcam","Foxp1","Btbd6","Ing2","Hs1bp3","Rwdd3","Wls","Zfp667","Chia","6430548M08Rik","Tgfbr3","4933406I18Rik","1110051M20Rik","Kirrel3","Bckdhb","Arsg","Ccdc28b","Gm9895","Phpt1","Zscan18","Ssx2ip","Ankrd13d","Cry2","Gm14718","Pou2af1","Plxnb1","Fam70a","Phactr2","Ttc8","Pgm2","Zfp518b","2610015P09Rik","Tmem42","Casd1","Ttn","Slc35f5","Lipe","Zan","Slc29a2","Jrk","B4galnt4","Rgs11","Nrbp2","Rnf144a","Slc27a1","Pkp4","Hspb11","Pcdhgc5","Pgcp","Ppp4r1l-ps","Gspt2","Plekha8","Nme4","Metap1d","Rpp38","Acaa1b","Kdm5b","Aldh5a1","Slc41a3","Endod1","Dnaja4","Ccdc91","Cd83","Naip5","BC035044","Myc","Naip2","Filip1l","Setbp1","Ptprs","Sertad2","Cbx6","Art2b","Padi1","Ppil6","Acsl1","Zfhx2","Ccdc67","Acer2","Plekho1","Zfp422","Vmac","Sec14l1","Snn","4930599N23Rik","Lrig1","Acpp","Tmie","Gemin8","Mapk8ip1","Zfp551","Creb3l2","Abtb2","Irak1bp1","Amigo2","Angptl2","Spats1","Mmachc","Ndrg4","Ramp3","Abhd4","Chd6","Mbd3","Apol9b","Mettl8","Ankrd46","Ccdc86","9230110C19Rik","Rdh9","Phf15","Wnt3","Tm7sf2","Isg15","Xkrx","Pla2g12a","Als2cr4","Sqrdl","Cln6","Tmem220","Tns1","Rere","Nphp1","Kif5c","Ryr1","Prpf40b","Ddx24","Tsr2","B9d1","Rapsn","Rnf41","Rph3al","Cd81","Ehbp1","Pus7","Rnf24","Meis3","Prnp","Crtc1","Bicd1","Slc14a1","Gm10548","Nlk","Epb4.1l5","Prdm5","Rplp1","Tagap1","Socs5","Itpr2","Cd320","Otud1","Reep2","Azi1","Zhx3","Bend3","1110032F04Rik","Fastkd1","Pik3r2","Pydc3","Fhdc1","Marcksl1","Gm88","Myo1e","Pank4","Vamp1","Cdk5r1","Bmp7","Zfp507","Akr1b8","AI854703","Aldoc","P2rx4","Tgfb3","Dgkg","Pydc4","2200002D01Rik","Snx30","Asb2","Zfr2","Ces2c","Slc25a1","Fgfr1op","Cryz","Fam43a","E130311K13Rik","Tmem38a","Hpcal1","1700001G11Rik","Socs3","Fbxw8","Osgin1","Ctns","Ldlrap1","Rdh13","Zmat1","Ifi27l2a","1500011H22Rik","Spnb1","Prmt3","Rad51l1","Mmab","Srm","Ociad2","Rps6kc1","Zc3h12d","Ms4a4c","Tbc1d5","Dnajc12","Dgka","Zfp1","Ttll9","Ftl1","Il1r2","Atmin","4932438H23Rik","Pomt1","Npepl1","Colq","Cul9","Hsh2d","Mx2","Lyl1","Ptrh1","Grhpr","Kdm6b","Akap7","2010107G23Rik","Msrb2","D730040F13Rik","Blcap","Polg","Kremen1","Rpgrip1","Lats2","Optn","Rpa1","2310022B05Rik","Slc26a11","Abcg1","Zfp446","Ehhadh","Ccdc126","Fam46c","Usp54","Msl3l2","Stx4a","Slc36a4","Prmt2","Zfp78","Bcl9","Clpb","Pfkfb4","Dock4","Rab11fip3","Gas2","Numb","Armc9","Pter","Taf4b","Actn1","Eea1","Lass6","Sp6","H2-Q2","Slc12a8","Cbr1","Acad12","Cd22","Slc46a1","Dnalc4","Fam101b","Erf","Pde2a","Klhl3","Senp8","Sostdc1","Wdr35","Cxcr4","Zfp318","D18Ertd653e","Sh3tc1","Tnfrsf8","Txlnb","Ncf2","Parp16","Slc43a1","Trib2","Serp1","Dnajc6","Usp11","Klhl35","Cand2","Tmem17","Slc15a1","Oasl2","Kat2a","Gaa","Zfp608","Tbc1d4","Tor3a","Iqcc","Tmem218","Thop1","Mansc1","Zbtb20","Gm5918","H2-Q1","Arhgap39","Spred2","Ly6d","1190007I07Rik","Zmym6","Ssh3","Amz2","Mccc1","Trio","Mrps6","Tnfrsf23","Fmnl3","Car12","Renbp","Slc1a4","Zfp523","Zfp770","2610301B20Rik","4933427D14Rik","Gns","Scn11a","2610528E23Rik","Zfp939","Gab2","Slfn5","Casc1","Wdpcp","Limd2","Arhgap29","Xlr3b","Itga9","Rassf2","Itgb3bp","Mttp","Rhod","Gcnt1","Plxna3","Ift122","Gprasp1","Lrp12","Ss18l1","Prkcz","Otud3","Decr2","Tcea2","Tpbg","Cstad","Ccdc157","Rab39b","Gpatch4","Rcl1","Cd74","Dusp23","Ly96","Penk","H2-Ke6","A530072M11Rik","Trim2","Arl5c","Dtx3","Ifit1","H2-Ob","4930417O13Rik","Zfp167","Cd82","Klhl26","Asah2","Cyth3","Eef2k","Mef2b","Ldoc1l","Sdc3","Map3k5","Mtm1","Ccny","Tcn2","Dbp","Nr4a1","2810454H06Rik","Zfp777","Il7","Rap1gap","Lrrc45","BC021614","Kif13a","3000002C10Rik","Fam131a","Pank1","Mrgpre","Aldh6a1","Rnf2","Adar","6430527G18Rik","Tesc","Tom1l2","Dap3","Art3","Daxx","Asns","Ctsl","Scrn2","Sat1","Dsel","Mtmr11","Gucy2e","Ppp1r15a","Cd160","D630037F22Rik","Gjb2","BC031353","Ppp1r14b","Chka","Ublcp1","Afap1","Krt222","Smox","Pex11a","Gpr146","Asap2","Wisp1","Pcsk4","Ralgds","Jarid2","Phtf1","Mgrn1","Zfp831","Ube2z","Tmprss11e","Mboat7","Dos","E130304F04Rik","4831440E17Rik","Arhgef10","Acot11","Stx16","Ccdc28a","Fmo5","Tmem30b","2310047B19Rik","Bphl","Tmem106a","Rhoh","Cited2","BC016423","P2rx7","Ccdc142","Ephx4","Mamdc4","Morn2","Mif","Dusp16","Ssh2","Trerf1","Fam26f","Bnip3","BC048403","Ly6k","1600016N20Rik","Fbxo32","Arhgef11","Amacr","Trpm6","D130040H23Rik","Cd3eap","Vipr1","Nucb2","Syt11","Wdr60","Sema7a","Decr1","Ptk2","Acad10","Taz","Bcat1","Ahi1","Rapgef4","Ecm1","Pacsin1","Ms4a6c","Sh3bp4","Cacna2d4","Tmem149","Ispd","Acy3","4930481A15Rik","Trim30d","6430562O15Rik","Sdc1","Fos","BC006779","Por","Pfn2","Rpp40","Dnajb14","Ppargc1b","Trub1","Stat2","3110043O21Rik","Urgcp","Nrip1","Tfpi","1110032A03Rik","Zdhhc23","Dkk3","Abcc3","Abcc5","Cttn","Ampd3","Rasa2","Cacna1b","Hspbap1","Mlh3","Ak4","Smyd5","Fam134b","Csad","Tpi1","Tshz2","Rabl5","Tbc1d30","Fbln1","Fcrl1","Trim8","Itga7","Pbx2","Clip3","Tada2b","1300018J18Rik","Rdh5","Csrnp2","Cables1","Akap2","Cxcr5","Nt5e","Sfxn5","Sccpdh","1700019L03Rik","Matk","Noc3l","Nelf","AW011738","Odc1","BC017647","Nanos1","Tnfrsf13b","Eomes","Csrnp1","Klc3","Sdr42e1","Tspan32","Ica1l","Nav2","Xpc","Naglu","5730409E04Rik","E330033B04Rik","Rasgrp1","Bcl7a","Usp33","Gtf2i","Trpm1","Kif3a","Il21","Lepre1","Il1r1","Icos","Kif1b","Tfdp2","Gm15708","Nfkbie","Evi5l","Slc45a4","Chchd4","Stx6","Cerk","Mettl1","Utrn","Ddx60","Zfp882","Gpx1","Znrf1","Magee2","Zfp235","Ttc3","Abi2","Adamtsl4","Flcn","Psen2","Aldh1b1","Pde6d","Arhgap32","Sall2","Zfand3","Larp1b","1700019E19Rik","Bbs2","Mtap2","Ubtd1","Gstt2","Pcbp3","Il5ra","Hspa1a","2610008E11Rik","Klrb1c","2210013O21Rik","Sox4","Cd27","Slc4a7","Rab6b","Tmem191c","Ms4a6d","Cdk20","Ibtk","Neu3","Tubg2","C1qtnf4","2310015A10Rik","Gm15645","Ankrd28","Zfp36l1","Grtp1","Thada","Akap12","Sh2d1a","Cd69","Cacna1d","Acp5","Xrcc6bp1","Ccdc102a","Pfkl","Tnfrsf22","1810043G02Rik","Rsad2","I830012O16Rik","Pold4","Zdhhc1","Pde11a","Col1a1","Shmt1","Galnt9","Folr4","Rhbdd2","Egln3","Rabl2","Iqsec1","Rgs10","4732471D19Rik","Rab24","Mapk14","Tnfaip8","Arl6","Ydjc","Ntf5","C78339","Immp2l","Slmo1","Hs6st1","Tlr2","Lipa","Zfp365","Zfp94","2010015L04Rik","2310016M24Rik","Tox2","Herc6","Josd2","Ptp4a2","Isyna1","Fam81a","Ccs","Dhdh","Cdc42ep4","Stxbp4","6330403K07Rik","Dynlt3","March8","Ccl6","Inadl","Luzp1","Tprg","Gpr19","1190002H23Rik","Ptcra","Zfp395","Egr1","Malat1","Sgpp2","Nop16","Fam109a","Efha2","Serpinf1","Otub2","Ikzf2","Oaz2","Klhdc3","Spry1","Tnfrsf14","Cnga1","Farp2","Rangrf","1700027F09Rik","Gm129","Mr1","Enox","Zfp467","Synj2","Cnr2","Zfp346","Cnnm3","Pdlim7","Inpp5a","Pard6g","Zfyve28","Elovl6","Ltbp3","Slc25a13","Gyk","Oxsr1","Gdpd5","Slc41a1","Mccc2","Bace2","Tha1","Rps6ka2","Cd200r3","Ifit3","Il6ra","Ncrna00086","Jak2","Tspan5","Fam116b","Gpr137","Zbtb4","Txnrd3","Acta2","Ift81","AU041133","Baz2b","Ephx1","Tmem38b","Fam118a","Prdm16","Neurl3","Gna13","Asap1","A130040M12Rik","Nhedc2","Pink1","Pgpep1","Ttc16","Dhx58","Rpap3","Cdh23","Bpgm","Zfp14","Dpp7","Ttc28","Siah2","Apba3","Wdr92","Trappc9","Prkar1b","Tmem9","Fam73b","Lrrc23","Fam98c","Tmem8b","Gnaq","Sipa1l1","Hyi","Stk3","P4ha1","Mxd4","Hook2","Siae","Gng3","Mest","Mex3b","Slc7a7","1110034G24Rik","Fggy","Pdzk1ip1","D830046C22Rik","Snx33","Pkd1","Pkn3","Kdm4b","Atp1b1","Olfr1033","Nit2","Tsga14","Fas","L2hgdh","Slc25a19","4930473A06Rik","C030006K11Rik","Brdt","Pck2","Aacs","Ikzf4","Clybl","Nlrc4","Pcdhgb6","Gimap3","Lrrc56","Crebl2","Bcl3","Nfkbia","Cd68","Slc39a13","Frat2","Snhg12","2610035D17Rik","Ggt7","Ppfibp1","Fgfbp3","Plekhg2","Tlcd1","Zfp28","Cask","Ccdc141","Mgst2","Sdhaf1","Plagl1","BC018242","Ptch1","Eif2c1","Smyd4","Acvr1b","Dok3","Ppp2r3a","Pnpla7","Sbno2","Gm9766","Tfcp2","Rhobtb2","1700010I14Rik","Ccm2","Vps13c","Angel1","Zfp783","Zfp36","Slc25a22","Cebpa","Ypel3","Dtwd1","Bag3"]

-- | Build the network from files containing the information of nodes and edges.
mkNetwork :: FilePath  -- ^ nodes
          -> FilePath  -- ^ edges
          -> IO (Graph 'D NetNode Double)
mkNetwork nodeFl edgeFl = do
    nodeMap <- M.fromList . map ((_node_name &&& id) . nodeFromLine) .
        tail . B.lines <$> B.readFile nodeFl
    
    gr <- runResourceT $ fromLabeledEdges' edgeFl (toEdge nodeMap)
    return $ changeWeights $ filterGraph (repressors <> both') gr
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

keepGroup1TFonly :: Graph 'D NetNode Double -> Graph 'D NetNode Double
keepGroup1TFonly gr = nfilter f gr
  where
    f (i, nd) = _node_name nd `S.member` group_1_genes || null (pre gr i)

changeWeights :: Graph 'D NetNode Double -> Graph 'D NetNode Double
changeWeights gr = nmap f gr
  where
    f (i, nd) | _node_name nd `S.member` group_4_genes = nd{_node_scaled_expression = Just 10}
              | otherwise = nd{_node_scaled_expression = Just 1}

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