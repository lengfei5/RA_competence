##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: search for the early bias/heterogeneity before RA, RA and/or RA withdraw
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jan 27 11:53:12 2026
##########################################################################
##########################################################################

rm(list = ls())

version.analysis = '_R16597_mNT_10xmultiome_reseq_20240517'

resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/'

source(paste0(functionDir,  'functions_scATAC.R'))
source(paste0(functionDir, 'functions_scRNAseq.R'))
source(paste0(functionDir, 'functions_Visium.R'))

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
require(SeuratObject)
library(data.table)
library(DropletUtils)
library(edgeR)
library(future)
library(tictoc)

options(future.globals.maxSize = 200 * 1024^3)
set.seed(1234)
mem_used()

species = 'mNT_multiome'

### Parameters for this project
levels = c("day0_beforeRA", "day1_beforeRA", 
           "day2_beforeRA",
           "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", 'day3.5_RA',
           "day4_RA", "day5_RA", "day6_RA",
           "day2.5_noRA", "day3_noRA", 'day3.5_noRA', "day4_noRA", "day5_noRA", "day6_noRA")

# manually set colors by Hannah
library(RColorBrewer)
library("viridis")
cols = rep(NA, length = 16)
names(cols) = levels
cols[grep('_beforeRA', names(cols))] = colorRampPalette((brewer.pal(n = 3, name ="Greys")))(3)
#cols[1:3] = viridis(3)
cols[grep('_noRA', names(cols))] = colorRampPalette((brewer.pal(n = 6, name ="Blues")))(6)
cols[grep('_RA', names(cols))] = colorRampPalette((brewer.pal(n = 7, name ="OrRd")))(7)


## subset the mutliome conditions
levels_sels = c("day2_beforeRA",
                "day2.5_RA", "day3_RA", "day3.5_RA",  "day4_RA", "day5_RA", 
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA")

cc = names(cols)
cc[which(cc == 'day3_RA.rep1')] = 'day3_RA'
names(cols) = cc

cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/mNTs_multiome_RA_',
                'searching_earlyBias')
if(!dir.exists(outDir)) dir.create(outDir)


##########################################
# features  
##########################################
sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

features_1 = unique(c('Sox2', 'Sox1', 'Tubb3', 'Elavl3', 
                      'Irx3', 'Irx5', 'Pax3', 'Pax7',
                      'Pax6', 'Olig2', 'Nkx2-9', 'Nkx2-2', 
                      'Nkx6-1', 'Foxa2', 'Arx', 'Shh'
)) # DV overview

features_2 = c('Dhrs3', 'Rarg', 'Cyp26a1',
               'Pou3f1', 'Hoxa1', 'Gas1', 'Spry4', 'Sox11', 
               'Cdh1', 'Cdh2', 'Shh', 'Rfx4', 'Zfp42', 'Tcf15', 'Prrx2', 'Gdf3',
               'Etv5', 'Fgf4', 'Otx2', 'Zscan10', 'Apoe', 'Peg10', 'Klf9', 'Tshz1', 'Skil', 'Zfp703')

features = unique(c(c('Pax6', 'Foxa2', 'Sox1', 'Sox2', 'Tubb3', 'Shh', 'Arx',
                      'Zfp703', 'Lef1', 'Irx5', 'Pou5f1', 'Otx2', 'Adgra2', 'Hoxb4', 
                      'Nkx2-2', 'Nkx2-9', 'Nkx6-1', 'Olig2', 'Pax3', 'Pax7', 'Cyp26a1', 'Dhrs3'), 
                    features_1,
                    features_2)) # marker wanted by Hannah

gene_examples = unique(c('Foxa2', 'Pax6', c('Zfp42', 'Tcf15', 'Skil', 'Lef1',
                                            'Sox2', 'Pou5f1', 'Sall4', 'Tdgf1', # pluripotency markers
                                            'Nanog', 'Nr5a2', #'Prdm14', 
                                            'Klf4', 'Fgf4', 'Esrrb', 'Tcf3', 'Tbx3'), # naive pluripotency
                         c('Zfp42', 'Tcf15', 'Skil',
                           'Fgf5', 'Otx2', 'Pou3f1', 'Lef1', 'Dnmt3b', 'Dnmt3a',	
                           'Foxd3', 'Utf1', 'Tcf15', 'Zic3', 'Rhox5', 'Etv5', 'Etv4',	
                           'Lin28b', 'Sox4', 'Sox3', 'Sox11'
                         ),
                         c('Lhx1','Eomes', 'Sox2', 'Hoxb4', 'Hoxb5', 'Hoxb6','Zfp703'),
                         c('Zfp42', 'Tcf15', 'Skil', 'Lef1', 'Dhrs3', 'Rarg', 'Cyp26a1'),
                         features,
                         features_1,
                         features_2
                         
))


########################################################
########################################################
# Section I : check gene examples of Pluripotent, RA response, AP identity genes
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir,
                           'seuratObj_RAmultiome_snRNA.normalized.umap_scATAC.normalized',
                           'noMatureNeurons.smallClusters_wnn.umap.clusters.rds'))

sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
# aa <- RunUMAP(aa, nn.name = "weighted.nn", 
#               reduction.name = "wnn.umap", 
#               reduction.key = "wnnUMAP_", min.dist = 0.1)

FeaturePlot(aa, features = c('Pax6', 'Foxa2'), reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_Pax6_Foxa2.pdf'), 
       height = 6, width = 17)

## # naive pluripotency
genes = c('Nanog', 'Nr5a2', 'Prdm14', 'Klf4', 'Fgf4', 'Esrrb', 'Tcf3', 'Tbx3')

FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_naivPluripotent.pdf'), 
       height = 12, width = 18)


## pluripotent markers
genes = c('Zfp42', 'Tcf15', 'Skil', 'Lef1','Sox2', 'Pou5f1', 'Sall4', 'Tdgf1')

FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_pluripotentMarkers.pdf'), 
       height = 12, width = 18)

## RA response
genes = c('Dhrs3', 'Rarg', 'Cyp26a1', 'Cyp26b1', 'Zfp703', 'Otx2', 'Stra8', 'Tead1', 'Tead2', 'Tead3', 'Tead4',
          'Yap1', 'Tbx1', 'Tbx2', 'T', 'Runx1', 'Rbp4', 'Rxra', 'Rxrb', 'Myb', 'Lefty1', 'Krt13', 
          'Igf2r', 'Gsk3b', 'Gata4', 'Gata6', 'Creb1', 'Asxl1', 'Ascl1', 'Bmp6')

FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_RAresponse.pdf'), 
       height = 36, width = 27)

## Hox genes
genes = rownames(aa)[grep('Hox', rownames(aa))]
genes = c('Hoxa1', 'Hoxa3', 'Hoxb4')

FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_HoxGenes.pdf'), 
       height = 12, width = 18)

## WNT genes
genes = c('Sfrp2', 'Zeb2', 'Zfp703', 'Sox17', 'Snai2', 'Sall1', 'Rspo2', 'Rspo4', 'Mdk', 'Lrrk1', 
          'Lmx1a', 'Lgr4', 'Lgr5', 'Dkk1', 'Dkk4', 'Tcf15', 'Wnt9a', 'Wnt4', 'Lef1')

FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_Wnt.pdf'), 
       height = 21, width = 27)

FeaturePlot(aa, features = c('Lef1', 'Tcf15', 'Zeb2', 'Zfp703', 'Axin2', 'Tcf7', 'Tcf7l1', 'Tcf7l2', 
                             'Dvl1', 'Dvl2'), 
            reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_Wnt_2.pdf'), 
       height = 16, width = 27)

genes = c('Gdf3',
          'Etv5', 'Fgf4', 'Otx2', 'Zscan10', 'Apoe', 'Peg10', 'Klf9', 'Tshz1')
FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_HannahsExamples.pdf'), 
       height = 16, width = 27)

## FGF genes
genes = c('Fgf4', 'Frs2', 'Akt1', 'Spry1', 'Spry2',
          'Etv5', 'Etv4', 'Dusp6', 'Egr1')
FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_FGFgenes.pdf'), 
       height = 16, width = 27)

## BMP and TGFb genes
genes = c('Bmp2', 'Bmp4', 'Bmp7', 'Gdf5', 'Gdf6', 'Id1', 'Id2', 'Id3', 'Msx1', 'Msx2', 'Fst', 'Smad1', 
          'Smad5', 'Smad4', 'Smad2', 'Smad3', 'Lefty1', 'Lefty2', 'Nog', 'Grem1')

p1 = FeaturePlot(aa, features = genes, reduction = 'wnn.umap')
plot(p1)
ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_BMP.TGFbeta.pdf'), 
       height = 21, width = 27)




########################################################
########################################################
# Section II: VarID method to identify the high expression noise
# 
########################################################
########################################################


