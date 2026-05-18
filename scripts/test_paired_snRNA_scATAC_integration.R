##########################################################################
##########################################################################
# Project: RA project 
# Script purpose: test integration methods for paired scRNA and scATAC
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon May 18 15:36:11 2026
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

library(RColorBrewer)
require(Matrix)
library(RaceID)
require(tictoc)
require(gridExtra)
require(cowplot)

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

outDir = paste0(resDir, '/mNTs_multiome_RA_scRNAscATAC_integraton_test')
if(!dir.exists(outDir)) dir.create(outDir)


##########################################
# features  
##########################################
sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

chromfactors = readRDS(file = '../data/annotations/human_chromatin_remodelers_Epifactors.database.rds')
chromfactors = as.character(unlist(sapply(chromfactors, firstup)))

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

sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))


########################################################
########################################################
# Section I : check gene examples of Pluripotent, RA response, AP identity genes
# test subclustering for each time points
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir,
                           'seuratObj_RAmultiome_snRNA.normalized.umap_scATAC.normalized',
                           'noMatureNeurons.smallClusters_wnn.umap.clusters.rds'))


DefaultAssay(aa) <- "RNA"
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000) # find subset-specific HVGs
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 100, min.dist = 0.1, 
              reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)


DefaultAssay(aa) <- "ATAC"
#aa <- FindTopFeatures(aa, min.cutoff = 'q5')
aa <- FindTopFeatures(aa, min.cutoff = 'q10')
aa <- RunTFIDF(aa, method = 1)
aa <- RunSVD(aa)

p1 <- ElbowPlot(aa, ndims = 30, reduction="lsi")
p2 <- DepthCor(aa, n = 30)
p1 | p2

aa <- RunUMAP(object = aa, reduction = 'lsi', 
              dims = 2:30, 
              n.neighbors = 30, 
              min.dist = 0.1, 
              reduction.name = "umap_atac",
              reduction.key = "UMAPATAC_"
)

DimPlot(object = aa, label = TRUE, reduction = 'umap_atac', group.by = 'condition') 

aa <- FindMultiModalNeighbors(aa, reduction.list = list("pca", "lsi"), 
                              dims.list = list(1:30, 2:30))

aa <- RunUMAP(aa, nn.name = "weighted.nn", 
              reduction.name = "wnn.umap", 
              reduction.key = "wnnUMAP_", min.dist = 0.1)

p1 <- DimPlot(aa, reduction = "umap_rna", group.by = "condition", 
              label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(aa, reduction = "umap_atac", group.by = "condition", label = TRUE, 
              label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(aa, reduction = "wnn.umap", group.by = "condition", label = TRUE, 
              label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3  & NoLegend() & theme(plot.title = element_text(hjust = 0.5))



