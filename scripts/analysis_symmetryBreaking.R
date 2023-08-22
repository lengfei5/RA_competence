##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: search for genes with asymmetric expression at early time points
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Mar 30 16:04:20 2023
##########################################################################
##########################################################################
library(Seurat)
library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)


# levels_sels = c("day3_RA.rep1", "day3.5_RA", "day4_RA")
# data_version = "_d3_d3.5_d4"

#levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1")
#data_version = "_d2_d2.5_d3"

#levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA")
#data_version = "_d2_d2.5_d3_d3.5_d4"

#levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA")
#data_version = "_d2_d2.5_d3_d3.5_d4"

# levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1",  "day3_RA.rep2",
#                 "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")

levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA")
data_version = "_d2_d2.5_d3_d3.5_d4_d5"

names(cols) = levels
cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/RA_symetryBreaking/sparse_featureSelection', data_version)
system(paste0('mkdir -p ', outDir))

##########################################
# tfs and sps annotations 
##########################################
sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

# all 16 samples
# aa = readRDS(file = paste0(RdataDir, 
#                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
#                            'cellCycleScoring_annot.v1_savedUMAP.v1_', species, version.analysis, '.rds'))

# only RA samples incl. dya2_beforeRA
aa = readRDS(file = paste0(RdataDir,
                           'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_savedUMAP.subs.v2_', species, version.analysis, '.rds'))

Idents(aa) = aa$condition
DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)

FeaturePlot(aa, features = c('Zfp42', 'Tcf15', 'Skil', 'Lef1'))

ggsave(filename = paste0(outDir, '/asymmetric_feature_expression.pdf'), 
       width = 14, height = 10)



##########################################
# 
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', data_version, '_',
                           species, version.analysis, '.rds'))

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

FeaturePlot(aa, features = c('Zfp42', 'Tcf15', 'Skil', 'Lef1'))

ggsave(filename = paste0(outDir, '/asymmetric_feature_expression.pdf'), 
       width = 14, height = 10)


aa <- RunUMAP(aa, 
              dims = NULL, features = c('Pax6', 'Foxa2', 'Sox2', 'Pou5f1', 'Sox1'),
              metric = "euclidean",
              a = 10, 
              b = 2.5,
              n.neighbors = 100, min.dist = 0.5, 
              reduction.key = 'UMAP.feature', 
              reduction.name = "umap_feature")

DimPlot(aa, label = TRUE, repel = TRUE, reduction = 'umap_feature',
        group.by = 'condition', raster=FALSE)




