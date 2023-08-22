##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: characterize the weired clusters in the RA samples and decide if discard them or not
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed May 10 16:18:44 2023
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

levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA")
data_version = "_d2_d2.5_d3_d3.5_d4_d5"

names(cols) = levels
cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/cluster_filtering_RAsamples/')
system(paste0('mkdir -p ', outDir))

##########################################
# import data 
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', data_version, '_',
                           species, version.analysis, '.rds'))

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

# quickly run clustering
aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
p1 + p2

ggsave(filename = paste0(outDir, '/UMAP_RA_condition_clusters.pdf'), 
       width = 14, height = 6)


##########################################
# find markers for cluster 7 and 8
##########################################
Idents(aa) = aa$clusters
cluster7.markers <- FindMarkers(aa, only.pos = TRUE, ident.1 = 7, ident.2 = c(0,1,3), 
                                min.pct = 0.25, logfc.threshold = 0.5)
head(cluster7.markers, n = 8)

FeaturePlot(aa, features = rownames(cluster7.markers)[1:9])
ggsave(filename = paste0(outDir, '/FeaturePlots_markers_cluster7.pdf'), 
       width = 10, height = 8)


cluster8.markers <- FindMarkers(aa, only.pos = TRUE, ident.1 = 8, ident.2 = c(0,1,3), 
                                min.pct = 0.25, logfc.threshold = 0.5)
head(cluster8.markers, n = 10)

FeaturePlot(aa, features = rownames(cluster8.markers)[1:9])
ggsave(filename = paste0(outDir, '/FeaturePlots_markers_cluster8.pdf'), 
       width = 10, height = 8)

all.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.4)
all.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10
DoHeatmap(aa, features = top10$gene) + NoLegend()

ggsave(filename = paste0(outDir, '/Heatmap_more.markers_cluster7_cluster8.pdf'), 
       width = 14, height = 20)

