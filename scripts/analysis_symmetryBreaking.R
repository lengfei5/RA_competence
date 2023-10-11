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
library(future)

options(future.globals.maxSize = 120000 * 1024^2)

source('functions_utility.R')
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

outDir = paste0(resDir, '/RA_symetryBreaking/subclustering_earlyTimePoints_', data_version)
system(paste0('mkdir -p ', outDir))


##########################################
# import the all data for RA treatment
##########################################
aa = readRDS(file = paste0(RdataDir,
                           'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_savedUMAP.subs.v2_', species, version.analysis, '.rds'))

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

ggsave(filename = paste0(outDir, '/UMAP_conditions.pdf'), 
       width = 10, height = 6)

FeaturePlot(aa, features = c('Zfp42', 'Tcf15', 'Skil', 'Lef1'))
ggsave(filename = paste0(outDir, '/asymmetric_feature_expression.pdf'), 
       width = 14, height = 10)

# aa <- RunUMAP(aa, 
#               dims = NULL, features = c('Pax6', 'Foxa2', 'Sox2', 'Pou5f1', 'Sox1'),
#               metric = "euclidean",
#               a = 10, 
#               b = 2.5,
#               n.neighbors = 100, min.dist = 0.5, 
#               reduction.key = 'UMAP.feature', 
#               reduction.name = "umap_feature")
# 
# DimPlot(aa, label = TRUE, repel = TRUE, reduction = 'umap_feature',
#         group.by = 'condition', raster=FALSE)


########################################################
########################################################
# Section : subclustering each time points
# 
########################################################
########################################################
p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
#p2 = DimPlot(aa, group.by = 'clusters', label = TRUE, repel = TRUE)

conditions = unique(aa$condition)
print(conditions)
bb = aa;
rm(aa)

conditions = c('day2_beforeRA', 'day2.5_RA', 'day3_RA.rep1')
noisyGenes = readRDS(file = paste0(RdataDir, 'topGenes_localVaribility.gene.expression_VarID2.rds'))
gene_examples = unique(c(gene_examples, noisyGenes))

for(cc in conditions)
{
  #cc = c('day2_beforeRA')
  # cc = c('day2.5_RA')
  # cc = c('day3_RA.rep1')
  # cc = c('day3.5_RA')
  # cc = c('day4_RA')
  # cc = c('day5_RA')
  #cc = c('day3_RA.rep1', 'day3.5_RA')
  cc = c('day2_beforeRA', 'day2.5_RA', 'day3_RA.rep1')
  
  outDir_cc = paste0(outDir, '/', paste0(cc, collapse = "_"), '/')
  system(paste0('mkdir -p ', outDir_cc))
  
  aa = subset(bb, cells = colnames(bb)[which(!is.na(match(bb$condition, cc)))])
  
  Idents(aa) = aa$condition
  table(aa$condition)
 
  DimPlot(aa, group.by = 'condition', label = TRUE)
  FeaturePlot(aa, features = c('Pax6', 'Foxa2'))
  
  ##########################################
  # UMAP and clustering without sparse feature selection
  ##########################################
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000) # find subset-specific HVGs
  variableGenes = VariableFeatures(object = aa)
  select.method = 'HVGs2000_pca.weighted_'
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  #variableGenes = setdiff(VariableFeatures(object = aa), 'Lypd2')
  # variableGenes = intersect(variableGenes, c(tfs, sps))
  
  #variableGenes = VariableFeatures(object = aa)
  cat(length(variableGenes), ' variable genes using \n ')
  
  aa <- RunPCA(aa, features = variableGenes, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) + NoLegend()
  p11 = FeaturePlot(aa, features = 'nCount_RNA')
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(outDir_cc, '/UMAP_RA_condition_clustering_cellcyclePhase_', select.method, '.pdf'), 
         width = 12, height = 6)
  
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(outDir_cc, '/featureExamples_originalUMAP.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = aa, 
                           features = gene_examples[!is.na(match(gene_examples, rownames(aa)))])
  
  dev.off()
  
  ##########################################
  # regress out the cell cycle  
  ##########################################
  xx = aa
  xx$CC.Difference <- xx$S.Score - xx$G2M.Score
  #xx <- ScaleData(xx, vars.to.regress = c('nCount_RNA', "CC.Difference"), features = variableGenes)
  
  xx <- ScaleData(xx, vars.to.regress = c('nCount_RNA', "S.Score", "G2M.Score"), 
                  features = variableGenes)
  
  # cell cycle effects strongly mitigated in PCA
  xx <- RunPCA(xx, features = variableGenes, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(xx, ndims = 50)
  
  xx <- RunUMAP(xx, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
  DimPlot(xx, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  
  xx <- FindNeighbors(xx, dims = 1:10)
  xx <- FindClusters(xx, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  p1 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) + NoLegend()
  p11 = FeaturePlot(xx, features = 'nCount_RNA')
  p2 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(outDir_cc, '/UMAP_RA_condition_clustering_cellcyclePhaseRegression_', 
                           select.method, '.pdf'), 
         width = 12, height = 6)
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(outDir_cc, '/featureExamples_UMAP_cellcycleRegression.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = xx, 
                           features = gene_examples[!is.na(match(gene_examples, rownames(aa)))])
  
  dev.off()
  
  rm(xx)
  
  ##########################################
  # Removing cell cycle-related genes in the HVGs
  # original code :
  # https://bioconductor.org/books/3.12/OSCA/cell-cycle-assignment.html#removing-cell-cycle-effects
  ##########################################
  library(scater)
  #require(batchelor)
  
  Idents(aa) = aa$condition
  
  # Identifying the likely cell cycle genes between phases,
  # using an arbitrary threshold of 5%.
  scaledMatrix = GetAssayData(aa, slot = c("scale.data"))
  
  diff <- getVarianceExplained(scaledMatrix, data.frame(phase = aa$Phase))
  diff = data.frame(diff, gene = rownames(diff))
  diff = diff[order(-diff$phase), ]
  
  hist(diff$phase, breaks = 100); abline(v = c(1:5), col = 'red')
  
  genes_discard = diff$gene[which(diff$phase>5)]
  cat(length(genes_discard), 'genes to discard \n')
  
  hvgs = VariableFeatures(FindVariableFeatures(subset(aa, features = setdiff(rownames(aa), genes_discard)), 
                                               selection.method = "vst", nfeatures = 2000))
  # 
  # discard <- diff > diff_cutoff
  # summary(discard)
  # cat('cut off -- ', diff_cutoff, '--', sum(discard), ' genes related to cell cycle \n')
  # genes.discard = rownames(diff)[which(discard)]
  aa <- RunPCA(aa, 
               features = hvgs, 
               verbose = FALSE, 
               weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) + NoLegend()
  p11 = FeaturePlot(aa, features = 'nCount_RNA')
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(outDir_cc, '/UMAP_RA_condition_clustering_cellcycle.rmCellCycleGenes_', 
                           select.method, '.pdf'),  width = 12, height = 6)
  
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(outDir_cc, '/featureExamples_UMAP_cellcycle.rmCellCycleGenes.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = aa, 
                           features = gene_examples[!is.na(match(gene_examples, rownames(aa)))])
  
  dev.off()
  
  # Idents(aa) = aa$seurat_clusters
  # all.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
  # all.markers %>%
  #   group_by(cluster) %>%
  #   top_n(n = 10, wt = avg_log2FC) -> top10
  # 
  # DoHeatmap(aa, features = top10$gene) + NoLegend()
  # ggsave(filename = paste0(outDir_cc, '/Heatmap_clusterMarkers_', select.method, '.pdf'), 
  #        width = 14, height = 20)
  
  
}
