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

options(future.globals.maxSize = 80000 * 1024^2)

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
# tfs and sps annotations 
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

for(cc in conditions){
  cc = c('day2_beforeRA')
  # cc = c('day2.5_RA')
  # cc = c('day3_RA.rep1')
  # cc = c('day3.5_RA')
  # cc = c('day4_RA')
  # cc = c('day5_RA')
  
  cc = c('day3_RA.rep1', 'day3.5_RA')
  
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
  plot_manyFeatures_seurat(seurat_obj = aa, features = gene_examples)
  
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
  plot_manyFeatures_seurat(seurat_obj = xx, features = gene_examples)
  
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
  plot_manyFeatures_seurat(seurat_obj = aa, features = gene_examples)
  
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

