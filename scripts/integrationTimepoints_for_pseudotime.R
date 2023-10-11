##########################################################################
##########################################################################
# Project: RA competence project
# Script purpose: try the data integration across time points to estimate the pseudotime in palantir or scanpy
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Oct 11 10:58:41 2023
##########################################################################
##########################################################################
suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome)
  library(eisaR)
  library(GenomicFeatures)
  library(SummarizedExperiment)
  library(tximeta)
  library(rjson)
  library(reticulate)
  library(SingleCellExperiment)
  library(scater)
})

names(cols) = levels

levels_sels = c("day2_beforeRA", 
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA", 'day6_RA')

cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/RA_symetryBreaking/dataIntegration_timePoints_4pseudotime/')
system(paste0('mkdir -p ', outDir))


##########################################
# import data and select the samples  
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d6_',
                           species, version.analysis, '.rds'))

Idents(aa) = aa$condition
p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)

p1 + p2

ggsave(filename = paste0(outDir, 'UMAP_RAtreatment_d2.to.d6_.pdf'), 
       width = 18, height = 8)

## remove the cluster 8 and 9 mainly mature neurons and also day6_RA
aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '8' & aa$celltypes != '9' & 
                                             aa$condition != 'day6_RA')])


# rerun the umap 
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)

ggsave(filename = paste0(outDir, 'UMAP_RAtreatment_d2.to.d5.noMatureNeurons.pdf'), 
       width = 10, height = 6)


Discard_cellCycle.corrrelatedGenes = TRUE
if(Discard_cellCycle.corrrelatedGenes){
  library(scater)
  Idents(aa) = aa$condition
  
  # Identifying the likely cell cycle genes between phases,
  # using an arbitrary threshold of 5%.
  scaledMatrix = GetAssayData(aa, slot = c("scale.data"))
  
  diff <- getVarianceExplained(scaledMatrix, data.frame(phase = aa$Phase))
  diff = data.frame(diff, gene = rownames(diff))
  diff = diff[order(-diff$phase), ]
  
  hist(diff$phase, breaks = 100); 
  abline(v = c(1:5), col = 'red')
  
  rm(scaledMatrix)
  
  genes_discard = diff$gene[which(diff$phase>5)]
  cat(length(genes_discard), 'genes to discard \n')
  
  print(intersect(genes_discard, gene_examples))
  
  aa = subset(aa, features = setdiff(rownames(aa), genes_discard))
  
}

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)
Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)


aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.3)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)


#DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)
#ggsave(filename = paste0(outDir, 'UMAP_RAsymmetryBreaking.onlyday3rep1_3000HVGs_noweighted.byvarPCA_',
#                         '30pcs_30neighbors_minDist0.1.pdf'), 
#       width = 10, height = 8)

# quickly run clustering
#ElbowPlot(aa, ndims = 50)
aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

ggsave(filename = paste0(outDir, 
                         'UMAP_RAsymmetryBreaking.onlyday3rep1_timePoints_clustering.res0.7_noRAday6',
                         version.analysis, '.pdf'), 
       width = 18, height = 8)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)

aa$clusters = aa$seurat_clusters


########################################################
########################################################
# Section :
# 
########################################################
########################################################

