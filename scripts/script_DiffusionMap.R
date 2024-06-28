##########################################################################
##########################################################################
# Project: run diffusion map
# Script purpose: 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Sep 29 10:59:04 2022
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R13547_10x_mNT_20220813'

resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R13547_10x'
functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(DropletUtils)
library(edgeR)
library(future)
options(future.globals.maxSize = 180000 * 1024^2)
mem_used()

species = 'mNT_scRNAseq'

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

##########################################
# main code:  test Diffusion Map reduction
##########################################
library(pheatmap)
library(RColorBrewer)
library(grid)
library(Seurat)
library(scater)
library(SingleCellExperiment)
library(scran)
library(destiny)

aa =  readRDS(file = paste0(RdataDir, 
                'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                'cellCycleScoring_annot.v1_', 
                species, version.analysis, '.rds'))

## filter the weird clusters identified in make_plots_v0.R
cells_2filter = readRDS(file = paste0(RdataDir, 'subObj_clusters_to_filter.rds'))
mm = match(colnames(aa), colnames(cells_2filter))
aa = subset(aa, cells = colnames(aa)[which(is.na(mm))])


Idents(aa) = factor(aa$condition, levels = levels)
aa$condition = factor(aa$condition, levels = levels)

# DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)
#aa = subset(aa, DF_out == 'Singlet')

levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA",
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA")

cols_sel = cols[match(levels(aa$condition), names(cols))]


aa = subset(aa, idents = levels_sels)
aa$condition = droplevels(aa$condition)

####
#### not use the processing steps from Seurat
####
#aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
#aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
#ElbowPlot(aa, ndims = 50)

#ll.pca = Embeddings(object = aa, reduction = "pca")

#aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
#ElbowPlot(aa, ndims = 50)

#ll.pca.notWeighted = Embeddings(object = aa, reduction = "pca")

#Idents(aa) = aa$condition
#aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.2)
#DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
#DimPlot(aa, reduction = 'pca', label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
metadata = aa@meta.data
sce = as.SingleCellExperiment(aa)
#rm(aa)
dec <- modelGeneVar(sce)


for(nb_features in c(3000, 5000, 8000))
{
  for(n_neighbors in c(100, 200, 300, 500))
  {
    top.hvgs <- getTopHVGs(dec, n=nb_features)
    
    sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 100)
    # reducedDimNames(sce)
    ll.pca = reducedDim(sce, 'PCA')[, c(1:50)]
    
    tic()
    dm <- DiffusionMap(ll.pca, sigma = 'local', k = n_neighbors, 
                       n_eigs = 50, distance = 'euclidean')
    
    saveRDS(dm, file = paste0(RdataDir, 'diffusionMap_firstBifurcation_RA.vs.noRA.d2tod5_sce.PCA_euclidean_',
                              'localSignal_nfeatures.', nb_features, '_nbNeighbors.', n_neighbors,
                              '.rds'))
    
    toc()
    
    tic()
    dm <- DiffusionMap(ll.pca, sigma = 'global', k = n_neighbors, n_eigs = 50, distance = 'euclidean')
    
    saveRDS(dm, file = paste0(RdataDir, 'diffusionMap_firstBifurcation_RA.vs.noRA.d2tod5_sce.PCA_euclidean_',
                              'globalSignal_nfeatures.', nb_features, '_nbNeighbors.', n_neighbors,
                              '.rds'))
    toc()
    
  }
}

