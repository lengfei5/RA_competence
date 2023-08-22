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
library(future)

options(future.globals.maxSize = 80000 * 1024^2)

species = 'mNT_scRNAseq'

### Parameters for this project
levels = c("day0_beforeRA", "day1_beforeRA", 
           "day2_beforeRA",
           "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", 'day3.5_RA',
           "day4_RA", "day5_RA", "day6_RA",
           "day2.5_noRA", "day3_noRA", 'day3.5_noRA', "day4_noRA", "day5_noRA", "day6_noRA")

##########################################
# main function 
##########################################
aa = readRDS(file = paste0(RdataDir, 
              'seuratObject_merged_cellFiltered_doubletFinderOut.v2_geneFiltered.15kGene_regressed.nCounts_', 
                           species, version.analysis, '.rds'))

Idents(aa) = factor(aa$condition, levels = levels)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)

aa = subset(aa, DF_out == 'Singlet')

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)


aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.2)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)

markers = FindAllMarkers(aa, only.pos = TRUE,  #min.pct = 0.25, 
                         logfc.threshold = 0.25)

saveRDS(markers, 
        file = paste0(RdataDir, 
  'seuratObject_merged_cellFiltered_doubletRemoved_geneFiltered.15kGene_regressed.nCounts_allMarkers_v1', 
  species, version.analysis, '.rds'))
