##########################################################################
##########################################################################
# Project: RA competence
# Script purpose: make a template of Seurat plotting 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Mar  2 16:15:57 2023
##########################################################################
##########################################################################
rm(list = ls())

##########################################
# specify input and output folders
##########################################
RdataDir = '/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/scRNAseq_mNT/saved_seuratObj/'

outDir = '/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/scRNAseq_mNT/Hannahs_analysis'

if(!dir.exists(outDir)) dir.create(outDir)

##########################################
# load packages and functions
##########################################
#source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
#source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')
library(pryr) # monitor the memory usage
require(ggplot2)
#require(dplyr)
#require(stringr)
#require(tidyr)
library(Seurat)
#library(DropletUtils)
#library(edgeR)
library(future)
options(future.globals.maxSize = 80000 * 1024^2)
mem_used()

### global parameters and color coding 
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
# import the data either the full dataset (16 samples)
# or only RA 
##########################################
# all 16 samples
# aa = readRDS(file = paste0(RdataDir,
#                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
#                            'cellCycleScoring_annot.v1_savedUMAP.v1_',  
#                            'mNT_scRNAseq_R13547_10x_mNT_20220813.rds'))

# only RA samples incl. dya2_beforeRA
aa = readRDS(file = paste0(RdataDir,
                           'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_savedUMAP.subs.v2_', 
                           'mNT_scRNAseq_R13547_10x_mNT_20220813.rds'))

Idents(aa) = aa$condition
DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)

##########################################
# main analysis step: 
# - subsetting cells
# - Find variable genes
# - run PCA 
# - run UMAP
# - run clustering
##########################################
Further.subset.samples = FALSE
if(Further.subset.samples){
  levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA")
  data_version = "_d2_d2.5_d3_d3.5_d4_d5"
  
  names(cols) = levels
  cols_sel = cols[match(levels_sels, names(cols))]
  
  aa = subset(aa, idents = levels_sels)
  cell_sels = colnames(aa)[which(aa$celltypes != 'Neurons'|is.na(aa$celltypes))]
  aa = subset(aa, cells = cell_sels)
  
  aa$condition = droplevels(as.factor(aa$condition))
  Idents(aa) = aa$condition
  DimPlot(aa, cols = cols_sel, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)
  # rerun the umap 
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.05)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  # quickly run clustering
  #ElbowPlot(aa, ndims = 50)
  aa <- FindNeighbors(aa, dims = 1:20)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
}

##########################################
# a bunch of plotting functions from Seurat
# more functions can be found in Seurat vignette
# https://satijalab.org/seurat/articles/visualization_vignette.html
##########################################
p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)
p1 + p2

ggsave(filename = paste0(outDir, '/UMAP_RA_symmetryBreaking.pdf'), 
       width = 14, height = 6)

FeaturePlot(aa, features = c('Foxa2', 'Pax6', 'Pou5f1', 'Sox11', 'Sox17', 'Sox4', 'Sox1', 'Sox2'))
ggsave(filename = paste0(outDir, '/FeaturePlot_geneExamples.pdf'), 
       width = 10, height = 8)



