rm(list = ls())

version.analysis = '_R13547_10x_mNT_20220813'

resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R13547_10x'
functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(future)
options(future.globals.maxSize = 80000 * 1024^2)

#plan()
# change the current plan to access parallelization
#plan("multicore", workers = 32)
#plan()

species = 'mNT_scRNAseq'

### Parameters for this project
levels = c("day0_beforeRA", "day1_beforeRA", 
           "day2_beforeRA",
           "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", 'day3.5_RA',
           "day4_RA", "day5_RA", "day6_RA",
           "day2.5_noRA", "day3_noRA", 'day3.5_noRA', "day4_noRA", "day5_noRA", "day6_noRA")

##########################################
# script
##########################################
# doulet removed already here
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_merged_cellFiltered_doublet.rm_', 
                           species, version.analysis, '_v3.rds'))
Idents(aa) = factor(aa$condition, levels = levels)

aa = subset(aa, DF_out == 'Singlet')
aa = subset(aa, features = rownames(aa)[grep('^Rp[sl]|^mt-', rownames(aa), invert = TRUE)])

aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(aa)

aa <- ScaleData(aa, features = all.genes, vars.to.regress = 'nCount_RNA')

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
# ElbowPlot(aa, ndims = 30)

Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_', 
                          species, version.analysis, '_cbe_v4.rds'))
