##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: remove the cell cycle effect 
# ideas from http://bioconductor.org/books/3.12/OSCA/cell-cycle-assignment.html#removing-cell-cycle-effects
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Nov 30 11:39:57 2022
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

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(future)
options(future.globals.maxSize = 120000 * 1024^2)

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
# import the seurat object
##########################################
aa =  readRDS(file = paste0(RdataDir, 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', 
                            species, version.analysis, '.rds'))

Idents(aa) = factor(aa$condition, levels = levels)
#DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)

s.genes <- firstup(cc.genes$s.genes)
s.genes[which(s.genes == 'Mlf1ip')] = 'Cenpu'
g2m.genes <- firstup(cc.genes$g2m.genes)

##########################################
# Regression both S.score and G2M score
##########################################
aa <- ScaleData(aa, vars.to.regress = c('nCount_RNA', "S.Score", "G2M.Score"), 
                features = rownames(aa))

# Now, a PCA on the variable genes no longer returns components associated with cell cycle
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
# ElbowPlot(aa, ndims = 30)
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
#aa <- RunPCA(aa, features = c(s.genes, g2m.genes))
#DimPlot(aa)

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_',
                          'regressout.nCounts.S.G2M.score_annot.v2_', 
                          species, version.analysis, '.rds'))

##########################################
# regress the difference between S.score and G2M score 
##########################################
# aa$CC.Difference <- aa$S.Score - aa$G2M.Score
# aa <- ScaleData(aa, vars.to.regress = c('nCount_RNA', "CC.Difference"), features = rownames(aa))
# 
# # cell cycle effects strongly mitigated in PCA
# aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
# # ElbowPlot(aa, ndims = 30)
# aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
# 
# # when running a PCA on cell cycle genes, actively proliferating cells remain distinct from G1
# # cells however, within actively proliferating cells, G2M and S phase cells group together
# #aa <- RunPCA(aa, features = c(s.genes, g2m.genes))
# #DimPlot(aa)
# 
# saveRDS(aa, file = paste0(RdataDir, 
#                           'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_',
#                           'regressout.nCounts.S.G2M.scoreDiff_annot.v2_', 
#                           species, version.analysis, '.rds'))

##########################################
# Removing cell cycle-related genes in the HVGs 
##########################################
