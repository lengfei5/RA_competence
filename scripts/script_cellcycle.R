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
# import the seurat object
##########################################
outDir = paste0(resDir, '/cellCycle_correction/')
system(paste0('mkdir -p ', outDir))

aa =  readRDS(file = paste0(RdataDir, 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', 
                            species, version.analysis, '.rds'))

Idents(aa) = factor(aa$condition, levels = levels)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)

s.genes <- firstup(cc.genes$s.genes)
s.genes[which(s.genes == 'Mlf1ip')] = 'Cenpu'
g2m.genes <- firstup(cc.genes$g2m.genes)

##########################################
# Regression both S.score and G2M score
##########################################
Regress.S.G2M.scores = FALSE
if(Regress.S.G2M.scores){
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
  
}

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_',
                           'regressout.nCounts.S.G2M.score_annot.v2_', 
                           species, version.analysis, '.rds'))


DimPlot(aa, label = FALSE, repel = TRUE, group.by = 'Phase', raster=FALSE)

# DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)
ggsave(filename = paste0(outDir, 'UMAP_cellcycle_regressed.S.G2M.scores.pdf'), 
       width = 12, height = 8)


##########################################
# regress the difference between S.score and G2M score 
##########################################
Regress.S.G2M.score.difference = FALSE
if(Regress.S.G2M.score.difference){
  aa$CC.Difference <- aa$S.Score - aa$G2M.Score
  aa <- ScaleData(aa, vars.to.regress = c('nCount_RNA', "CC.Difference"), features = rownames(aa))
  
  # cell cycle effects strongly mitigated in PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
  # ElbowPlot(aa, ndims = 30)
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
  
  # when running a PCA on cell cycle genes, actively proliferating cells remain distinct from G1
  # cells however, within actively proliferating cells, G2M and S phase cells group together
  #aa <- RunPCA(aa, features = c(s.genes, g2m.genes))
  #DimPlot(aa)
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_',
                            'regressout.nCounts.S.G2M.scoreDiff_annot.v2_', 
                            species, version.analysis, '.rds'))
  
}

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_',
                           'regressout.nCounts.S.G2M.scoreDiff_annot.v2_', 
                           species, version.analysis, '.rds'))

DimPlot(aa, label = FALSE, repel = TRUE, group.by = 'Phase', raster=FALSE)

# DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)
ggsave(filename = paste0(outDir, 'UMAP_cellcycle_regressed.S.G2M.score.Diff.pdf'), 
       width = 12, height = 8)


##########################################
# Removing cell cycle-related genes in the HVGs 
##########################################
#library(scRNAseq)
#require(scran)
library(scater)
#require(batchelor)

Idents(aa) = aa$condition
p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)

# Identifying the likely cell cycle genes between phases,
# using an arbitrary threshold of 5%.

scaledMatrix = GetAssayData(aa, slot = c("scale.data"))

diff <- getVarianceExplained(scaledMatrix, data.frame(phase = aa$Phase))

hist(diff, breaks = 100); abline(v = c(1:5), col = 'red')
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)

for(diff_cutoff in c(seq(0.1, 0.8, by=0.2), 1:5))
{
  
  discard <- diff > diff_cutoff
  summary(discard)
  
  cat('cut off -- ', diff_cutoff, '--', sum(discard), ' genes related to cell cycle \n')
  genes.discard = rownames(diff)[which(discard)]
  
  aa <- RunPCA(aa, 
               features = setdiff(VariableFeatures(object = aa), genes.discard), 
               verbose = FALSE, 
               weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.2)
  
  p1 = DimPlot(aa, label = FALSE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  
  p0 + p1
  
  # DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)
  ggsave(filename = paste0(outDir, 'UMAP_cellcycleGened_discard_percentVariance_', diff_cutoff, '.pdf'), 
         width = 20, height = 8)
  
}
