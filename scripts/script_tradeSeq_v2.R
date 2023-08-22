##########################################################################
##########################################################################
# Project:
# Script purpose: run the tradeSeq 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Sat Oct 22 10:20:28 2022
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
options(future.globals.maxSize = 80000 * 1024^2)
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

names(cols) = levels
levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", "day3.5_RA",   "day4_RA", 
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA")

cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/RA.vs.noRA_firstBifurcation/')
system(paste0('mkdir -p ', outDir))


##########################################
# load the input files and test tradeSeq
##########################################
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(BiocParallel)

load(file = paste0(outDir, '/counts_pseudotime_cellWeights_for_tradeSeq.Rdata'))

#RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
#data(countMatrix, package = "tradeSeq")
counts <- as.matrix(counts)
#rm(countMatrix)
#data(crv, package = "tradeSeq")
#data(celltype, package = "tradeSeq")

set.seed(7)
index_sub = sample(c(1:ncol(counts)), 5000)


### Downstream of any trajectory inference method using pseudotime and cell weights

BPPARAM <- BiocParallel::bpparam()
BPPARAM # lists current options
BPPARAM$workers <- 32 # use 2 cores

# slow, but still ok. with 60K cells, it takes 10-15 mins for each knot.
tic()
set.seed(7)
icMat <- evaluateK(counts = counts[, index_sub], 
                   pseudotime = pseudotime[index_sub, ], 
                   cellWeights = cellWeights[index_sub, ],
                   k=3:20, 
                   nGenes = 500, 
                   verbose = TRUE, 
                   plot = FALSE,
                   parallel=FALSE
                   #BPPARAM = BPPARAM
                   )

saveRDS(icMat, file = paste0(outDir, 'tradeSeqDE_evaluateK_res_k3.20_500genes_5k.cells_v7.rds'))

toc()

