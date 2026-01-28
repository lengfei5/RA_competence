rm(list = ls())

version.analysis = '_R16597_mNT_10xmultiome_reseq_20240517'

resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/'

source(paste0(functionDir,  'functions_scATAC.R'))
source(paste0(functionDir, 'functions_scRNAseq.R'))
source(paste0(functionDir, 'functions_Visium.R'))

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
require(SeuratObject)
library(data.table)
library(DropletUtils)
library(edgeR)
library(future)
library(tictoc)

options(future.globals.maxSize = 200 * 1024^3)
set.seed(1234)
mem_used()

species = 'mNT_multiome'

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


## subset the mutliome conditions
levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA", "day3.5_RA",  "day4_RA", "day5_RA", 
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA")

cc = names(cols)
cc[which(cc == 'day3_RA.rep1')] = 'day3_RA'
names(cols) = cc

cols_sel = cols[match(levels_sels, names(cols))]


srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                                'DFoutSinglet_cell.gene.filtered.rds'))


##########################################
# main script
##########################################
srat_cr = NormalizeData(srat_cr, normalization.method = "LogNormalize", scale.factor = 10000)

srat_cr <- FindVariableFeatures(srat_cr, selection.method = "vst", nfeatures = 5000)

srat_cr$CC.Difference <- srat_cr$S.Score - srat_cr$G2M.Score
#srat_cr <- ScaleData(srat_cr, vars.to.regress = "CC.Difference", features = rownames(srat_cr))
all.genes <- rownames(srat_cr)
srat_cr <- ScaleData(srat_cr, features = all.genes, 
                     vars.to.regress = c('nFeature_RNA',  "percent.mt", "CC.Difference"))

saveRDS(srat_cr, file = paste0(RdataDir, 
                               'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                               'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT.ccDifference.rds'))

