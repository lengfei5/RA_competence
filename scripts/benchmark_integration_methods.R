##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: to compare and benchmark the integration and projection (mapping methods)
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Nov 23 11:59:22 2023
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_MouseGastrulationData/'

#resDir = paste0("../results/dataset_scRNAseq", version.analysis)
resDir = paste0('../results/scRNAseq_R13547_10x_mNT_20220813/mapping_to', version.analysis)
RdataDir = paste0(resDir, 'Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source(paste0(functionDir, '/functions_scRNAseq.R'))
source(paste0(functionDir, '/functions_Visium.R'))
source(paste0(functionDir, '/functions_dataIntegration.R'))

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(DropletUtils)
library(edgeR)
library(future)
library(data.table)
library(tidyverse)

options(future.globals.maxSize = 160000 * 1024^2)
mem_used()

species = '_mouseGastrulation'

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

## subset our scRNA-seq data 
levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA",
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA")

#levels_sels = c("day2_beforeRA",  "day2.5_RA", "day3_RA.rep1", "day3.5_RA",
#                "day4_RA", "day5_RA", "day6_RA")

cols_sel = cols[match(levels_sels, names(cols))]


##########################################
# import the mNT scRNAseq data and reference
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition.rds'))

ref = readRDS(file = paste0(RdataDir,  
                            'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v2.rds'))

cols_mouse = sapply(ref$colour, function(x) {paste0('#', x, collapse = '')})
names(cols_mouse) =ref$celltype
cols_mouse = cols_mouse[match(unique(names(cols_mouse)), names(cols_mouse))]

data_version = 'mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes'

##########################################
# test Seurat data projection
##########################################
mapping_method = "seurat_projection"

outDir = paste0(resDir,  data_version, '/', mapping_method, '/')
system(paste0('mkdir -p ', outDir))

features.common = intersect(rownames(aa), rownames(ref))

aa = subset(aa, features = features.common)
ref = subset(ref, features = features.common)

aa$dataset = 'mNT'
aa$stage = aa$condition
aa$sequencing.batch = 'mNT'
ref$dataset = 'ref'

aa$celltype = paste0('mNT_', aa$condition)

##########################################
# test projection method (slight different from integration)
# https://satijalab.org/seurat/articles/integration_mapping.html (original code)
##########################################
pancreas.anchors <- FindTransferAnchors(reference = pancreas.ref, query = pancreas.query, dims = 1:30,
                                        reference.reduction = "pca")

pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")


