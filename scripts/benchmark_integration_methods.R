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
                            'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v3.rds'))

cols_mouse = sapply(ref$colour, function(x) {paste0('#', x, collapse = '')})
names(cols_mouse) =ref$celltype
cols_mouse = cols_mouse[match(unique(names(cols_mouse)), names(cols_mouse))]

data_version = 'mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes'


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
# did not work well and don't know the reason
##########################################
Test_Seurat_projection = FALSE 
if(Test_Seurat_projection){
  mapping_method = "seurat_projection"
  
  outDir = paste0(resDir,  data_version, '/', mapping_method, '/')
  system(paste0('mkdir -p ', outDir))
  
  
  ElbowPlot(ref, ndims = 50, reduction = 'pca')
  ref = RunUMAP(ref, reduction = "pca", dims = 1:30, n.neighbors = 30, 
                min.dist = 0.1, return.model = TRUE) 
  
  DimPlot(ref, reduction = "umap", 
          group.by = "celltype", label = TRUE,
          repel = TRUE, raster=FALSE, cols = cols_mouse) 
  
  anchors <- FindTransferAnchors(reference = ref, 
                                 query = aa, 
                                 dims = 1:50,
                                 normalization.method = "LogNormalize",
                                 reference.reduction = "pca",
                                 reduction = "rpca"
  )
  
  query <- MapQuery(anchorset = anchors, 
                    reference = ref, 
                    query = aa,
                    #refdata = list(celltype = "celltype"), 
                    reference.reduction = "pca", 
                    reduction.model = "umap")
  
  p1 <- DimPlot(ref, reduction = "umap", group.by = "celltype", 
                label = TRUE, label.size = 3,
                repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
  p2 <- DimPlot(query, reduction = "ref.umap", group.by = "condition", label = TRUE,
                label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
  p1 + p2
  
}

##########################################
# test seurat rpca but trying to use the orignal mnn reudction 
##########################################
Test_Seurat_rpca_bc = FALSE
if(Test_Seurat_rpca_bc){
  mapping_method = "seurat_rpca_ref.mnn"
  
  outDir = paste0(resDir,  data_version, '/', mapping_method, '/')
  system(paste0('mkdir -p ', outDir))
  
  refs.merged = merge(aa, y = ref, add.cell.ids = c("mNT", "mouseGastrulation"), project = "RA_competence")
  ref.list <- SplitObject(refs.merged, split.by = "dataset")
  
  rm(list = c('refs.merged')) # remove big seurat objects to clear memory
  
  # normalize and identify variable features for each dataset independently
  ref.list <- lapply(X = ref.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize")
    x <- FindVariableFeatures(x, selection.method = "vst")
    
  })
  
  # select features that are repeatedly variable across datasets for integration run PCA on each
  # dataset using these features
  features <- SelectIntegrationFeatures(object.list = ref.list)
  
  ref.list <- lapply(X = ref.list, FUN = function(x) {
    x <- ScaleData(x, features = features.common, verbose = TRUE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    
  })
  
  ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                        anchor.features = features, 
                                        #reference = c(2),
                                        #reduction = "cca", 
                                        reduction = 'rpca',
                                        k.anchor = 5,
                                        dims = 1:50)
  
  rm(ref.list)
  
  # this command creates an 'integrated' data assay
  ref.combined <- IntegrateData(anchorset = ref.anchors, features.to.integrate = features.common, 
                                dims = 1:50) ## take ~100G memory
  
  rm(ref.anchors)
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(ref.combined) <- "integrated"
  
  ref.combined <- ScaleData(ref.combined, verbose = FALSE)
  ref.combined <- RunPCA(ref.combined, npcs = 50, verbose = FALSE)
  
  ElbowPlot(ref.combined, ndims = 50)
  
  kk = which(ref.combined$dataset == 'mNT') 
  ref.combined$celltype[kk] = paste0('mNT_', ref.combined$condition[kk])
  names(cols_sel) = paste0('mNT_', names(cols_sel))
  
  cols_used = c(cols_mouse, cols_sel)
  #ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
  #ref.combined <- FindClusters(ref.combined, resolution = 0.2)
  
  ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:50, n.neighbors = 50, 
                          min.dist = 0.2) 
  
  
}




