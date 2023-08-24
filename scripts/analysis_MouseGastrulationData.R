##########################################################################
##########################################################################
# Project: RA competence project
# Script purpose: check the Pax6 and Fox2 co-expression in mouse data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Mar 29 11:02:25 2023
##########################################################################
##########################################################################

########################################################
########################################################
# Section I : import the data 
# the data from the following paper:
# https://www.nature.com/articles/s41586-019-0933-9
# R pacakge of the data were found: 
# https://bioconductor.org/packages/release/data/experiment/html/MouseGastrulationData.html
########################################################
########################################################

#BiocManager::install("MouseGastrulationData")
rm(list = ls())

version.analysis = '_MouseGastrulationData/'

resDir = paste0("../results/dataset_scRNAseq", version.analysis)
RdataDir = paste0(resDir, 'Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)


functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source(paste0(functionDir, '/functions_scRNAseq.R'))
source(paste0(functionDir, '/functions_Visium.R'))

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
options(future.globals.maxSize = 80000 * 1024^2)
mem_used()

species = '_mouseGastrulation'

library(MouseGastrulationData)

##########################################
# make plots of processed data
##########################################
head(AtlasSampleMetadata, n = 3)

sce <- EmbryoAtlasData(type = 'processed', samples = NULL)
sce

saveRDS(sce, file = paste0(RdataDir, 'EmbryoAtlasData_all36sample.rds'))

head(rowData(sce))
rownames(sce) = rowData(sce)$SYMBOL

head(colData(sce))

#exclude technical artefacts
singlets <- which(!(colData(sce)$doublet | colData(sce)$stripped))

plot(
  x = reducedDim(sce, "umap")[singlets, 1],
  y = reducedDim(sce, "umap")[singlets, 2],
  col = EmbryoCelltypeColours[colData(sce)$celltype[singlets]],
  pch = 19,
  xaxt = "n", yaxt = "n",
  xlab = "UMAP1", ylab = "UMAP2"
)

sce = sce[, singlets]
EmbryoCelltypeColours = EmbryoCelltypeColours[colData(sce)$celltype[singlets]]

sce = scuttle::logNormCounts(sce)
rownames(sce) = make.unique(rownames(sce))
srat = as.Seurat(sce, counts = "counts", data = "logcounts")

rm(sce)

DimPlot(srat, reduction = 'umap', cols = EmbryoCelltypeColours, group.by = 'celltype', 
        label = TRUE, repel = TRUE,
        raster=FALSE)

ggsave(filename = paste0(resDir, '/umap_celltypes.pdf'), width = 20, height = 10)

saveRDS(srat, file = paste0(RdataDir, 'seuratObject_EmbryoAtlasData_all36sample.rds'))

##########################################
# check the Pax6 and FoxA2 co-expression
##########################################
scrat = readRDS(file = paste0(RdataDir, 'seuratObject_EmbryoAtlasData_all36sample.rds'))
aa = readRDS('/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/scRNAseq_mNT/saved_seuratObj/seuratObject_EmbryoAtlasData_all36sample.rds')

DimPlot(aa, reduction = 'umap', 
        #cols = EmbryoCelltypeColours, 
        group.by = 'celltype', 
        label = TRUE, repel = TRUE,
        raster=FALSE)

p1 = VlnPlot(srat, features = c('Pax6'), group.by = 'celltype') + NoLegend()
p2 = VlnPlot(srat, features = 'Foxa2', group.by = 'celltype') + NoLegend()

p1 / p2

ggsave(filename = paste0(resDir, '/Vlnplot_Pax6_Foxa2_celltypes.pdf'), width = 14, height = 10)


FeaturePlot(srat, features = c('Pax6', 'Foxa2'), blend = TRUE, raster = TRUE, order = TRUE)

ggsave(filename = paste0(resDir, '/Featureplots_Pax6_Foxa2_blended_ordered.pdf'), width = 20, height = 8)


## subset potential clusters
Idents(srat) = as.factor(srat$celltype)

celltypes_sels = c('Surface ectoderm', 'Spinal cord', 'Rostral neurectoderm',
                   'NMP', 'Forebrain/Midbrain/Hindbrain', 
                   'Caudal neurectoderm', 'Caudal epiblast')

sub.obj = subset(srat, idents = celltypes_sels)

Idents(sub.obj) = as.factor(sub.obj$celltype)

DimPlot(sub.obj, reduction = 'umap', group.by = 'celltype', 
        label = TRUE, repel = TRUE,
        raster=FALSE)


sub.obj <- RunUMAP(sub.obj, reduction = "pca.corrected", dims = 1:50, 
                   n.neighbors = 30, min.dist = 0.1)

p1 = DimPlot(sub.obj, reduction = 'umap', group.by = 'celltype', 
        label = TRUE, repel = TRUE,
        raster=FALSE)
p2 = DimPlot(sub.obj, reduction = 'umap', group.by = 'stage', 
             label = TRUE, repel = TRUE,
             raster=FALSE)

p1 + p2

ggsave(filename = paste0(resDir, '/umap_subsetting_celltypes_stages.pdf'), width = 25, height = 8)

FeaturePlot(sub.obj, features = c('Pax6', 'Foxa2'), blend = TRUE, raster = TRUE, order = TRUE)
ggsave(filename = paste0(resDir, '/Featureplots_Pax6_Foxa2_blended_ordered_subsetting.pdf'), 
       width = 20, height = 8)




