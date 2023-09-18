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


##########################################
# make plots of processed data
##########################################
Load_process_MouseGastrulation = FALSE
if(Load_process_MouseGastrulation){
  library(MouseGastrulationData)
  head(AtlasSampleMetadata, n = 3)
  
  sce <- EmbryoAtlasData(type = 'processed', samples = NULL)
  sce
  
  saveRDS(sce, file = paste0(RdataDir, 'EmbryoAtlasData_all36sample.rds'))
  
  sce = readRDS(paste0('../results/dataset_scRNAseq_MouseGastrulationData/Rdata/', 
                       'EmbryoAtlasData_all36sample.rds'))
  
  head(rowData(sce))
  rownames(sce) = rowData(sce)$SYMBOL
  
  head(colData(sce))
  
  #exclude technical artefacts
  singlets <- which(!(colData(sce)$doublet | colData(sce)$stripped))
  
  plot(
    x = reducedDim(sce, "umap")[singlets, 1],
    y = reducedDim(sce, "umap")[singlets, 2],
    col = sce$colour[singlets],
    pch = 19,
    xaxt = "n", yaxt = "n",
    xlab = "UMAP1", ylab = "UMAP2"
  )
  
  sce = sce[, singlets]
  EmbryoCelltypeColours = EmbryoCelltypeColours[colData(sce)$celltype[singlets]]
  
  sce = scuttle::logNormCounts(sce)
  rownames(sce) = make.unique(rownames(sce))
  srat = as.Seurat(sce, counts = "counts",  assay = NULL)
  
  rm(sce)
  
  DimPlot(srat, reduction = 'umap', cols = EmbryoCelltypeColours, group.by = 'celltype', 
          label = TRUE, repel = TRUE,
          raster=FALSE)
  
  ggsave(filename = paste0(resDir, '/umap_celltypes.pdf'), width = 20, height = 10)
  
  saveRDS(srat, file = paste0(RdataDir, 'seuratObject_EmbryoAtlasData_all36sample.rds'))
  
}

########################################################
########################################################
# Section I: process the mouse gastrulation data and first check Pax6-FoxA2 positive cells
# 
########################################################
########################################################
Modify.Assay.Name = FALSE
if(Modify.Assay.Name){
  srat = readRDS(file = paste0('../results/dataset_scRNAseq_MouseGastrulationData/Rdata/',
                               'seuratObject_EmbryoAtlasData_all36sample.rds'))
  
  ## change assay name
  adt.data <- GetAssayData(object =  srat[['originalexp']], slot = 'counts')
  srat[["RNA"]] <- CreateAssayObject(counts = adt.data )
  DefaultAssay(srat) <- "RNA"
  srat[['originalexp']] = NULL
  
  srat <- NormalizeData(srat, normalization.method = "LogNormalize")
  srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 3000)
  srat <- ScaleData(srat, verbose = FALSE)
  srat <- RunPCA(srat, verbose = FALSE)
  
  saveRDS(srat, file = paste0(RdataDir,  'seuratObject_EmbryoAtlasData_all36sample_RNAassay.rds'))
}

srat = readRDS(file = paste0(RdataDir,  'seuratObject_EmbryoAtlasData_all36sample_RNAassay.rds'))
xx = readRDS(file = paste0('../results/dataset_scRNAseq_MouseGastrulationData/Rdata/',
                      'seuratObject_EmbryoAtlasData_all36sample.rds'))

umap.embedding = xx@reductions$umap@cell.embeddings
umap.embedding = umap.embedding[match(colnames(srat), rownames(umap.embedding)), ]
srat[['umap']] = Seurat::CreateDimReducObject(embeddings=umap.embedding,
                                                 key='UMAP_',
                                                 assay='RNA')
rm(xx)
rm(umap.embedding)

srat = readRDS(file = paste0(RdataDir,  'seuratObject_EmbryoAtlasData_all36sample_RNAassay.rds'))
xx = readRDS(file = paste0('../results/dataset_scRNAseq_MouseGastrulationData/Rdata/',
                           'seuratObject_EmbryoAtlasData_all36sample.rds'))

umap.embedding = xx@reductions$umap@cell.embeddings
umap.embedding = umap.embedding[match(colnames(srat), rownames(umap.embedding)), ]
srat[['umap']] = Seurat::CreateDimReducObject(embeddings=umap.embedding,
                                              key='UMAP_',
                                              assay='RNA')
rm(xx)
rm(umap.embedding)

## filter unlikely celltypes in the reference
sels = grep('Erythroid|Blood|Allantois|mesoderm|Haemato|Cardiomy|Endothelium|Mesenchyme|ExE', srat$celltype, 
            invert = TRUE)
srat = subset(srat, cells = colnames(srat)[sels])

saveRDS(srat, file = paste0(RdataDir,  
                            'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes.rds'))

srat = readRDS(file = paste0(RdataDir,  
                             'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes.rds'))


sels = grep('Parietal', srat$celltype, 
            invert = TRUE)
srat = subset(srat, cells = colnames(srat)[sels])

saveRDS(srat, file = paste0(RdataDir,  
                            'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v2.rds'))



p1 = DimPlot(srat, reduction = 'umap', 
        #cols = EmbryoCelltypeColours, 
        group.by = 'celltype', 
        label = TRUE, repel = TRUE,
        raster=FALSE) 

p2 = DimPlot(srat, reduction = 'umap', 
        #cols = EmbryoCelltypeColours, 
        group.by = 'stage', 
        label = TRUE, repel = TRUE,
        raster=FALSE) 

p1 /p2

ggsave(filename = paste0(resDir, '/MouseGastrulation_celltypes_stage.pdf'), width = 18, height = 20)

##########################################
# check the Pax6 and FoxA2 co-expression
##########################################
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


########################################################
########################################################
# Section II: Explore the mapping between the mouse data and our scRNA-seq data
# 
########################################################
########################################################

##########################################
# test Seurat data integration
##########################################
data_version = "subsettingRef_mNT.noRA.RA.d2_d5_Harmony"

outDir = paste0(resDir, 'dataMapping_', data_version)
system(paste0('mkdir -p ', outDir))

srat = readRDS(file = paste0(RdataDir,  
                             'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v2.rds'))

cols_mouse = sapply(srat$colour, function(x) {paste0('#', x, collapse = '')})
names(cols_mouse) =srat$celltype
cols_mouse = cols_mouse[match(unique(names(cols_mouse)), names(cols_mouse))]

aa =  readRDS(file = paste0('../results/Rdata/',  
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', 'mNT_scRNAseq',
                            '_R13547_10x_mNT_20220813', '.rds'))

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

Idents(aa) = factor(aa$condition)
aa = subset(aa, idents = levels_sels)

# downsample for each condition
aa = subset(x = aa, downsample = 1000)

features.common = intersect(rownames(aa), rownames(srat))
aa = subset(aa, features = features.common)
srat = subset(srat, features = features.common)

ref = srat;
rm(srat)

aa$dataset = 'mNT'
aa$stage = aa$condition
aa$sequencing.batch = 'mNT'
ref$dataset = 'ref'

refs.merged = merge(aa, y = ref, add.cell.ids = c("mNT", "mouseGastrulation"), project = "RA_competence")


