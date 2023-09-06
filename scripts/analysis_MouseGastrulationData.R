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
# Section I: Explore the mapping between the mouse data and our scRNA-seq data
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


##########################################
# test Seurat data integration
##########################################
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
sels = grep('Erythroid|Blood|Allantois|mesoderm|Haemato|Cardiomy|Endothelium|Mesenchyme', srat$celltype, 
            invert = TRUE)
srat = subset(srat, cells = colnames(srat)[sels])


aa =  readRDS(file = paste0('../results/Rdata/',  
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', 'mNT_scRNAseq',
                            '_R13547_10x_mNT_20220813', '.rds'))

## subset our scRNA-seq data 
levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA",
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA")

levels_sels = c("day2_beforeRA",  "day2.5_RA", "day3_RA.rep1", "day3.5_RA",
                "day4_RA", "day5_RA", "day6_RA")

data_version = "subsettingRef_mNT.RA.d2_d6_Harmony"
outDir = paste0(resDir, 'dataMapping_', data_version)
system(paste0('mkdir -p ', outDir))

Idents(aa) = factor(aa$condition)
aa = subset(aa, idents = levels_sels)

aa = subset(x = aa, downsample = 2000)

aa$dataset = 'mNT'
aa$stage = aa$condition
srat$dataset = 'ref'


#refs <- CreateAssayObject(counts = srat@assays[["originalexp"]]@counts)
#metadata = srat@meta.data
#refs = AddMetaData(refs,  metadata, col.name = NULL) 

#aa$annot.ref = aa$my_annot
#cms$annot.ref = cms$CellType
features.common = intersect(rownames(aa), rownames(srat))

refs.merged = merge(aa, y = srat, add.cell.ids = c("mNT", "mouseGastrulation"), project = "RA_competence")

Run_Harmony = FALSE
if(!Run_Harmony){
  refs.merged = NormalizeData(refs.merged, normalization.method = "LogNormalize")
  refs.merged <- FindVariableFeatures(refs.merged, selection.method = "vst")
  refs.merged <- ScaleData(refs.merged, verbose = TRUE)
  refs.merged <- RunPCA(refs.merged, verbose = TRUE)
  
  ## original code 
  # from http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/advanced.html
  #V <- harmony::cell_lines$scaled_pcs
  #V_cos <- cosine_normalize(V, 1)
  #meta_data <- harmony::cell_lines$meta_data  
  V = refs.merged@reductions$pca@cell.embeddings
  meta_data = refs.merged@meta.data
  
  set.seed(1)
  harmony_embeddings <- harmony::HarmonyMatrix(
    data_mat = V, ## PCA embedding matrix of cells
    meta_data = meta_data, ## dataframe with cell labels
    #vars_use = 'dataset',
    theta = 0.5, ## cluster diversity enforcement
    vars_use = 'dataset', ## variable to integrate out
    npcs = 30,
    nclust = 5, ## number of clusters in Harmony model
    max.iter.harmony = 10, ## stop after initialization
    return_object = FALSE, ## return the full Harmony model object
    do_pca = FALSE ## don't recompute PCs
    
  )
  
  p1 <- do_scatter(harmony_embeddings, meta_data, 'dataset') + 
    labs(title = 'Colored by dataset')
  p2 <- do_scatter(harmony_embeddings, meta_data, 'cell_type') + 
    labs(title = 'Colored by cell type')
  cowplot::plot_grid(p1, p2, nrow = 1)
  
  
  refs.merged[['harmony']] = Seurat::CreateDimReducObject(embeddings=harmony_embeddings,
                                                key='HARMONY_',
                                                assay='RNA')
  
  ref.combined <- RunUMAP(refs.merged, reduction = "pca", dims = 1:50, n.neighbors = 50, 
                          min.dist = 0.1) 
  
  
}


ref.list <- SplitObject(refs.merged, split.by = "dataset")

rm(list = c('refs.merged', 'aa', 'srat')) # remove big seurat objects to clear memory

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
  x <- RunPCA(x, features = features, verbose = TRUE)
  
})

ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                      anchor.features = features, 
                                      reference = c(2),
                                      reduction = "cca", 
                                      #k.anchor = 5,
                                      dims = 1:50)

rm(ref.list)

# this command creates an 'integrated' data assay
ref.combined <- IntegrateData(anchorset = ref.anchors, features.to.integrate = features.common, 
                              dims = 1:50) ## take ~100G memory

rm(ref.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(ref.combined) <- "integrated"

xx = DietSeurat(ref.combined, counts = TRUE, data = TRUE, scale.data = TRUE, assays = 'integrated')
xx@assays$integrated@counts = ref.combined@assays$RNA@counts

#saveRDS(xx, file = paste0(outDir, 
#                          '/Seurat.obj_mouseGastrulation_mNT_integrated.rds'))

# Run the standard workflow for visualization and clustering
ref.combined = readRDS(file =paste0(outDir, 
                                    '/Seurat.obj_mouseGastrulation_mNT_integrated.rds'))


ref.combined <- ScaleData(ref.combined, verbose = FALSE)
ref.combined <- RunPCA(ref.combined, npcs = 50, verbose = FALSE)

ElbowPlot(ref.combined, ndims = 50)

#ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
#ref.combined <- FindClusters(ref.combined, resolution = 0.2)
ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:50, n.neighbors = 50, 
                        min.dist = 0.1) 

#DimPlot(ref.combined, reduction = "umap")

kk = which(ref.combined$dataset == 'mNT') 
ref.combined$celltype[kk] = ref.combined$condition[kk]


# Visualization
p1 <- DimPlot(ref.combined, reduction = "umap", group.by = "dataset", raster=FALSE)
p2 <- DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
              repel = TRUE, raster=FALSE) + NoLegend() 

p1 / p2 

ggsave(paste0(outDir, '/Integration_dataset_celltypes.pdf'), 
       width = 16, height = 24)

DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE, split.by = 'dataset',
        repel = TRUE, raster=FALSE) + NoLegend()

ggsave(paste0(outDir, '/Integration_celltypes_split.dataset.pdf'), 
       width = 24, height = 10)

DimPlot(ref.combined, reduction = "umap", group.by = "stage", label = TRUE,
        repel = TRUE, raster=FALSE)

ggsave(paste0(outDir, '/Integration_stage.pdf'), 
       width = 16, height = 10)

