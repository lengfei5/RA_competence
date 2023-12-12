##########################################################################
##########################################################################
# Project: RA competence project
# Script purpose: check the Pax6 and Fox2 co-expression in mouse data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Mar 29 11:02:25 2023
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
# import mNT scRNA-seq data
##########################################
aa =  readRDS(file = paste0('../results/Rdata/',  
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', 'mNT_scRNAseq',
                            '_R13547_10x_mNT_20220813', '.rds'))


Idents(aa) = factor(aa$condition)
aa = subset(aa, idents = levels_sels)

# downsample for each condition
aa = subset(x = aa, downsample = 1000)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition.rds'))

########################################################
########################################################
# Section I: Explore the mapping our mNT scRNA-seq data to mouse gastrulation atlas
# Reference choices
# 1) full atlas from Marioni2019 or Chan2019 or integrated Marioni.Chan
# 2) selected relevant cell types from those full atlas 
# Methods choices: seurat_CCA, seurat_RPCA, harmony, fastMNN, scanorama, scVI
########################################################
########################################################

##########################################
# import the mNT data and reference
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition_reclustered.rds'))

mapping_method = "seurat_rpca"


ref = readRDS(file = paste0(RdataDir,  
                             'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v2.rds'))

cols_mouse = sapply(ref$colour, function(x) {paste0('#', x, collapse = '')})
names(cols_mouse) =ref$celltype
cols_mouse = cols_mouse[match(unique(names(cols_mouse)), names(cols_mouse))]

data_version = 'mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes'

outDir = paste0(resDir, mapping_method, '/',  data_version, '/')
system(paste0('mkdir -p ', outDir))

##########################################
# test Seurat data integration
##########################################

features.common = intersect(rownames(aa), rownames(ref))

aa = subset(aa, features = features.common)
ref = subset(ref, features = features.common)

aa$dataset = 'mNT'
aa$stage = aa$condition
aa$sequencing.batch = 'mNT'
ref$dataset = 'ref'

aa$celltype = paste0('mNT_', aa$condition)

##########################################
# calculate similarity before data integration
##########################################
# aa = FindVariableFeatures(aa, selection.method = "vst")
# aa <- ScaleData(aa, verbose = FALSE)
# aa <- RunPCA(aa, npcs = 50, verbose = FALSE)
# ElbowPlot(aa, ndims = 50)
# 
# aa <- RunUMAP(aa, reduction = "pca", dims = 1:30, n.neighbors = 50, 
#               min.dist = 0.2) 
# DimPlot(aa, group.by = 'condition')
# 
# aa <- FindNeighbors(aa, reduction = "pca", dims = 1:20)
# aa <- FindClusters(aa, resolution = 0.7)
# aa$clusters = aa$seurat_clusters
#saveRDS(aa, file = paste0(RdataDir, 
#                          'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition_reclustered.rds'))

p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)
p1 + p2

ggsave(paste0(outDir, '/umap_query__conditions_clusters.pdf'), 
       width = 16, height = 6)

FeaturePlot(aa, features = 'Foxa2')

source(paste0(functionDir, '/functions_dataIntegration.R'))

ref = FindVariableFeatures(ref, selection.method = 'vst', nfeatures = 1000)

cc = c(3, 6, 5, 1, 7, 4, 11, 10)

for(n in 1:length(cc))
{
  # n = 4
  cat(n, ' -- ', cc[n], '\n')
  
  subs = subset(aa, cells = colnames(aa)[which(aa$clusters == cc[n])]);
  px = calculate_similarity_query_ref(query = subs, 
                                      ref = ref, 
                                      nHVGs = 1000, 
                                      method = c("spearman"),
                                      group.by = 'celltype')
  
  pdfname = paste0(outDir, '/spearman_similarity_withRefCelltypes_', cc[n], '.pdf')
  
  pdf(pdfname, width=16, height = 8)
  plot(px)
  
  dev.off()
  
}

##########################################
# test integration method
##########################################
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

saveRDS(cols_used, file = paste0(outDir, '/integrated_mNT_mouseGastrulation_colorsUsed.rds'))

# Visualization
DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
        repel = TRUE, raster=FALSE, cols = c(cols_mouse, cols_sel)) 

ggsave(paste0(outDir, '/Integration_mNT_celltypes_noref.batchcorrection.pdf'), 
       width = 16, height = 8)

#DimPlot(ref.combined, reduction = "umap")
saveRDS(ref.combined, file = paste0(outDir, '/integrated_mNT_mouseGastrulation_SeuratRPCA.rds'))

ref.combined = readRDS(file = paste0(outDir, '/integrated_mNT_mouseGastrulation_SeuratRPCA.rds'))
cols_used = readRDS(file = paste0(outDir, '/integrated_mNT_mouseGastrulation_colorsUsed.rds'))

# Visualization
DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
        repel = TRUE, raster=FALSE, cols = cols_used) 

ggsave(paste0(outDir, '/Integration_mNT_celltypes_noref.batchcorrection.pdf'), 
       width = 16, height = 8)

DimPlot(ref.combined, reduction = "umap", group.by = "dataset", raster=FALSE)
ggsave(paste0(outDir, '/Integration_dataset.pdf'), 
       width = 16, height = 8)

ggs = c('Pax6', 'Foxa2', 'Pou5f1', 'Sox17', 'Sox1', 'Sox2', 
        "Ifitm1", "T",
"Krt8", 'Tuba1a', "Eno1", 'Krt18','Foxj1', 'Sp5', 'Noto', 'Pou5f1','Foxa2','Cdx2', 'Gsta4', 'Sox9','Fst',
'Nog', 'Shh', 'Slc2a1', 'Cxx1b','Igfbp2','Epcam', 'Lhx1','Ptch1',
'Flrt3','Foxp1', 'Hoxb1', 'Dnmt3b', 'Nr6a1','Cdh2','Rspo3','Zfp503','Fgf5')

ggs = ggs[which(!is.na(match(ggs, rownames(ref.combined))))]

pdf(paste0(outDir, '/FeaturePlot_Markers.pdf'),
    width =10, height = 8, useDingbats = FALSE)
for(n in 1:length(ggs))
{
  cat(n, '--', ggs[n], '\n')
  p1 = FeaturePlot(ref.combined, features = ggs[n], min.cutoff = 'q5')
  #FeaturePlot(ref.combined, features = 'Foxa2', min.cutoff = 'q5')
  #FeaturePlot(ref.combined, features = 'Sox17', min.cutoff = 'q5')
  plot(p1)
  
}

dev.off()

DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE, split.by = 'dataset', 
        cols = cols_used,
        repel = TRUE, raster=FALSE) + NoLegend()

ggsave(paste0(outDir, '/Integration_celltypes_split.dataset.pdf'), 
       width = 24, height = 8)

DimPlot(ref.combined, reduction = "umap", group.by = "stage", label = TRUE,
        repel = TRUE, raster=FALSE)

ggsave(paste0(outDir, '/Integration_stage.pdf'), 
       width = 16, height = 8)

##########################################
# test the similarity calculation between mNT cells and cell types in the reference
# after data integration
# 
##########################################
DefaultAssay(ref.combined) = 'integrated'
#ref.combined = FindVariableFeatures(ref.combined, selection.method = 'vst', nfeatures = 1000)

refs_subs = subset(ref.combined, cells = colnames(ref.combined)[grep('mouseGastrulation', 
                                                                    colnames(ref.combined))]);
cc = c(3, 6, 5, 1, 7, 4, 11, 10)
source(paste0(functionDir, '/functions_dataIntegration.R'))

for(n in 1:length(cc))
{
  # n = 1
  cat('cluster -- ', cc[n], '\n')
  
  cells = colnames(aa)[which(aa$clusters == cc[n])]
  mm1 = match(colnames(ref.combined), paste0('mNT_', cells))
  subs = subset(ref.combined, cells = colnames(ref.combined)[!is.na(mm1)]);
  
  px = calculate_similarity_query_ref(query = subs, 
                                      ref = refs_subs, 
                                      assay_use = 'integrated',
                                      find_hvg = FALSE,
                                      method = c("spearman"),
                                      group.by = 'celltype')
  
  pdfname = paste0(outDir, '/spearman_similarity_withRefCelltypes_dataIntegration_', cc[n], '_test_1.pdf')
  
  pdf(pdfname, width=16, height = 8)
  plot(px)
  
  dev.off()
  
}

##########################################
# clustering the combined data
##########################################
ref.combined = readRDS(file = paste0(outDir, '/integrated_mNT_mouseGastrulation_SeuratRPCA.rds'))

ElbowPlot(ref.combined, ndims = 50)
ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
ref.combined <- FindClusters(ref.combined, resolution = 1.0)

DimPlot(ref.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
        repel = TRUE, raster=FALSE)

ggsave(paste0(outDir, '/integrated_ref_mNT_27clusters_resolution1.0.pdf'), 
       width = 16, height = 10)

cluster19.markers <- FindMarkers(ref.combined, ident.1 = 19)
head(cluster19.markers, n = 10)

cluster.markers <- FindMarkers(ref.combined, ident.1 = 19, ident.2 = c(16, 24))
head(cluster.markers, n = 10)

markers <- FindAllMarkers(ref.combined, only.pos = TRUE)

save(cluster.markers, cluster19.markers, markers, 
     file = paste0(outDir, '/marker_Genes.Rdata'))

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 30) %>%
  ungroup() -> top10

DoHeatmap(ref.combined, features = top10$gene) + NoLegend()
ggsave(paste0(outDir, '/markerGenes_all.clusters.for.cluster19_top30.pdf'), 
       width = 16, height = 30)



########################################################
########################################################
# Section II: first check Pax6-FoxA2 positive cells
# 
########################################################
########################################################
srat = readRDS(file = paste0(RdataDir,  
                             'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v2.rds'))

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


