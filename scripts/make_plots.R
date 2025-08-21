##########################################################################
##########################################################################
# Project: RA competence
# Script purpose: make plots for figures 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jan 19 10:59:56 2023
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R13547_10x_mNT_20240522'
species = 'mNT_scRNAseq'

resDir = paste0("../results/figures_tables", version.analysis)
figureDir = '/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/plots_manuscript_1/'
RdataDir = '../results/Rdata/'
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

library(pryr) # monitor the memory usage
require(ggplot2)

require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(viridis)

#library(DropletUtils)
library(future)
options(future.globals.maxSize = 160000 * 1024^2)
mem_used()

levels = c("day0_beforeRA", "day1_beforeRA", 
           "day2_beforeRA",
           "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", 'day3.5_RA',
           "day4_RA", "day5_RA", "day6_RA",
           "day2.5_noRA", "day3_noRA", 'day3.5_noRA', "day4_noRA", "day5_noRA", "day6_noRA")

cols = readRDS(file = '../results/Rdata/color_scheme_4scRNAseq.rds')
load(file = '../results/Rdata/tfs_sps_geneExamples_4scRNAseq.Rdata')

feat_cols = c("#F0F0F0", "#EFFAB6", "#69C6BE", "#007BB7", "#121D60")

cols_cluster = c( "#7F7F7F", "#FFC000", "#70AD47", "#337f01", "#CD00CF", "#265401","#800080")

########################################################
########################################################
# Section 0: the starting data object 
# saved in the following directory and 
# also in "../results/Rdata/"
########################################################
########################################################
aa = readRDS(paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/',
                    'scRNAseq_mNT/saved_seuratObj/',
                    'RA_noRAsamples_d2_d6_UMAP_selectedParam.rds'))

DimPlot(aa, group.by = 'condition')

########################################################
########################################################
# Section I: overview of scRNA-seq data and feature highlight
# 
########################################################
########################################################
outDir = paste0(resDir, '/UMAP_allSamples/')
if(!dir.exists(outDir)) dir.create(outDir)

Save_weirdClusters = FALSE
if(Save_weirdClusters){
  bb = readRDS(file = paste0('../results/Rdata/', 
                             'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                             'cellCycleScoring_annot.v2_newUMAP_clusters_time_',
                             species, '_R13547_10x_mNT_20220813', '.rds'))
  
  p1 = DimPlot(bb, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(bb, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  p1 + p2
  
  clusters_weird = colnames(bb)[which(bb$clusters == '7'|bb$clusters == '8')]
  saveRDS(clusters_weird, file = paste0(RdataDir, 'cellNames_weirdClusters_7_8.rds'))
  
}


##########################################
# Identify clusters to discard from all samples 
##########################################
aa =  readRDS(file = paste0('../results/Rdata/', 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', 
                            species, '_R13547_10x_mNT_20220813', '.rds'))
Idents(aa) = factor(aa$condition, levels = levels)
table(aa$condition)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)

ggsave(filename = paste0(resDir, '/UMAP_overview_initial_beforeFiltering.pdf'), 
       width = 12, height = 8)

## reproduce the same umap
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.1)
DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)

# quickly run clustering
aa <- FindNeighbors(aa, dims = 1:30)

aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 1.0)
p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

ggsave(filename = paste0(outDir, '/UMAP_RA_condition_clusters_res2.0.pdf'), 
       width = 14, height = 6)

table(aa$seurat_clusters)

jj = match(cells_2filter, colnames(aa))
xx = table(aa$seurat_clusters[jj])
yy = xx/table(aa$seurat_clusters)
xx[order(-xx)]
yy[order(-yy)]

cells =  colnames(aa)[which(!is.na(match(aa$seurat_clusters, c(15))))]

mm = match(cells, cells_2filter)
length(which(!is.na(mm)))

DimPlot(aa, label = TRUE, repel = TRUE,  raster=FALSE,
        cells.highlight = cells,
        sizes.highlight = 0.2)

clusters_sel = c(15, 8, 19, 0, 9, 21, 3)
cells = colnames(aa)[which(!is.na(match(aa$seurat_clusters, clusters_sel)))]

# discard first the cells in cluster 23
aa = subset(aa, cells = cells)                               
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)


aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.1)
DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)

# quickly run clustering
aa <- FindNeighbors(aa, dims = 1:30)

aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 2.0)
p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

ggsave(filename = paste0(outDir, '/UMAP_subset_condition_clusters_res2.0.pdf'), 
       width = 14, height = 6)

clusters_sel = c(13, 18, 25, 30)
cells = colnames(aa)[which(!is.na(match(aa$seurat_clusters, clusters_sel)))]

jj = match(cells_2filter, colnames(aa))
xx = table(aa$seurat_clusters[jj])
yy = xx/table(aa$seurat_clusters)
xx[order(-xx)]
yy[order(-yy)]

xx = subset(aa, cells = cells)
saveRDS(xx, file = paste0(RdataDir, 'subObj_clusters_to_filter.rds'))

##########################################
# filter the clusters and explore UMAP parameters 
##########################################
aa =  readRDS(file = paste0('../results/Rdata/', 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', 
                            species, '_R13547_10x_mNT_20220813', '.rds'))
Idents(aa) = factor(aa$condition, levels = levels)
table(aa$condition)

Filter_weirdCluster = FALSE
if(Filter_weirdCluster){
  
  cells_2filter = readRDS(file = paste0(RdataDir, 'subObj_clusters_to_filter.rds'))
  mm = match(colnames(aa), colnames(cells_2filter))
  
  aa = subset(aa, cells = colnames(aa)[which(is.na(mm))])
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.1)
  DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)
  
  saveRDS(aa, file = paste0(RdataDir, 'seuratObj_clustersFiltered_umapOverview.rds'))
  
}

## explore the umap parameters
aa = readRDS(file = paste0(RdataDir, 'seuratObj_clustersFiltered_umapOverview.rds'))

source(paste0(functionDir, '/functions_scRNAseq.R'))
explore.umap.params.combination(sub.obj = aa, resDir = outDir, 
                                pdfname = 'UMAP_test_allSamples_selectedParameters.pdf',
                                use.parallelization = FALSE,
                                group.by = 'condition',
                                cols = cols, 
                                weight.by.var = TRUE,
                                nfeatures.sampling = c(5000),
                                nb.pcs.sampling = c(30), 
                                n.neighbors.sampling = c(100),
                                min.dist.sampling = c(0.1)
)


DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE, cols = cols)

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
## test Hannah's umap parameter: nfeatures = 5000;nb.pcs=30;n.neighbors = 100; min.dist = 0.1
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.1, spread = 1)

DimPlot(aa, label = FALSE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)

ggsave(filename = paste0(resDir, '/scRNAseq_overview_allSamples.noLabels', version.analysis, '.pdf'), 
       width = 10, height = 8)

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObj_clustersFiltered_umapOverview_selectedUmapParams.rds'))


Idents(aa) = factor(aa$condition, levels = levels)
aa$condition = factor(aa$condition, levels = levels)

##########################################
# plot some features
##########################################
outDir = paste0(resDir, '/scRNAseq_featurePlots/')
if(!dir.exists(outDir)) dir.create(outDir)

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_clustersFiltered_umapOverview_selectedUmapParams.rds'))

Idents(aa) = aa$condition

DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)


features = c('Pou5f1', 'Sox2', 'Lef1', 'Otx2', 'Zfp703', 'Pax6', 'Foxa2', 'Shh', 'Nkx6-1', 'Nkx2-2', 'Olig2', 
             'Sox1', 'Tubb3', 'Bmp4', 'Bmp7', 'Nog', 'Pax3', 'Pax7', 'Arx')
features = c(features, c('Pou5f1', 'Sox2', 'Nanog', 'Zfp42', 'Klf4', 'Esrrb', 'Nr0b1', 'Dazl', 
             'Pou3f1', 'Otx2', 'Klf2', 'Klf5', 'Etv4', 'Etv5'))

features = c(features, c('Foxa2', 'Shh', 'Arx', 'Pax6', 'Sox1', 'Nkx2-2', 'Nkx6-1', 'Olig2', 'Tubb3', 
             'Map2', 'Rbfox3', 'Irx3', 'Irx5', 'Sp8', 'Dbx2', 'Msx2', 'Lmx1b', 'Pax3', 'Pax7', 
             'Sox2', 'Elavl3', 'Sox10', 'Tlx2', 'Six1', 'Bmp4', 'Bmp7', 'Wnt1', 'Olig3'))

features = c(features, c('Foxa2', 'Shh', 'Arx', 'Ferd3l', 'Lmx1b', 'Sox2', 'Nkx6-1', 
             'Tubb3', 'Elavl3')) # floor plate markers

features = c(features, c('Foxa2', 'Shh', 'Arx', 'Pax6', 'Sox1', 'Nkx2-2', 'Nkx6-1', 'Olig2', 'Tubb3', 
             'Map2', 'Rbfox3', 'Irx3', 'Irx5', 'Sp8', 'Dbx2', 'Msx2', 'Lmx1b', 'Pax3', 'Pax7', 
             'Sox2', 'Elavl3', 'Sox10', 'Tlx2', 'Six1', 'Bmp4', 'Bmp7', 'Wnt1', 'Olig3'))

features = unique(c(features, c('Sox2', 'Sox1', 'Tubb3', 'Elavl3', 
                    'Irx3', 'Irx5', 'Pax3', 'Pax7',
                    'Pax6', 'Olig2', 'Nkx2-9', 'Nkx2-2', 
                    'Nkx6-1', 'Foxa2', 'Arx', 'Shh'
))) # DV overview

for(n in 1:length(features))
{
  # n = 1
  cat(n, ' -- gene ', features[n], '\n')
  FeaturePlot(aa, features = features[n]) + 
    scale_color_viridis_c() +
    theme(axis.text.x = element_text(angle = 0, size = 14), 
          axis.text.y = element_text(angle = 0, size = 14), 
          axis.title =  element_text(size = 14),
          legend.text = element_text(size=12),
          legend.title = element_text(size = 14)
          #legend.position=c(0.2, 0.8),
          #plot.margin = margin()
          #legend.key.size = unit(1, 'cm')
          #legend.key.width= unit(1, 'cm')
    )
  
  ggsave(paste0(outDir, '/FeaturePlot_allSamples_', features[n], '.pdf'),  width=8, height = 6)
  
  
}

DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE) +
  theme(axis.text.x = element_text(angle = 0, size = 14), 
        axis.text.y = element_text(angle = 0, size = 14), 
        axis.title =  element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)
        #legend.position=c(0.2, 0.8),
        #plot.margin = margin()
        #legend.key.size = unit(1, 'cm')
        #legend.key.width= unit(1, 'cm')
  )

ggsave(paste0("../results/plots_MondaySeminar", 
              '/UMAP_condition_toUse.pdf'),  width=8, height = 6) 


########################################################
########################################################
# Section : RA samples
# 
########################################################
########################################################
outDir = paste0(resDir, '/UMAP_allSamples/')
if(!dir.exists(outDir)) dir.create(outDir)

aa =  readRDS(file = paste0('../results/Rdata/',
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_',
                            species, '_R13547_10x_mNT_20220813', '.rds'))

Idents(aa) = factor(aa$condition, levels = levels)
table(aa$condition)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)


Filter_weirdCluster_7_8_old = FALSE
if(Filter_weirdCluster_7_8){
  
  cells_2filter = readRDS(paste0(RdataDir, 'cellNames_weirdClusters_7_8.rds'))
  mm = match(colnames(aa), cells_2filter)
  aa = subset(aa, cells = colnames(aa)[which(is.na(mm))])
  
}

Filter_weirdCluster_new = FALSE
if(Filter_weirdCluster){
  
  cells_2filter = readRDS(file = paste0(RdataDir, 'subObj_clusters_to_filter.rds'))
  mm = match(colnames(aa), colnames(cells_2filter))
  
  aa = subset(aa, cells = colnames(aa)[which(is.na(mm))])
  
  # aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
  # ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  # aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
  # ElbowPlot(aa, ndims = 50)
  # 
  # Idents(aa) = aa$condition
  # 
  # aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.1)
  # DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)
  #saveRDS(aa, file = paste0(RdataDir, 'seuratObj_clustersFiltered_umapOverview.rds'))
  
}

##########################################
# redo umap for RA d2.5-d6 
##########################################
#bb = readRDS(file = paste0(RdataDir, 
#                          'seuratObj_clustersFiltered_umapOverview_selectedUmapParams.rds'))

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_clustersFiltered_umapOverview_selectedUmapParams.rds'))

names(cols) = levels
levels_sels = c("day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
cols_sel = cols[match(levels_sels, names(cols))]

Idents(aa) = factor(aa$condition, levels = levels)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)

aa = subset(aa, idents = levels_sels)

source(paste0(functionDir, '/functions_scRNAseq.R'))

outDir = paste0(resDir, '/UMAP_RAsamples/')
if(!dir.exists(outDir)) dir.create(outDir)

explore.umap.params.combination(sub.obj = aa, resDir = outDir, 
                                pdfname = 'UMAP_test_RASamples_selectedParameters_seurat5_v4.pdf',
                                use.parallelization = FALSE,
                                group.by = 'condition',
                                cols = cols, 
                                weight.by.var = TRUE,
                                nfeatures.sampling = c(3000),
                                nb.pcs.sampling = c(50), 
                                n.neighbors.sampling = c(30, 50), 
                                min.dist.sampling = c(0.1)
)


saveRDS(aa, file = paste0(RdataDir, 'seuratObj_clustersFiltered_umap_RAsamples_4save.rds'))

##########################################
# UMAP visualization of the RA sample 
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObj_clustersFiltered_umap_RAsamples_4save.rds'))

levels_sels = c("day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
cols_sel = cols[match(levels_sels, names(cols))]

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE,
             npcs = 50, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 50, min.dist = 0.1, spread = 1)

DimPlot(aa, label = FALSE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)  

ggsave(filename = paste0(resDir, '/scRNAseq_overviewUMAP_RAsamples_seurat5.0.2_final_v1.pdf'), 
       width = 10, height = 8)

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam_s5.rds'))

##########################################
# Fig 1C clustering of RA samples
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam_s5.rds'))

### test some clustering options
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.8)
DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)

DimPlot(aa, group.by = "Phase", label = TRUE, repel = TRUE, raster=FALSE)

knownGenes =  c('Pou5f1', 'Sox2', 'Zfp42', 'Utf1',
                'Otx2', 'Cyp26a1', 'Stra8',
                'Hoxa1', 'Hoxa3', 'Hoxb4',
                'Sox1', 'Pax6', 'Tubb3', 'Elavl4', 'Neurod4',
                'Foxa2', 'Shh', 'Arx', 'Vtn', "Spon1", 'Slit2', "Ntn1",
                'Nkx2-2', 'Olig2', 'Pax3', 'Pax7', 'Nkx6-1')

Discard_cellCycle.corrrelatedGenes = TRUE
if(Discard_cellCycle.corrrelatedGenes){
  library(scater)
  Idents(aa) = aa$condition
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  
  # Identifying the likely cell cycle genes between phases,
  # using an arbitrary threshold of 5%.
  scaledMatrix = GetAssayData(aa, slot = c("scale.data"))
  
  diff <- getVarianceExplained(scaledMatrix, data.frame(phase = aa$Phase))
  diff = data.frame(diff, gene = rownames(diff))
  diff = diff[order(-diff$phase), ]
  
  hist(diff$phase, breaks = 100); 
  abline(v = c(1:5), col = 'red')
  
  rm(scaledMatrix)
  
  genes_discard = diff$gene[which(diff$phase>5)]
  cat(length(genes_discard), 'genes to discard \n')
  
  tfs_sels = intersect(genes_discard, gene_examples)
  print(tfs_sels)
  
  knownGenes_sels = intersect(genes_discard, knownGenes)
  print(knownGenes_sels)
  
  if(length(knownGenes_sels)>0) genes_discard = setdiff(genes_discard, knownGenes_sels)
  
  tfs_sels = intersect(genes_discard, gene_examples)
  print(tfs_sels)
  
  aa = subset(aa, features = setdiff(rownames(aa), genes_discard))
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000) # find subset-specific HVGs
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- FindNeighbors(aa, dims = 1:30)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)
  
  aa$clusters = aa$seurat_clusters
  aa$clusters[which(aa$seurat_clusters == 3)] = 1
  aa$clusters[which(aa$seurat_clusters == 0)] = 2
  aa$clusters[which(aa$seurat_clusters == 1)] = 3
  aa$clusters[which(aa$seurat_clusters == 6)] = 4
  aa$clusters[which(aa$seurat_clusters == 2)] = 5
  aa$clusters[which(aa$seurat_clusters == 4)] = 6
  aa$clusters[which(aa$seurat_clusters == 5 | aa$clusters == 8)] = 7
  
}

bb = readRDS(file = paste0(RdataDir, 
                          'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam',
                          '_clustered.discardCellcycle.corrrelatedGenes.rds'))

aa$clusters = bb$clusters[match(colnames(aa), colnames(bb))]

cols_cluster = c( "#7F7F7F", "#FFC000", "#70AD47", "#337f01", "#CD00CF", "#265401","#800080")
# = c( "#7f7f7f", "#ffc000", "#70ad47", "#337f01", "#cd00cf", "#265401","#800080")

DimPlot(aa, label = FALSE, group.by =  'clusters', repel = TRUE, raster=FALSE, 
        cols = cols_cluster) 

ggsave(filename = paste0(resDir, '/scRNAseq_overview_RAsamples_clustering_final_v1.pdf'), 
       width = 10, height = 8)

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam_clusters_4save.rds'))


##########################################
# timpe point contribution of clusters 
##########################################
pcts = table(aa$clusters, aa$condition)
pcts = pcts[!is.na(match(rownames(pcts), c(1:7))), ]
for(n in 1:nrow(pcts)) pcts[n, ] = pcts[n, ]/sum(pcts[n, ])
pcts = as.matrix(pcts)

clusters = rep(rownames(pcts), each = 6)
samples = rep(colnames(pcts), nrow(pcts))
freqs = as.numeric(t(pcts))

pcts = data.frame(clusters, samples, freqs)

library(ggplot2)
levels_sels = c("day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
cols_sel = cols[match(levels_sels, names(cols))]


ggplot(pcts, aes(x = factor(clusters, levels = c(7:1)), y = freqs, fill = samples)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  scale_fill_manual(values=cols_sel) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  labs( x = 'Clusters', y = '% of cells from each time points' )

ggsave(filename = 
         paste0(resDir, '/scRNAseq_overview_RAsamples_clusters_timePointsContribution_v1.pdf'), 
       width = 8, height = 8)


##########################################
# highlight the co-expression of FoxA2 and Pax6 
##########################################
feature_blend = FeaturePlot(aa, features = c('Foxa2', 'Pax6'),
                            blend = TRUE, blend.threshold = 0.1, alpha = 0.4, order = TRUE, 
                            cols = c("#F0F0F0", "green", "magenta"))


blend_only = feature_blend[[3]] + theme_void() + theme(legend.position = "none")
blend_legend = feature_blend[[4]] 
#rel_widths = c(0.7, 0.3), nrow = 1)

ggsave(filename = file.path(resDir, "FeaturePlot_RAonly_FoxA2-Pax6-blend_LEGEND.pdf"), 
       plot = blend_legend, width = 4, height = 4, units = "in")

ggsave(filename = file.path(resDir, "FeaturePlot_RAonly_FoxA2-Pax6-blend_UMAP.pdf"), 
       plot = blend_only, width = 8, height = 6, units = "in")

#Cairo::CairoPNG(filename = file.path(resDir, "FeaturePlot_RAonly_FoxA2-Pax6-blend_UMAP.png"), 
#                width = 8, height = 6, units = "in", dpi = 300)
#print(blend_only)
#dev.off()


##########################################
# highlight the marker genes or specified genes
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam_clusters_4save.rds'))

Run_allMarkerGene_searching = FALSE
if(Run_allMarkerGene_searching){
  Idents(aa) = factor(aa$clusters)
  allMarker <- FindAllMarkers(aa, only.pos = TRUE)
  
  saveRDS(allMarker, file = paste0(RdataDir, 
                                   'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam',
                                   '_clustered.discardCellcycle.corrrelatedGenes_markerGenes.rds'))
  
}else{
  allMarker = readRDS(file = paste0(RdataDir, 
                                    'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam',
                                    '_clustered.discardCellcycle.corrrelatedGenes_markerGenes.rds'))
  
}

Idents(aa) = factor(aa$condition)
xx = subset(x = aa, downsample = 2000)

allMarker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top10

write.csv(allMarker, file = paste0(resDir, '/Table_S1_clusters_markerGenes_all.csv'),
          row.names = FALSE)

write.csv(top10, file = paste0(resDir, '/Table_S1_clusters_markerGenes_top20_FigS1D.csv'),
          row.names = FALSE)

xx$clusters = droplevels(xx$clusters)

xx$clusters = factor(xx$clusters, levels = c(1, 2, 3, 4, 6, 5, 7))

cols_cluster = c( "#7F7F7F", "#FFC000", "#70AD47", "#337f01", "#265401", "#CD00CF", "#800080")
#cols_cluster = c("#F0F0F0", "#EFFAB6", "#69C6BE", "#007BB7", "#121D60")
DoHeatmap(xx, group.by = 'clusters', features = top10$gene, draw.lines = TRUE, disp.min = -2.,
          group.colors = cols_cluster, angle = 0, size = 5) + 
  theme(text = element_text(size = 6)) +
  #NoLegend() +
  #scale_fill_gradientn(colors = feat_cols)
  scale_fill_viridis_c(option = "magma")
  #scale_fill_viridis_c() 
  #scico::scale_fill_scico(palette = "vik")
  #scale_fill_viridis(option = "D") 
ggsave(filename = paste0(resDir, '/scRNAseq_RAsamples_clustering_allMarkers_v2.pdf'), 
       width = 10, height = 6)

# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
#DotPlot(pbmc3k.final, features = features, split.by = "groups") + RotatedAxis()
require(data.table)
## show the marker genes https://ouyanglab.com/singlecell/clust.html
Idents(aa) = factor(aa$clusters)
oupMarker <- FindAllMarkers(aa, features = knownGenes, min.pct = 0.01, logfc.threshold = 0.01)
oupMarker <- data.table(oupMarker)
oupMarker$pct.diff = oupMarker$pct.1 - oupMarker$pct.2
oupMarker <- oupMarker[, c("cluster","gene","avg_log2FC","pct.1","pct.2",
                           "pct.diff","p_val","p_val_adj")]
#fwrite(oupMarker, sep = "\t", file = "images/clustMarkers.txt")
aa@misc$marker <- oupMarker      # Store markers into Seurat object

# Get top genes for each cluster and do dotplot / violin plot
#oupMarker$cluster = factor(oupMarker$cluster, levels = reorderCluster)
oupMarker = oupMarker[order(cluster, -avg_log2FC)]

knownGenes =  c('Zfp42', 'Utf1', 'Pou5f1',  
                 'Cyp26a1', 'Stra8','Hoxa1', 'Hoxa3', 'Hoxb4', 'Sox2',
                'Pax6', 'Sox1', 'Pax3',   'Irx3', 'Irx5',
                'Foxa2', 'Shh', 'Arx', "Spon1", 'Slit2', "Ntn1",
                'Nkx6-1', 'Nkx2-2',  'Olig2',  
                 'Tubb3',  'Neurod4', 'Elavl4', 'Pou3f2')
genes.to.plot <- knownGenes

library(RColorBrewer)
feat_cols = c("#F0F0F0", "#EFFAB6", "#69C6BE", "#007BB7", "#121D60")
#to set the gradient of the colour scale from min to max
feat_values = c(-1, -0.5, 0.5, 1.5)
#colGEX = c("grey85", brewer.pal(7, "Reds"))
aa$clusters = factor(aa$clusters, levels = c(7, 5, 6, 4, 3, 2, 1))

DotPlot(aa, group.by = "clusters", features = genes.to.plot) + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  #scale_colour_viridis(option="magma") +
  #guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  scale_color_gradientn(colors = feat_cols) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(axis.text.x = element_text(angle = 90, size = 12, vjust = 0.),
        axis.text.y = element_text(angle = 0, size = 14)) +
  labs( x = '', y = 'Clusters' )

ggsave(filename = paste0(resDir, '/scRNAseq_overview_RAsamples_clustering_markers_dotplot_sorted_v2.pdf'), 
       width = 8, height = 5)


nClust = 7
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(nClust)
VlnPlot(aa, group.by = "clusters", fill.by = "ident", cols = cols_cluster,
              features = genes.to.plot, stack = TRUE, flip = TRUE)

ggsave(filename = paste0(resDir, '/scRNAseq_overview_RAsamples_clustering_markers_Vlnplot.pdf'), 
       width = 8, height = 6)


##########################################
# highlight the DE genes in cluster 3 and 5 to characterize the cell states in salt-and-pepper  
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam',
                           '_clustered.discardCellcycle.corrrelatedGenes.rds'))

DimPlot(aa, group.by = 'clusters', label = TRUE)

Idents(aa) = factor(aa$clusters)

cluster5.markers <- FindMarkers(aa, ident.1 = 5, ident.2 = c(3), 
                                logfc.threshold = 0.2, only.pos = TRUE)

cluster3.markers <- FindMarkers(aa, ident.1 = 3, ident.2 = c(5), logfc.threshold = 0.2, only.pos = TRUE)

head(cluster5.markers, n = 10)
head(cluster3.markers, n = 10)

marker_clusters_3_5 = c(rownames(cluster5.markers)[1:20], rownames(cluster3.markers)[1:20])

library(RColorBrewer)
colGEX = c("grey85", brewer.pal(7, "Reds"))
DotPlot(aa, group.by = "clusters", idents = c(3, 5),
        features = marker_clusters_3_5) + 
  coord_flip() + scale_color_gradientn(colors = colGEX) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0))

ggsave(filename = paste0(resDir, '/scRNAseq_overview_RAsamples_clustering_markers_dotplot_sorted.pdf'), 
       width = 8, height = 4)

#nClust = 7
#colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(nClust)
VlnPlot(aa, group.by = "clusters", fill.by = "ident", 
        idents = c(3, 5), 
        cols = cols_cluster[c(3, 5)],
        features = marker_clusters_3_5, stack = TRUE, flip = TRUE)

ggsave(filename = paste0(resDir, '/scRNAseq_RAsamples_clusterMarkers_3vs5_Vlnplot.pdf'), 
       width = 4, height = 8)

##########################################
# double positive cell in RA samples with scatter plots  
##########################################
p1 = FeatureScatter(aa, feature1 = "Pax6", feature2 = "Foxa2", group.by = "condition", cols = cols_sel,
                    smooth = FALSE, shuffle = TRUE)

p2 = FeatureScatter(aa, feature1 = "Pax6", feature2 = "Foxa2", group.by = "clusters", cols = cols_cluster,
                    smooth = FALSE)

## see help from https://github.com/satijalab/seurat/issues/3544
data <- FetchData(aa, vars = c('Pax6', 'Foxa2', 'Sox1'))
# Use position = 'jitter' to apply a slight jitter to the points, makes it look nicer
# c('lightgrey', 'blue') are the standard FeaturePlot colors
# cowplot::theme_cowplot() is the the standard theme used in Seurat
p3 = ggplot(data = data) +
  geom_point(mapping = aes(x = Pax6, y = Foxa2, color = Sox1), position = 'jitter') +
  scale_color_gradientn(
    colors = c('lightgrey', "maroon4")
    #colors = viridis::magma(9)
    ) +
  cowplot::theme_cowplot()

p1 + p2 + p3

ggsave(filename = paste0(resDir, '/scRNAseq_RAsamples_FoxA2_Pax6_coexpression_shuffled.pdf'), 
       width = 20, height = 6)

### counting the double positive cells 
cc = unique(aa$condition)
counts = rep(NA, length(cc))
names(counts) = cc

kk = match(c('Foxa2', 'Pax6'), rownames(aa))
for(n in 1:length(counts))
{
  jj = which(aa$condition == cc[n])
  counts[n] = length(which(aa@assays$RNA@data[kk[1], jj] >0 & aa@assays$RNA@data[kk[2], jj] >0))/length(jj)
}

df = data.frame(condition = names(counts), counts)

ggplot(data=df, aes(x=condition, y=counts, fill=condition)) +
  geom_bar(stat="identity")+
  #geom_text(aes(label=condition), vjust=1.6, color="white", size=3.5)+
  scale_fill_manual(values=cols_sel) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 10, vjust = 0.4),
        axis.text.y = element_text(angle = 45, size = 10)) +
  labs( x = '', y = '% of FoxA2+ & Pax6+' )

ggsave(filename = paste0(resDir, '/scRNAseq_RAsamples_FoxA2_Pax6_percentage.pdf'), 
       width = 8, height = 5)


########################################################
########################################################
# Section : FACS data diffusion map for RA samples
# 
########################################################
########################################################
library(CATALYST)
library(SingleCellExperiment)

RdataDir_facsData = paste0('../results/FACS_analysis_clusteringWT', '/Rdata')

sce = readRDS(file = paste0(RdataDir_facsData, '/sce_FACS_RA_nod2_DiffusionMap_saved.rds'))

sce$clusters = factor(sce$clusters)


p0 = plotDR(sce, "DiffusionMap", color_by = "condition") +
  theme_classic() +
  scale_colour_manual(values = c('#fee8c8', '#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f',
                                 '#b30000', '#7f0000')) +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12))

p1 = plotDR(sce, "DiffusionMap", color_by = "clusters") +
  theme_classic() + 
  scale_colour_manual(values = c("#464646", "#7F7F7F", "#CD00CF", "#FFC000", "#70AD47"))+
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12))

p0 + p1 

ggsave(paste0(resDir, '/FACS_RA_DiffusionMap_timePoints_clusters_colors.pdf'), width=14, height = 6)


p2 = plotDR(sce, "DiffusionMap", color_by = "Oct4") + theme_classic() + scale_color_gradientn(colors = feat_cols)
p3 = plotDR(sce, 'DiffusionMap', color_by = 'Pax6') + theme_classic() + scale_color_gradientn(colors = feat_cols)
p4 = plotDR(sce, 'DiffusionMap', color_by = 'FoxA2') + theme_classic() + scale_color_gradientn(colors = feat_cols)
p5 = plotDR(sce, 'DiffusionMap', color_by = 'Sox1') + theme_classic() + scale_color_gradientn(colors = feat_cols)
p6 = plotDR(sce, 'DiffusionMap', color_by = 'Sox2') + theme_classic() + scale_color_gradientn(colors = feat_cols)

p2 + p3 + p4 + p5 + p6 

ggsave(paste0(resDir, '/FACS_RA_DiffusionMap_geneExpression_v2.pdf'), width=16, height = 8)


########################################################
########################################################
# Section : plot the feature selection outcome 
# 
########################################################
########################################################
Merge_sparseFeatures_mutlipleMethods = FALSE
if(Merge_sparseFeatures_mutlipleMethods){
  dataDir = paste0("../results/scRNAseq_R13547_10x_mNT_20220813/",
                   "RA_symetryBreaking/sparse_featureSelection_d2_d2.5_d3_d3.5_d4_d5/")
  
  aa = readRDS(paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/',
                      'RA_competence/scRNAseq_mNT/saved_seuratObj/',
                      'RAsamples_d2_d6_UMAP_selectedParam.rds'))
  
  
  ggs = c()
  xx = read.csv(file = paste0(dataDir, 'geneBasis_top50_tfs_sps_final.csv'), sep = ";")
  ggs = unique(c(ggs, xx[, 2]))
  
  #write.csv(xx, file = paste0(resDir, '/geneBasis_top50_tfs_sps_final.csv'), row.names = FALSE)
  
  xx = read.csv(file = paste0(dataDir, 'dubstep_sparse_74.tfs.sps_final.csv'))
  ggs = unique(c(ggs, xx[, 1]))
  
  #write.csv(xx, file = paste0(resDir, '/dubstep_sparse_74.tfs.sps_final.csv'), row.names = FALSE)
  
  xx = c('Mdk', 'Sox11', 'Cdh1', 'Shh', 'Pou3f2', 'Rfx4', 'Cyp26a1', 'Apoe', 'Nedd4', 'Meis2', 'Cdh2', 'Foxa1',
         'Tshz1', 'Peg10', 'Pax6', 'Pou3f3', 'Tfap2c', 'Cebpb', 'Dll1', 'Gpc3', 'Peg3', 'Vcp', 'Ddx3x', 'Zic1',
         'Wnt6', 'Pfdn5', 'Rack1', 'Foxa2', 'Nkx6-2', 'Prrx2', 'Plagl1')
  
  match(xx, rownames(aa))
  
  ggs = unique(c(ggs, xx))
    
  ggs = ggs[which(!is.na(match(ggs, rownames(aa))))]
  saveRDS(ggs, file = paste0(resDir, '/geneList_sparseFeatureSelection_v2.rds'))
  
  ### save the sparse features for Table_S1
  ggs = readRDS(file = paste0(resDir, '/geneList_sparseFeatureSelection_v2.rds'))
  
  ggs = data.frame(features = ggs, geneBasis = rep(NA, length(ggs)), 
                   dubstep = rep(NA, length(ggs)), SMD = rep(NA, length(ggs)))
  
  ggs$SMD[!is.na(match(ggs$features, xx))] = 'yes' 
  
  xx = read.csv(file = paste0(dataDir, 'geneBasis_top50_tfs_sps_final.csv'), sep = ";")
  
  ggs$geneBasis[!is.na(match(ggs$features, xx$gene))] = 'yes'
  
  xx = read.csv(file = paste0(dataDir, 'dubstep_sparse_74.tfs.sps_final.csv'))
  ggs$dubstep[!is.na(match(ggs$features, xx[,1]))] = 'yes'
  
  rm(xx)
  
  nb_pre = c()
  for(n in 1:nrow(ggs)) nb_pre = c(nb_pre, length(which(ggs[n, c(2:4)] == 'yes')))
  
  o1 = order(-nb_pre)
  
  ggs = ggs[o1, ]
  
  write.csv(ggs, file = paste0(resDir, '/Table_S1_sparseFeatures.csv'),
            row.names = FALSE)
  
  
  #### umap of samples only using selected sparse features
  ggs = readRDS(file = paste0(resDir, '/geneList_sparseFeatureSelection_v2.rds'))
  
  levels_sels = c("day2_beforeRA",  
                  "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA")
  
  #levels_sels = unique(aa$condition)
  cols_sel = cols[match(levels_sels, names(cols))]
  Idents(aa) = as.factor(aa$condition)
  
  aa = subset(aa, idents = levels_sels)
  
  aa$condition = droplevels(aa$condition)
  
  aa <- RunPCA(aa, features = ggs, verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  
  outDir = paste0(resDir, '/UMAP_RAsamples_sparseFeatures/')
  if(!dir.exists(outDir)) dir.create(outDir)
  
  source(paste0(functionDir, '/functions_scRNAseq.R'))
  explore.umap.params.combination(sub.obj = aa, 
                                  resDir = outDir, 
                                  pdfname = 'UMAP_test_RASamples_sparseFeatures_weight.by.var.FALSE.pdf',
                                  use.parallelization = FALSE,
                                  group.by = 'condition',
                                  cols = cols, 
                                  weight.by.var = FALSE,
                                  #nfeatures.sampling = c(3000, 5000),
                                  features = ggs,
                                  nb.pcs.sampling = c(10, 20), 
                                  n.neighbors.sampling = c(20, 30, 50), 
                                  min.dist.sampling = c(0.01, 0.05, 0.1)
  )
  
  
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 30, min.dist = 0.01)
  DimPlot(aa, group.by = 'condition', label = TRUE, cols = cols_sel)
  
  ggsave(filename = paste0(resDir, '/RA_umap_with_selectedFeatures_v2.pdf'), width = 8, height = 6)
  
  
}

##########################################
# Visualising RA-specific expression pattern / kinetics (feature plots / exp vs pseudotime?  
# by Hannah
##########################################
aa = readRDS(paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/',
                    'scRNAseq_mNT/saved_seuratObj/',
                    'RA_noRAsamples_d2_d6_UMAP_selectedParam.rds'))

DimPlot(aa, group.by = 'condition', label = TRUE)

levels_sels = c("day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA", "day6_RA", 
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA", "day6_noRA")

Idents(aa) = as.factor(aa$condition)
aa = subset(aa, idents = levels_sels)

levels_sels = unique(aa$condition)
cols_sel = cols[match(levels_sels, names(cols))]

DimPlot(aa, cols = cols_sel)

aa$groups = NA
aa$groups[grep('_RA', aa$condition)] = 'RA'
aa$groups[grep('_noRA', aa$condition)] = 'noRA'

aa$condition = droplevels(aa$condition)

aa$time = sapply(aa$condition, function(x) {unlist(strsplit(as.character(x), '_'))[1]})

features =c('Pou5f1', 'Foxa2', 'Pax6',  'Sox1')

aa$groups = factor(aa$groups, levels = c('noRA', 'RA'))

VlnPlot(aa, features = features, split.by = "groups", group.by = 'time', log = FALSE, 
        same.y.lims = TRUE, ncol = 1, pt.size = 0, 
        cols = c("#6BAED6", "#EF6548"))

ggsave(filename = paste0(resDir, '/RAspecific_geneExpression_kinetics_noDot_v3.pdf'), 
       width = 6, height = 14) 


##########################################
# Dotplot for RA vs noRA
##########################################
# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
aa = readRDS(paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/',
                    'scRNAseq_mNT/saved_seuratObj/',
                    'RA_noRAsamples_d2_d6_UMAP_selectedParam.rds'))

ggs = readRDS(file = paste0(resDir, '/geneList_sparseFeatureSelection_v2.rds'))

levels_sels = c("day2_beforeRA",
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA", "day6_RA", 
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA", "day6_noRA")

Idents(aa) = as.factor(aa$condition)
aa = subset(aa, idents = levels_sels)

levels_sels = unique(aa$condition)
cols_sel = cols[match(levels_sels, names(cols))]

DimPlot(aa, cols = cols_sel)

aa$condition = droplevels(aa$condition)

levels_reorder = c("day6_noRA", "day5_noRA", "day4_noRA", 
                   "day3.5_noRA", "day3_noRA",  "day2.5_noRA",
                   "day2_beforeRA", 
                   "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   
                   "day4_RA", "day5_RA", "day6_RA")

aa$condition = factor(as.character(aa$condition), 
                      levels = levels_reorder
)


clustering_DotPlot = FALSE
if(clustering_DotPlot){
  ## this pacakge is installed in R_4.3.2 and Seurat_5.0.2
  # can't be install in lower versions
  require(scCustomize)  
  aa$orig.ident = factor(aa$condition, levels = levels_reorder)
  
  Idents(aa) = aa$condition
  #DotPlot_scCustom(aa, features = ggs[c(1:20)], split.by = "groups",  flip_axes = TRUE)
  
  pdf(paste0(resDir, '/RA_vs_noRA_sparseFeatures_clusteredDotPlot.pdf'), 
      height = 24, width =10, useDingbats = FALSE)
  Clustered_DotPlot(aa, features = ggs,
                    colors_use_exp = viridis::plasma(n = 20, direction = 1),
                    plot_km_elbow = FALSE,
                    cluster_ident = FALSE,
                    column_label_size = 12,
                    row_label_size = 10,
                    colors_use_idents = cols[match(levels_reorder, names(cols))],
                    column_names_side = "top",
                    #legend_position = "bottom", #legend_orientation = "horizontal", 
                    show_ident_legend = TRUE)
  
  dev.off()
  
  #ggsave(filename = paste0(resDir, '/RA_vs_noRA_sparseFeatures_clusteredDotPlot.pdf'), 
  #       width = 14, height = 20)
  
}else{
  
  #Idents(aa) = factor(aa$time)
  gene_sels = c('Zfp42', 'Pou5f1', 'Utf1', 'Pou3f1', 'Rarg', 'Tcf15', 'Otx2', 'Sox2', 'Pax6', 
                'Cdh1', 'Cdh2', 'Fbn2', 'Spry4', 'Fgf4', 'Etv5', 'Id1', 'Cyp26a1', 'Hoxa1', 'Zfp703', 
                'Tshz1', 'Msx1', 'Sox17', 'Foxa2', 'Foxa1', 'Nkx6-1', 'Shh')
  gene_sels = gene_sels[!is.na(match(gene_sels, ggs))]
  gene_sels = gene_sels[c(length(gene_sels):1)]
  
  feat_cols = c("#F0F0F0", "#EFFAB6", "#69C6BE", "#007BB7", "#121D60")
  
  DotPlot(aa, group.by = "condition", features = gene_sels) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
    scale_color_gradientn(colors = feat_cols) +
    guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
    RotatedAxis() + 
    coord_flip()  +
    theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
    theme(axis.text.x = element_text(angle = 90, size = 12, vjust = -0.4),
          axis.text.y = element_text(angle = 0, size = 12)) +
    labs( x = 'Sparse features', y = '' )
  
  ggsave(filename = paste0(resDir, '/RA_vs_noRA_sparseFeatures_v2.pdf'), width = 8, height = 10)
  
  
}

########################################################
########################################################
# Section : Quantify the cell hetereogeneity with RA and without RA  
# pseudo time calculation with palantir
# the pseudo time used here is from palantir
########################################################
########################################################

##########################################
# calculate the pseudo-time with palantir  
##########################################
Use_pseudotime_palantir = TRUE
if(Use_pseudotime_palantir){
  #aa = readRDS(file = paste0(RdataDir, 
  #                           'seuratObj_clustersFiltered_umapOverview_selectedUmapParams.rds'))
  
  #Idents(aa) = aa$condition
  #DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)
  
  outDir = paste0('../results/scRNAseq_R13547_10x_mNT_20220813/RA_symetryBreaking/',
                  'dataIntegration_timePoints_4pseudotime/')
  
  aa = readRDS(file = paste0(outDir, 
                             'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                             'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d5.noNeurons.rds'))
  
  Idents(aa) = factor(aa$condition)
  
  levels_sels = c("day2_beforeRA",  
                  "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", "day3.5_RA",   "day4_RA", "day5_RA",
                  "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA")
  
  cols_sel = cols[match(levels_sels, names(cols))]
  
  #Idents(aa) = aa$condition
  table(aa$condition)
  
  DimPlot(aa, group.by = 'condition', label = TRUE, cols = cols_sel)
  
  ## add DM components from palantir
  dm = read.csv(paste0(outDir, 'DM_components_palantir.csv'), row.names = c(1))
  dm = as.matrix(dm)
  colnames(dm) = paste0('DM_', c(1:ncol(dm)))
  
  aa[["dm"]] <- CreateDimReducObject(embeddings = dm, key = "DM_", assay = DefaultAssay(aa))
  
  DimPlot(aa, reduction =  'dm', group.by = 'condition', label = TRUE, cols = cols_sel)
  
  #aa = subset(aa, idents = levels_sels)
  
  pst = read.csv(paste0(outDir, 'palantir_pseudotime_d2_d5.csv'), row.names = c(1))
  xx = read.csv(paste0(outDir, 'palantir_fate_probabilities_d2_d5.csv'), row.names = c(1))
  mm = match(rownames(pst), rownames(xx))
  pst = data.frame(pst, prob_noRA = xx[mm, 1], prob_RA = xx[mm, 2], stringsAsFactors = FALSE)
  xx = read.csv(paste0(outDir, 'geosketch_N10k_d2_d5.csv'), header = FALSE)
  pst$sketch = NA
  pst$sketch[xx[,1]] = 1
  
  aa = AddMetaData(aa, metadata =  pst)
  
  aa$pseudot_noRA = NA
  aa$pseudot_RA = NA 
  aa$pseudot = NA
  
  FeaturePlot(aa, features = c('prob_noRA', 'prob_RA', 'palantir_pseudotime'))
  
  hist(aa$palantir_pseudotime[which(aa$condition == 'day2_beforeRA'| 
                                      aa$condition == 'day2.5_noRA')], breaks = 100)
  
  abline(v = 0.05, col = 'darkred')
  
  jj = which(aa$condition == 'day2_beforeRA' & aa$palantir_pseudotime <0.05)
  tt = aa$palantir_pseudotime[jj]
  tt = scales::rescale(tt, to = c(0, max(tt)))
  aa$pseudot_noRA[jj] = tt
  aa$pseudot_RA[jj] = tt
  
  jj = intersect(grep('_RA', aa$condition), which(aa$prob_RA > 0.6 | aa$condition == 'day2.5_RA'))
  tt = aa$palantir_pseudotime[jj]
  aa$pseudot_RA[jj] = scales::rescale(tt, to = c(0.05, 1))
  
  jj = intersect(grep('_noRA', aa$condition), which(aa$prob_noRA > 0.6))
  tt = aa$palantir_pseudotime[jj]
  aa$pseudot_noRA[jj] = scales::rescale(tt, to = c(0.0, 1))
  
  cell_beforeRA = colnames(aa)[which(!is.na(aa$pseudot_noRA) & aa$condition == 'day2_beforeRA')]
  cell_noRA = colnames(aa)[which(!is.na(aa$pseudot_noRA) & grepl('_noRA', aa$condition))]
  cell_RA = colnames(aa)[which(!is.na(aa$pseudot_RA) & grepl('_RA', aa$condition))]
  
  kk = match(cell_beforeRA, colnames(aa)) # cell before RA shared by two trajectories
  cat(length(which(is.na(kk))), 'missing cells \n')
  aa$pseudot[kk] = aa$pseudot_noRA[kk]
  max_pst_shared = max(aa$pseudot_noRA[kk])
  
  kk1 = match(cell_noRA, colnames(aa)) # cell in noRA trajectory
  cat(length(which(is.na(kk1))), 'missing cells \n')
  tt = aa$pseudot_noRA[kk1]
  aa$pseudot[kk1] =  tt
  
  kk2 = match(cell_RA, colnames(aa)) # cell in RA trajectory
  cat(length(which(is.na(kk2))), 'missing cells \n')
  tt = aa$pseudot_RA[kk2]
  aa$pseudot[kk2] = tt
  
  aa = subset(aa, cells = colnames(aa)[which(!is.na(aa$pseudot))])
  
  VlnPlot(aa, features = 'pseudot', group.by = 'condition', cols = cols_sel, pt.size = 0.0)
  
  ggsave(filename = paste0(outDir, 'pseudotime_vs_realTime_palantir.pdf'), width = 10, height = 8)
  
  
}else{
  #aa = readRDS(file = paste0(RdataDir, 
  #                           'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
  #                           'cellCycleScoring_annot.v1_reduction.DM_princurves_pseudotime_',
  #                           species, '_R13547_10x_mNT_20220813', '.rds'))
  
  Idents(aa) = factor(aa$condition)
  
  levels_sels = c("day2_beforeRA",  
                  "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", "day3.5_RA",   "day4_RA", "day5_RA",
                  "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA")
  
  
  cols_sel = cols[match(levels_sels, names(cols))]
  
  Idents(aa) = aa$condition
  table(aa$condition)
  aa = subset(aa, idents = levels_sels)
  aa = subset(aa, cells = colnames(aa)[which(!is.na(aa$pseudot))])
  
}

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, group.by = 'condition', label = TRUE, cols = cols_sel)

p1 = DimPlot(aa, group.by = 'condition', label = TRUE, cols = cols_sel)
p2 = FeaturePlot(aa, features = 'pseudot', cols = c('lightgray', 'blue')) +
  ggtitle(label = 'pseudo time') +
  scale_color_viridis_c(direction = -1)

p1 + p2


DimPlot(aa, group.by = 'condition', label = TRUE, cols = cols_sel)

ggsave(filename = paste0(figureDir, 'RA_noRA_d2_d5_condition_pseudotime.pdf'), width = 16, height = 6) 

saveRDS(aa, file = paste0(RdataDir, 'RA_noRA_d2_d5_condition_pseudotimePalantir_saved4heterogeity.rds'))


### make plot for palantir pseudotime
levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA",
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA")

cols_sel = cols[match(levels_sels, names(cols))]


p1 = DimPlot(aa, group.by = 'condition', label = FALSE, cols = cols_sel) +
  scale_y_reverse() +
  scale_x_reverse()

p2 = FeaturePlot(aa, features = 'pseudot', cols = c('lightgray', 'blue')) +
  ggtitle(label = 'pseudo time') +
  scale_color_viridis_c(direction = -1) +
  scale_y_reverse() +
  scale_x_reverse()

p1 + p2

ggsave(filename = paste0(resDir, '/RA_noRA_d2_d5_condition_pseudotimePalantir_v3.pdf'), 
       width = 10, height = 6) 

VlnPlot(aa, features = 'pseudot', group.by = 'condition', pt.size = 0.00, cols = cols_sel) +
  geom_hline(yintercept = c(0.05), col = 'red') +
  ggtitle(label = 'pseudotime distribution per condition')

ggsave(filename = paste0(resDir, '/RA_noRA_d2_d5_condition_pseudotime_v2.pdf'), width = 10, height = 6) 

## plot the pseudo-time in the original UMAP
aa = readRDS(paste0(RdataDir, 'RA_noRAsamples_d2_d6_UMAP_selectedParam.rds'))

DimPlot(aa, group.by = 'condition')

xx = readRDS(file = paste0(RdataDir, 'RA_noRA_d2_d5_condition_pseudotimePalantir_saved4heterogeity.rds'))

FeaturePlot(xx, features = 'pseudot', cols = c('lightgray', 'blue')) +
  ggtitle(label = 'pseudo time') +
  scale_color_viridis_c(direction = -1)

aa$pseudot = NA
mm = match(colnames(aa), colnames(xx))
jj = which(!is.na(mm))
aa$pseudot[jj] = xx$pseudot[mm[jj]]

FeaturePlot(aa, cells = colnames(aa)[which(!is.na(aa$pseudot))], 
            features = 'pseudot', cols = c('lightgray', 'blue')) +
  ggtitle(label = 'pseudo time') +
  #scale_color_viridis_c(option = "magma")
  #scico::scale_color_scico(palette = "vik")
  scale_color_viridis_c(direction = -1)

ggsave(filename = paste0(resDir, '/Palantir_pseudotime_RAnonRA_umap.pdf'), 
       width = 8, height = 6) 


##########################################
###### test heterogeneity quantification by binning pseudotime
##########################################
aa = readRDS(file = paste0(RdataDir, 'RA_noRA_d2_d5_condition_pseudotimePalantir_saved4heterogeity.rds'))


## specify the bins manually
#jj = intersect(which(!is.na(aa$pseudot_RA)))

pdf(paste0(resDir, "/pseudotime_distribution_binning_RA_tinyBins.pdf"),
    height = 4, width =6, useDingbats = FALSE)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3, 4, 2, 1), tcl = -0.1)

hist(aa$pseudot_RA, breaks = 50, main = 'distribution of pseudotime by RA')
abline(v = c(seq(0.15, 0.25, by = 0.05), seq(0.5, 0.9, by = 0.05)), col = 'red')
#abline(v = seq(0.1, 0.9, length.out = 50), col = 'red')

dev.off()

pdf(paste0(resDir, "/pseudotime_distribution_binning_noRA.pdf"),
    height = 4, width =6, useDingbats = FALSE)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3, 4, 2, 1), tcl = -0.1)

hist(aa$pseudot_noRA, breaks = 50)
abline(v = c(seq(0.15, 0.45, by = 0.05), seq(0.7, 0.9, by = 0.05)), col = 'red')
#abline(v = seq(0.1, 0.9, length.out = 50), col = 'red')

dev.off()

aa$pst_group = NA
#bins = seq(0.1, 0.9, length.out = 16)

bin_interval = 0.05
bins = c(seq(0.15, 0.25, by = bin_interval), seq(0.5, 0.9, by = bin_interval))
for(n in 1:(length(bins)-1))
{
  if(bins[n+1]<0.3){
    jj1 = which(aa$pseudot_RA < bins[(n+1)] & aa$pseudot_RA >= bins[n])
    aa$pst_group[jj1] = paste0('RA_groupPst_', n)
    cat(n, '-- ', length(jj1), ' RA cells -- pst interval :', bins[n], ' --- ', bins[n+1],  '\n')
  }
  if(bins[n] > 0.4){
    jj1 = which(aa$pseudot_RA < bins[(n+1)] & aa$pseudot_RA >= bins[n])
    aa$pst_group[jj1] = paste0('RA_groupPst_', (n-1))
    cat(n, '-- ', length(jj1), ' RA cells -- pst interval :', bins[n], ' --- ', bins[n+1],  '\n')
  }
}

bins = c(seq(0.15, 0.45, by = bin_interval), seq(0.7, 0.9, by = bin_interval))

for(n in 1:(length(bins)-1))
{
  if(bins[n+1] < 0.5){
    jj1 =  which(aa$pseudot_noRA < bins[(n+1)] & aa$pseudot_noRA >= bins[n])
    aa$pst_group[jj1] = paste0('noRA_groupPst_', n)
    cat(n, '-- ', length(jj1), ' noRA cells -- pst interval :', bins[n], ' --- ', bins[n+1],  '\n')
  }
  if(bins[n] > 0.6){
    jj1 = which(aa$pseudot_noRA < bins[(n+1)] & aa$pseudot_noRA >= bins[n])
    aa$pst_group[jj1] = paste0('noRA_groupPst_', (n-1))
    cat(n, '-- ', length(jj1), ' noRA cells -- pst interval :', bins[n], ' --- ', bins[n+1],  '\n')
  }
  
  #jj2 = which(aa$pseudot_noRA < bins[(n+1)] & aa$pseudot_noRA >= bins[n])
  #aa$pst_group[jj2] = paste0('noRA_groupPst_', n)
  
  #cat(n, '-- ', length(jj1), ' RA cells -- ', length(jj2), ' noRA cells \n')
}


# nb_bins = 8
# aa$pst_group = NA
# 
# bins = seq(0.6, 1, length.out = (nb_bins + 1))
# for(n in 1:(length(bins)-1))
# {
#   jj = which(aa$pseudot_RA < bins[(n+1)] & aa$pseudot_RA > bins[n])
#   aa$pst_group[jj] = paste0('RA_groupPst_', n)
# }
# 
# bins = seq(0.25, 1, length.out = (nb_bins + 1))
# for(n in 1:(length(bins)-1))
# {
#   jj = which(aa$pseudot_noRA < bins[(n+1)] & aa$pseudot_noRA > bins[n])
#   aa$pst_group[jj] = paste0('noRA_groupPst_', n)
# }

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.1)
DimPlot(aa, group.by = 'condition', label = TRUE, cols = cols_sel)


bb = subset()

p1 = DimPlot(aa, group.by = 'condition', label = TRUE, cols = cols_sel)
p2 = FeaturePlot(aa, features = 'pseudot', cols = c('lightgray', 'blue')) +
  ggtitle(label = 'pseudo time') +
  scale_color_viridis_c(direction = -1)

p3 = DimPlot(aa, group.by = 'pst_group', label = FALSE, repel = TRUE)

p1 + p2 + p3

ggsave(filename = paste0(resDir, '/RA_noRA_d2_d5_condition_pseudotime_palantir_bins_v2.pdf'), 
       width = 20, height = 6) 


source('functions_utility.R')
hete = calc_heterogeneity_RA.noRA(seuratObj = aa, 
                                  method = 'pairwiseDist_pca', 
                                  subsample.cells = 200, 
                                  #subsample.sketch = FALSE,
                                  nb_features = 1000, 
                                  nb_pcs = 20)

hete$condition = factor(hete$condition, 
                        levels = paste0(c('noRA', 'RA'), '_groupPst_', rep(c(1:(length(bins)-1)), 
                                                                           each = 2)))

hete$treament = 'RA'
hete$treament[grep('noRA_', hete$condition)] = 'noRA'
hete$treament = factor(hete$treament, levels = c('RA', 'noRA'))

ggplot(hete, aes(x= condition, y=dists, fill = treament)) + 
  geom_violin() + 
  #scale_fill_manual(values = ) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 10)) +
  labs( x = '', y = 'pairwise distance' )

ggsave(filename = paste0(resDir, '/pairwise_distance_pst.palantir.Bins_', length(bins), 
                         '.bins.pdf'), width = 12, height = 8) 


########################################################
########################################################
# Section: RA vs noRA using palantir pseudotime
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 
                           'RA_noRA_d2_d5_condition_pseudotimePalantir_saved4heterogeity.rds'))

# aa = readRDS(file = paste0(RdataDir, 
#                            'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
#                            'cellCycleScoring_annot.v1_reduction.DM_princurves_pseudotime_',
#                            'mNT_scRNAseq_R13547_10x_mNT_20220813.rds'))

candidates = readRDS(file = paste0('../results/scRNAseq_R13547_10x_mNT_20220813', 
                                   '/RA.vs.noRA_firstBifurcation/', 
                                   'DElist_1932genes_pairwiseComaprison_v2.rds'))
candidates = unique(candidates$gene)

DimPlot(aa, cols = cols_sel, group.by = 'condition', reduction = 'dm')
#FeaturePlot(aa, features = candidates[1])

cols = readRDS(file = '../results/Rdata/color_scheme_4scRNAseq.rds')
levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", "day3.5_RA",   "day4_RA", "day5_RA",
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA")
cols_sel = cols[match(levels_sels, names(cols))]

source('functions_utility.R')

plot_genes_branched_heatmap(seuratObj = aa, 
                            outDir = resDir,
                            gene_subset = candidates,
                            nbCell_condition = 50,
                            cols_sel = cols_sel,
                            #hmcols = viridis(20, option = "D", direction = -1),
                            plotName = 'pheatmap_RA_noRA_branchinng_expression_pseudotime_test_v2'
                            )


########################################################
########################################################
# Section : signaling pathways for regulative cell proportion
# gloabl pathway analysis and cellcht for LR anlaysis
########################################################
########################################################
library(data.table)
library(plyr)
library(ggplot2)
library(scales)

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam_clusters_4save.rds'))

sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))


levels_sels = unique(aa$condition)
names(cols) = levels
cols_sel = cols[match(levels_sels, names(cols))]

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE, cols = cols_sel)

p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE, cols = cols_cluster) +
  #scale_colour_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) +
  ggtitle('clusters 3 and 5 for signaling pathway analysis')
#labs( x = '', y = '% of FoxA2+ & Pax6+' )

p1 + p2

ggsave(filename = paste0(resDir, '/scRNAseq_RAsamples_signalingPathways_clusterAnnot.pdf'), 
       width = 14, height = 5)


Use_scFates_selectClusters = FALSE
if(Use_scFates_selectClusters){
  ##########################################
  # select the time points and clusters 
  ##########################################
  
  dataDir = paste0("../results/scRNAseq_R13547_10x_mNT_20220813/RA_symetryBreaking/TF_modules/",
                   'd2.5_d5_TFs_SPs_regressed.CellCycle_v1/')
  
  aa = readRDS(file = paste0(dataDir, 
                             'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered',
                             '_regressout.nCounts_cellCycleScoring_annot.v2_newUMAP_clusters_time_',
                             'd2.5_d5_regressed.CellCycle_v1.rds'))
  
  pst = read.csv(file = paste0(dataDir, 'annData_pseudotime_segments_milestones.csv'), header = TRUE,
                 row.names = c(1))
  
  levels_sels = unique(aa$condition)
  
  names(cols) = levels
  cols_sel = cols[match(levels_sels, names(cols))]
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p1 + p2
  
  mm = match(colnames(aa), rownames(pst))
  aa = AddMetaData(aa, metadata = pst[mm, ])
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE, cols = cols_sel)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seg', raster=FALSE)
  p1 + p2
  
  assign = read.csv(file = paste0(dataDir, 'annData_cellAssignment_to_nonIntersectingWindows_8.csv'),
                    header = TRUE, row.names = c(1))
  #assign = t(assign)
  assignment = rownames(assign)[apply(assign, 2, which.max)] 
  names(assignment) = gsub('[.]','-', colnames(assign))
  
  aa$windows = assignment[match(colnames(aa), names(assignment))]
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'milestones', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'windows', raster=FALSE)
  p1 + p2
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE, cols = cols_sel)
  
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'windows', raster=FALSE) +
    scale_colour_brewer(palette = "Set1") +
    theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
          axis.text.y = element_text(angle = 0, size = 12)) +
    ggtitle('scFates windows')
  #labs( x = '', y = '% of FoxA2+ & Pax6+' )
  
  p1 + p2
  
  ggsave(filename = paste0(resDir, '/scRNAseq_RAsamples_signalingPathways_scFatesWindows.pdf'), 
         width = 14, height = 5)
  
  xx = aa
  
  aa = readRDS(file = paste0(RdataDir, 
                             'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam',
                             '_clustered.discardCellcycle.corrrelatedGenes.rds'))
  
  aa = subset(aa, cells = colnames(xx))
  
  aa$windows = xx$windows[match(colnames(aa), colnames(xx))]
  rm(xx)
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam',
                            '_clustered.discardCellcycle.corrrelatedGenes_forSignalingPathways.rds'))
  
}

##########################################
# test cellchat for pathway activity 
# https://htmlpreview.github.io/?https://github.com/
# jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
##########################################
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#outDir = paste0(resDir, '/test/')
#if(!dir.exists(outDir)) dir.create(outDir)
# aa = readRDS(file = paste0(RdataDir, 
#                            'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam',
#                            '_clustered.discardCellcycle.corrrelatedGenes_forSignalingPathways.rds'))
# 
# aa$windows = paste0('w', aa$windows)
# aa$windows = factor(aa$windows)
aa$windows = factor(aa$clusters)

Idents(aa) = aa$windows

data.input <- aa[["RNA"]]@data # normalized data matrix
# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(aa)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat <- createCellChat(object = aa, group.by = "windows", assay = "RNA")

#cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", 
#                           key = "annotation") # use Secreted Signaling
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

cellchat <- identifyOverExpressedGenes(cellchat,)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692

cellchat <- computeCommunProb(cellchat, 
                              #type = "triMean", 
                              type =  "truncatedMean", trim = 0.1)

cellchat <- filterCommunication(cellchat, min.cells = 5, rare.keep = TRUE)

# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat, signaling = c("WNT", "TGFb", 'HH', "ncWNT", 'BMP', 'FGF', 'TGFb'),
                         remove.isolate = FALSE)


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways

# # check the order of cell identity to set suitable vertex.receiver
# levels(cellchat@idents)
# vertex.receiver = seq(1,4)
# 
# for (i in 1:length(pathways.show.all)) {
#   
#   i = 5
#   # Visualize communication network associated with both signaling pathway and individual L-R pairs
#   netVisual_circle(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, 
#             layout = "hierarchy")
#   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#   #gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#   ggsave(filename=paste0(outDir, pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, 
#          width = 3, height = 2, units = 'in', dpi = 300)
#   
# }

# Compute the network centrality scores
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
# Visualize the computed centrality scores using heatmap, 
# allowing ready identification of major signaling roles of cell groups

pathways.show = 'FGF'

netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, 
                                  height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("BMP"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


pathways_sel = c('HH', 'WNT', 'FGF', 'BMP')
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", 
                                         signaling = pathways_sel)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                         signaling = pathways_sel)
ht3 = netAnalysis_signalingRole_heatmap(cellchat, pattern = "all",
                                        signaling = pathways_sel)

ht1 + ht2 + ht3

pdf(paste0(resDir, '/scRNAseq_RAsamples_signalingPathways_CellChat_all_v2.pdf'), 
    height = 6, width = 8, useDingbats = FALSE)

netAnalysis_signalingRole_heatmap(cellchat, pattern = "all",
                                  signaling = pathways_sel)

dev.off()


##########################################
# test global pathway activities with ReactomeGSA 
##########################################
