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

##########################################
# redo umap for RA d2.5-d6 
##########################################
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
                                pdfname = 'UMAP_test_RASamples_selectedParameters.pdf',
                                use.parallelization = FALSE,
                                group.by = 'condition',
                                cols = cols, 
                                weight.by.var = TRUE,
                                nfeatures.sampling = c(3000, 5000),
                                nb.pcs.sampling = c(30, 50, 100), 
                                n.neighbors.sampling = c(30, 50, 100, 200), 
                                min.dist.sampling = c(0.1)
                                
)

saveRDS(aa, file = paste0(RdataDir, 'seuratObj_clustersFiltered_umap_RAsamples_nod2.rds'))

##########################################
# cluster the RA sample 
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObj_clustersFiltered_umap_RAsamples_nod2.rds'))

levels_sels = c("day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
cols_sel = cols[match(levels_sels, names(cols))]

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 30, min.dist = 0.1, spread = 1, seed.use = 45)

DimPlot(aa, label = FALSE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE) +
  coord_flip() + scale_x_reverse()


ggsave(filename = paste0(resDir, '/scRNAseq_overviewUMAP_RAsamples', version.analysis, '.pdf'), 
       width = 10, height = 8)


saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam.rds'))

### test some clustering options
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.8)
DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)

DimPlot(aa, group.by = "Phase", label = TRUE, repel = TRUE, raster=FALSE)

knownGenes =  c('Pou5f1', 'Sox2', 'Zfp42',
                'Otx2', 'Cyp26a1', 'Stra8',
                'Hoxa1', 'Hoxa3', 'Hoxb4',
                'Sox1', 'Pax6', 'Tubb3', 
                'Foxa2', 'Shh', 'Arx', 
                'Nkx2-2', 'Olig2', 'Pax3')

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

cols_cluster = c( "#ff7f00", "#08519c", "#762a83", "chartreuse4", "#e31a1c", "#20b2aa","#b30024")
DimPlot(aa, label = TRUE, group.by =  'clusters', repel = TRUE, raster=FALSE, 
        cols = cols_cluster) 

ggsave(filename = paste0(resDir, '/scRNAseq_overview_RAsamples_clustering.pdf'), 
       width = 10, height = 8)

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam',
                          '_clustered.discardCellcycle.corrrelatedGenes.rds'))

##########################################
# highlight the marker genes or specified genes
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam',
                           '_clustered.discardCellcycle.corrrelatedGenes.rds'))

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

xx$clusters = droplevels(xx$clusters)

DoHeatmap(xx, group.by = 'clusters', features = top10$gene, draw.lines = TRUE, disp.min = -2.,
          group.colors = cols_cluster, angle = 0, size = 5) + 
  theme(text = element_text(size = 6)) +
  #NoLegend() +
  scale_fill_viridis_c(option = "magma")

  #scale_fill_viridis_c() 
  #scico::scale_fill_scico(palette = "vik")
  #scale_fill_viridis(option = "D") 
ggsave(filename = paste0(resDir, '/scRNAseq_RAsamples_clustering_allMarkers.pdf'), 
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
genes.to.plot <- knownGenes


library(RColorBrewer)
colGEX = c("grey85", brewer.pal(7, "Reds"))
DotPlot(aa, group.by = "clusters", features = genes.to.plot[length(genes.to.plot):1]) + 
  coord_flip() + scale_color_gradientn(colors = colGEX) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0))

ggsave(filename = paste0(resDir, '/scRNAseq_overview_RAsamples_clustering_markers_dotplot_sorted.pdf'), 
       width = 8, height = 4)

nClust = 7
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(nClust)
VlnPlot(aa, group.by = "clusters", fill.by = "ident", cols = cols_cluster,
              features = genes.to.plot, stack = TRUE, flip = TRUE)

ggsave(filename = paste0(resDir, '/scRNAseq_overview_RAsamples_clustering_markers_Vlnplot.pdf'), 
       width = 8, height = 6)


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
                                 '#b30000', '#7f0000'))

p1 = plotDR(sce, "DiffusionMap", color_by = "clusters") +
  theme_classic() + 
  scale_colour_manual(values = c("#464646", "#7F7F7F", "#CD00CF", "#FFC000", "#70AD47"))

p0 + p1 

ggsave(paste0(figureDir, 'FACS_RA_DiffusionMap_timePoints_clusters_colors.pdf'), width=14, height = 6)


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
  
  ggs = c()
  xx = read.csv(file = paste0(dataDir, 'geneBasis_top50_tfs_sps_final.csv'), sep = ";")
  ggs = unique(c(ggs, xx[, 2]))
  
  write.csv(xx, file = paste0(resDir, '/geneBasis_top50_tfs_sps_final.csv'), row.names = FALSE)
  
  xx = read.csv(file = paste0(dataDir, 'dubstep_sparse_74.tfs.sps_final.csv'))
  ggs = unique(c(ggs, xx[, 1]))
  
  write.csv(xx, file = paste0(resDir, '/dubstep_sparse_74.tfs.sps_final.csv'), row.names = FALSE)
  
  xx = read.csv(file = paste0(dataDir, 'SMD_sparse_31tfs.sps_final.csv'))
  ggs = unique(c(ggs, xx[c(1:31), 1]))
  write.csv(xx[c(1:31), ], 
            file = paste0(resDir, '/SMD_sparse_31tfs.sps_final.csv'), row.names = FALSE)
  
  aa = readRDS(paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/',
                      'RA_competence/scRNAseq_mNT/saved_seuratObj/',
                      'RAsamples_d2_d6_UMAP_selectedParam.rds'))
  
  ggs = ggs[which(!is.na(match(ggs, rownames(aa))))]
  
  
  levels_sels = c("day2_beforeRA",  
                  "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA")
  
  #levels_sels = unique(aa$condition)
  cols_sel = cols[match(levels_sels, names(cols))]
  Idents(aa) = as.factor(aa$condition)
  
  aa = subset(aa, idents = levels_sels)
  aa <- RunPCA(aa, features = ggs, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 100, min.dist = 0.1)
  DimPlot(aa, group.by = 'condition', label = TRUE, cols = cols_sel)
  
  ggsave(filename = paste0(resDir, '/RA_umap_with_selectedFeatures.pdf'), width = 8, height = 6) 
  
  saveRDS(ggs, file = paste0(resDir, '/geneList_sparseFeatureSelection.rds'))
  
}

##########################################
# Visualising RA-specific expression pattern / kinetics (feature plots / exp vs pseudotime?  
# by Hannah
##########################################
aa = readRDS(paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/',
                    'scRNAseq_mNT/saved_seuratObj/',
                    'RA_noRAsamples_d2_d6_UMAP_selectedParam.rds'))

ggs = readRDS(file = paste0(resDir, '/geneList_sparseFeatureSelection.rds'))

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
aa$groups = factor(aa$groups, levels = c('RA', 'noRA'))

aa$condition = droplevels(aa$condition)

aa$time = sapply(aa$condition, function(x) {unlist(strsplit(as.character(x), '_'))[1]})

features =c('Foxa2', 'Pax6', 'Pou5f1', 'Sox1')
VlnPlot(aa, features = features, split.by = "groups", group.by = 'time', log = FALSE, 
        same.y.lims = TRUE, ncol = 2, pt.size = 0, 
        cols = c("#EF6548", "#6BAED6"))

ggsave(filename = paste0(resDir, '/RAspecific_geneExpression_kinetics_noDot_v2.pdf'), 
       width = 16, height = 8) 


##########################################
# Dotplot for RA vs noRA
##########################################
# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
aa = readRDS(paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/',
                    'scRNAseq_mNT/saved_seuratObj/',
                    'RA_noRAsamples_d2_d6_UMAP_selectedParam.rds'))

ggs = readRDS(file = paste0(resDir, '/geneList_sparseFeatureSelection.rds'))

levels_sels = c("day2_beforeRA",
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA", "day6_RA", 
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA", "day6_noRA")

Idents(aa) = as.factor(aa$condition)
aa = subset(aa, idents = levels_sels)

levels_sels = unique(aa$condition)
cols_sel = cols[match(levels_sels, names(cols))]

DimPlot(aa, cols = cols_sel)


aa$condition = droplevels(aa$condition)

aa$condition = factor(as.character(aa$condition), 
                      levels = c("day6_noRA", "day5_noRA", "day4_noRA", 
                                 "day3.5_noRA", "day3_noRA",  "day2.5_noRA",
                                 "day2_beforeRA", 
                                 "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   
                                 "day4_RA", "day5_RA", "day6_RA")
)

#Idents(aa) = factor(aa$time) 
DotPlot(aa, features = ggs, 
        group.by = "condition"
) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) +
  RotatedAxis() + 
  coord_flip() 

ggsave(filename = paste0(resDir, '/RA_vs_noRA_sparseFeatures.pdf'), width = 8, height = 28) 

#DotPlot_scCustom(aa, features = ggs[c(1:20)], split.by = "groups",  flip_axes = TRUE)


test_scCutomize_clustering = FALSE
if(test_scCutomize_clustering){
  require(scCustomize)
  aa$orig.ident = factor(aa$condition)
  Clustered_DotPlot(aa, features = c('Pax6', 'Sox1'), split.by = 'condition',
                    plot_km_elbow = FALSE)
  
}



########################################################
########################################################
# Section II: RA vs noRA
# 
########################################################
########################################################
saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                          'cellCycleScoring_annot.v1_reduction.DM_princurves_',
                          species, version.analysis, '.rds'))

saveRDS(res_kmean, file = paste0(RdataDir, 'DM_princurves_clusterCenter.rds'))

##########################################
# try to merge two principle curves 
##########################################
pcurve_noRA = readRDS(file = paste0(outDir, 'principle_curve_noRA_v2.rds'))
pcurve_RA = readRDS(file = paste0(outDir, 'principle_curve_RA_borderCellsSelected_v3.rds'))
res_kmean = readRDS(file = paste0(RdataDir, 'DM_princurves_clusterCenter.rds'))

dcs_all = data.frame(aa[['DC']]@cell.embeddings[, c(1,2)])

cluster_sels = unique(aa$dc_clusters[!is.na(aa$dc_clusters)])
mm = which(!is.na(match(aa$dc_clusters, cluster_sels)))
dcs = dcs_all[mm, ]


#cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[aa$dc_clusters[mm]]
mm = match(rownames(dcs), colnames(aa))
cols = cols_sel[match(aa$condition[mm], names(cols_sel))]

pdf(paste0(outDir, "Two_principleCurves_noRA_RA_v3.pdf"),
    height = 8, width =10, useDingbats = FALSE)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3, 4, 2, 1), tcl = -0.1)

plot(dcs, col = cols, cex = 0.1, main = 'RA vs. noRA trajectories')
lines(pcurve_RA$s[order(pcurve_RA$lambda),], lty=1,lwd=4,col="red",type = "l")

xx = pcurve_noRA$s[order(pcurve_noRA$lambda), ]
#pt = pseudotime.scaling(pcurve_noRA$lambda[order(pcurve_noRA$lambda)])
jj = which(xx[, 1] <  min(pcurve_RA$s))

lines(xx[jj, ], lty=1,lwd=4,col="black",type = "l")
lines(xx[-jj, ], lty=1,lwd=4,col="blue",type = "l")


dev.off()


##########################################
# plot heatmap  
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_reduction.DM_princurves_pseudotime_',
                           species, version.analysis, '.rds'))

candidates = readRDS(file = paste0(outDir, 'DElist_1932genes_pairwiseComaprison_v2.rds'))

candidates = unique(candidates$gene)

DimPlot(aa, cols = cols_sel, group.by = 'condition', reduction = 'DC')
FeaturePlot(aa, features = candidates[1])

source('functions_utility.R')
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_reduction.DM_princurves_pseudotime_',
                           species, version.analysis, '.rds'))

candidates = readRDS(file = paste0(outDir, 'DElist_1932genes_pairwiseComaprison_v2.rds'))

candidates = unique(candidates$gene)

DimPlot(aa, cols = cols_sel, group.by = 'condition', reduction = 'DC')
FeaturePlot(aa, features = candidates[1])

source('functions_utility.R')
# seuratObj = aa; nbCell_condition = 100;scale_max=3; scale_min=-3;hmcols = NULL; Get.Smooth.Curve = TRUE;
# gene_subset = candidates;hclust_method = "ward.D2";num_clusters = 6
library(VGAM) # an example code from https://online.stat.psu.edu/stat504/lesson/8/8.2/8.2.2
library(MASS)
library(tidyr)
require(colorRamps)
require(pheatmap)

cell_beforeRA = colnames(seuratObj)[which(!is.na(seuratObj$pseudot) & 
                                            seuratObj$condition == 'day2_beforeRA')]
cell_noRA = colnames(seuratObj)[which(!is.na(seuratObj$pseudot) & grepl('_noRA', seuratObj$condition))]
cell_RA = colnames(seuratObj)[which(!is.na(seuratObj$pseudot) & grepl('_RA', seuratObj$condition))]

cat('subsampling ', nbCell_condition, ' cells\n')
cell.sels = c()
cc = unique(seuratObj$condition)
cc = cc[which(cc != "day3_RA.rep2")]

for(n in 1:length(cc))
{
  cell.sels = c(cell.sels, sample(colnames(seuratObj)[which(seuratObj$condition == cc[n] & 
                                                              !is.na(seuratObj$pseudot))], 
                                  size = nbCell_condition, 
                                  replace = FALSE))
}

subs = subset(seuratObj, cells = cell.sels)

get_smooth_curve_spline = function(x, t, newt, downsample = TRUE)
{
  # x = as.numeric(cds[1, jj]); t = Pseudotime; newt = pseudot_comomon;
  if(downsample){
    nb_t = min(5000, length(t))
    nn = sample(1:length(t), size = nb_t, replace = FALSE)
    t = t[nn]
    x = x[nn]
  }
  
  fit_sel = smooth.spline(t, x, df = 3)
  
  #plot(Pseudotime, cds_sel, cex = 0.5)
  #lines(fit_sel, col = 'red', lwd =2.0)
  newx = predict(fit_sel, newt)
  return(newx$y)
  #VGAM::vglm(~sm.ns(Pseudotime, df=3), family = 'gaussian', data = cds_sel)
  
}

if(Get.Smooth.Curve){
  cds <- seuratObj@assays$RNA@scale.data
  cds = cds[which(!is.na(match(rownames(cds), gene_subset))), ]
  cat(' -- smoothing the single cell data for subsampled cells -- \n')  
  
  # before RA
  jj = match(cell_beforeRA, colnames(seuratObj))
  jj = jj[which(!is.na(seuratObj$pseudot[jj]))]
  Pseudotime = as.numeric(seuratObj$pseudot[jj])
  kk_common = which(subs$condition == 'day2_beforeRA')
  kk_common = kk_common[order(subs$pseudot[kk_common])]
  pseudot_comomon = subs$pseudot[kk_common]
  
  common_ancestor_cells = t(apply(cds[ ,jj], 1, get_smooth_curve_spline, t = Pseudotime, newt = pseudot_comomon))
  
  jj = grep('_RA$|_RA.rep1', seuratObj$condition)
  jj = jj[which(!is.na(seuratObj$pseudot[jj]))]
  Pseudotime = as.numeric(seuratObj$pseudot[jj])
  
  kk_BrachA = grep('_RA$|_RA.rep1', subs$condition)
  kk_BrachA = kk_BrachA[order(subs$pseudot[kk_BrachA])]
  pseudot_BrachA = subs$pseudot[kk_BrachA]
  
  BranchA_exprs <- t(apply(cds[ ,jj], 1, get_smooth_curve_spline, t = Pseudotime, newt = pseudot_BrachA))
  
  jj = grep('_noRA', seuratObj$condition)
  jj = jj[which(!is.na(seuratObj$pseudot[jj]))]
  Pseudotime = as.numeric(seuratObj$pseudot[jj])
  
  kk_BrachB = grep('_noRA', subs$condition)
  kk_BrachB = kk_BrachB[order(subs$pseudot[kk_BrachB])]
  pseudot_BrachB = subs$pseudot[kk_BrachB]
  
  BranchB_exprs <- t(apply(cds[ ,jj], 1, get_smooth_curve_spline, t = Pseudotime, newt = pseudot_BrachB))
  
}

col_gap_ind <- c(length(kk_BrachB), length(kk_BrachB) + length(kk_common), 
                 length(kk_BrachB) + 2*length(kk_common))

heatmap_matrix <- cbind(BranchB_exprs[, ncol(BranchB_exprs):1], 
                        common_ancestor_cells[, ncol(common_ancestor_cells):1],
                        common_ancestor_cells,
                        BranchA_exprs)

indexs = c(kk_BrachB[ncol(BranchB_exprs):1], 
           kk_common[ncol(common_ancestor_cells):1], 
           kk_common,
           kk_BrachA)

heatmap_matrix=heatmap_matrix[!apply(heatmap_matrix, 1, sd)==0,]
heatmap_matrix=Matrix::t(scale(Matrix::t(heatmap_matrix), center=TRUE))
heatmap_matrix=heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE, ]
heatmap_matrix[is.nan(heatmap_matrix)] = 0
heatmap_matrix[heatmap_matrix>scale_max] = scale_max
heatmap_matrix[heatmap_matrix<scale_min] = scale_min

saveRDS(heatmap_matrix, file = paste0(outDir, "/heatmap_matrix_forPlot.rds"))
#heatmap_matrix_ori <- heatmap_matrix
#heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & is.finite(heatmap_matrix[, col_gap_ind]), ] #remove the NA fitting failure genes for each branch 

row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1

exp_rng <- range(heatmap_matrix) #bks is based on the expression range

bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
if(is.null(hmcols)) {
  hmcols <- blue2green2red(length(bks) - 1)
}

# prin  t(hmcols)
ph <- pheatmap(heatmap_matrix, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=TRUE, 
               show_rownames=F, 
               show_colnames=F, 
               #scale="row",
               clustering_distance_rows=row_dist, #row_dist
               clustering_method = hclust_method,
               cutree_rows=num_clusters,
               silent=TRUE,
               #filename=NA,
               breaks=bks,
               color=hmcols
               #color=hmcols#
)

annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)),
                             pseudot = subs$pseudot[indexs],
                             condition = subs$condition[indexs])
#names(annotation_colors$`Cell Type`) = c('Pre-branch', branch_labels)
feature_label <- row.names(heatmap_matrix)
row_ann_labels <- row.names(annotation_row)

row.names(heatmap_matrix) <- feature_label
row.names(annotation_row) <- row_ann_labels

##########################################
# DE TFs and signaling pathways 
##########################################
tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))
sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v2.rds'))
targets = unique(c(tfs, sps$gene))
xx = read.table('../data/annotations/GO_term_summary_RAR_signalingPathway.txt', header = TRUE, sep = '\t', 
                row.names = NULL)
targets = unique(c(targets, xx[,2]))
xx = read.table('../data/annotations/GO_term_summary_TGFb.txt', header = TRUE, sep = '\t', 
                row.names = NULL)
targets = unique(c(targets, xx[,2]))
#sps = toupper(unique(sps$gene))
#sps = setdiff(sps, tfs)

sels = which(!is.na(match(rownames(heatmap_matrix), targets)))
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix[sels,])))/2)
row_dist[is.na(row_dist)] <- 1
annotation_rowSel = data.frame(Cluster = annotation_row[sels, ])
rownames(annotation_rowSel) = rownames(annotation_row)[sels]

pheatmap(heatmap_matrix[sels, ], #ph$tree_row$order
         useRaster = T,
         cluster_cols=FALSE, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         show_colnames=FALSE, 
         scale='none',
         clustering_distance_rows=row_dist, #row_dist
         clustering_method = hclust_method, #ward.D2
         cutree_rows=num_clusters,
         # cutree_cols = 2,
         annotation_row=annotation_rowSel,
         annotation_col=annotation_col,
         annotation_colors=annotation_colors,
         gaps_col = col_gap_ind,
         treeheight_row = 30, 
         breaks=bks,
         fontsize = 6,
         color=hmcols, 
         border_color = NA,
         silent=TRUE, 
         filename=paste0(outDir, "/expression_pseudotime_pheatmap_allDEgenes_TF.SP.pdf"),
         width = 6, height = 12
)

pheatmap(heatmap_matrix[sels, ], #ph$tree_row$order
         useRaster = T,
         cluster_cols=FALSE, 
         cluster_rows=TRUE, 
         show_rownames=TRUE,
         show_colnames=FALSE, 
         scale='none',
         clustering_distance_rows=row_dist, #row_dist
         clustering_method = hclust_method, #ward.D2
         cutree_rows=num_clusters,
         # cutree_cols = 2,
         annotation_row=annotation_rowSel,
         annotation_col=annotation_col,
         annotation_colors=annotation_colors,
         gaps_col = col_gap_ind,
         treeheight_row = 30, 
         breaks=bks,
         fontsize_row = 4,
         #fontsize = 4,
         color=hmcols, 
         border_color = NA,
         silent=TRUE, 
         filename=paste0(outDir, "/expression_pseudotime_pheatmap_allDEgenes_TF.SP",
                         "_withgeneNames.pdf"),
         width = 8, height = 20
)


########################################################
########################################################
# Section :  Quantify the cell hetereogeneity with RA and without RA  
# 
########################################################
########################################################
Use_pseudotime_palantir = FALSE
if(Use_pseudotime_palantir){
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
  
  ## add DM components
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

saveRDS(aa, file = paste0(RdataDir, 'RA_noRA_d2_d5_condition_pseudotime_saved4heterogeity.rds'))

##########################################
###### test heterogeneity quantification by binning pseudotime
##########################################
aa = readRDS(file = paste0(RdataDir, 'RA_noRA_d2_d5_condition_pseudotime_saved4heterogeity.rds'))

levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA",
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA")

cols_sel = cols[match(levels_sels, names(cols))]


p1 = DimPlot(aa, group.by = 'condition', label = TRUE, cols = cols_sel)
p2 = FeaturePlot(aa, features = 'pseudot', cols = c('lightgray', 'blue')) +
  ggtitle(label = 'pseudo time') +
  scale_color_viridis_c(direction = -1)

p1 + p2

ggsave(filename = paste0(resDir, '/RA_noRA_d2_d5_condition_pseudotime_palantir_v2.pdf'), 
       width = 10, height = 6) 


VlnPlot(aa, features = 'pseudot', group.by = 'condition', pt.size = 0.00, cols = cols_sel) +
  geom_hline(yintercept = c(0.05), col = 'red') +
  ggtitle(label = 'pseudotime distribution per condition')

ggsave(filename = paste0(resDir, '/RA_noRA_d2_d5_condition_pseudotime_v2.pdf'), width = 10, height = 6) 

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

