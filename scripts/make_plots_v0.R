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

resDir = paste0("../results/figures_talbes", version.analysis)
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

aa = readRDS(file = paste0(RdataDir, 'seuratObj_clustersFiltered_umapOverview.rds'))


source(paste0(functionDir, '/functions_scRNAseq.R'))

explore.umap.params.combination(sub.obj = aa, resDir = outDir, 
                                pdfname = 'UMAP_test_allSamples.pdf',
                                use.parallelization = FALSE,
                                group.by = 'condition',
                                cols = cols, 
                                weight.by.var = TRUE,
                                nfeatures.sampling = c(2000, 3000, 5000),
                                nb.pcs.sampling = c(20, 30, 50), 
                                n.neighbors.sampling = c(30, 50, 100),
                                min.dist.sampling = c(0.1, 0.3)
)


p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2


aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)

ggsave(filename = paste0(outDir, 'UMAP_rmDoublet_rmRiboMT_regressed.nCounts_annot.v1_',
                         'subsetting.RAsymmetryBreaking.onlyday3rep1',
                         version.analysis, '.pdf'), 
       width = 10, height = 8)

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                          'cellCycleScoring_annot.v2_',
                          species, version.analysis, '.rds'))



Idents(aa) = factor(aa$condition, levels = levels)
aa$condition = factor(aa$condition, levels = levels)


aa =  readRDS(file = paste0(RdataDir, 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', species, version.analysis, '.rds'))

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

#aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 50, min.dist = 0.1)
#aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.1)
# aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 50, min.dist = 0.2)
# aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 50, min.dist = 0.3)
#aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 100, min.dist = 0.2)
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)


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


DimPlot(aa, cols = cols, group.by = 'condition', label = FALSE, repel = TRUE, raster = FALSE) +
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
              '/UMAP_condition_toUse_noLabel.pdf'),  width=8, height = 6) 


saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                          'cellCycleScoring_annot.v1_savedUMAP.v1_', species, version.analysis, '.rds'))

##########################################
# features and TF activity overlaying UMAP 
##########################################
levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1",  'day3_RA.rep2',
                "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
data_version = "_d2_d2.5_d3_d3.5_d4_d5"

names(cols) = levels
cols_sel = cols[match(levels_sels, names(cols))]

Idents(aa) = factor(aa$condition, levels = levels)
subs = subset(aa, idents = levels_sels)

DimPlot(subs, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)

library(viridis) 
(FeaturePlot(subs, features = c('Foxa2', 'Pax6'), cols = c("lightgrey", "blue"))) 
#&
#  scale_colour_viridis_c(option = "D")

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_Foxa2_Pax6.pdf'),  width=12, height = 6) 

# rerun the umap 
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(subs, ndims = 50)
Idents(subs) = subs$condition

subs <- RunUMAP(subs, dims = 1:30, n.neighbors = 100, min.dist = 0.2)

DimPlot(subs, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE, pt.size = 1.5)

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_updatedUMAP.pdf'),  width=8, height = 6) 

DimPlot(subs, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE, pt.size = 1.5) +
  NoLegend()

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_updatedUMAP.pdf'),  width=7, height = 6) 


DimPlot(subs, cols = cols, group.by = 'condition', label = FALSE, repel = TRUE, raster = FALSE, pt.size = 1.5) +
  NoLegend()

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_updatedUMAP.pdf'),  width=7, height = 6) 

FeaturePlot(subs, features = c('Foxa2', 'Pax6'))

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_Foxa2_Pax6_updatedUMAP.pdf'),  width=12, height = 6) 

FeaturePlot(object = subs,  features = c('Foxa2', 'Pax6'), blend = TRUE,
            cols =  c("lightgrey","#00ff00",  "magenta"))

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_Foxa2_Pax6_updatedUMAP_blend.pdf'),  width=24, height = 6) 

saveRDS(subs, file = paste0(RdataDir, 
              'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
              'cellCycleScoring_annot.v1_savedUMAP.subs.v2_', species, version.analysis, '.rds'))

ggs = c('Cdh1', 'Crabp2', 'Zeb2', 'Rbpj', 'Gli1', 'Shh', 'Tead1', 'Tead2', 'Tead3', "Tead4",
        "Gli2", "Esr1", "Kdm5b", 'Lef1', 'Otx2','Pou5f1', 'Sox1', 'Nanog', 'Nkx6-1', 
        'Nkx2-2', "Olig2", "Pax3", "Pax7", "Arx", "Sox10", "Pbx3", "Prdm1",  "Sox2", "Tcf7l1",
        "Tcf7l2", "Zeb1", "Xbp1", "Tfdp1","Pax6", "Foxa2", "Smad1", "Smad4", "Smad2", "Smad5", "Smad3", "Smad7", 
        "Smad9", "Bmp4", 'Bmpr2', 'Bmpr1b', 'Bmpr1a', "Acvr2a", "Acvr2b", "Smo", "Notch1", 
        "Lrp5", "Lrp6", "Fzd4", 'Fzd8', 'Wnt6', "Wnt5a", "Wnt5b", "Wnt4", "Wnt1", "Wnt3a", "Ptch1",
        "Ptch2", "Lrp2", "Gpc1", "Gas1", "Cdon", "Hhip", "Fgf8", "Cdh2")
mm = match(ggs, rownames(subs))
ggs[which(is.na(mm))]

ggs = ggs[which(!is.na(mm))]

for(g in ggs)
{
  # g = "Foxa2"
  cat(g, "-- \n")
  FeaturePlot(object = subs,  features = g)
  
  ggsave(paste0("../results/plots_MondaySeminar", 
                '/umap_d2.to.d6_Foxa2_Pax6_updatedUMAP_geneExpression.example_', g, '_withLabel.pdf'),
         
         width=8, height = 6) 
  
  FeaturePlot(object = subs,  features = g) + NoLegend()
  
  ggsave(paste0("../results/plots_MondaySeminar", 
                '/umap_d2.to.d6_Foxa2_Pax6_updatedUMAP_geneExpression.example_', g, '.pdf'),
         width=7, height = 6) 
  
}


aa = readRDS(file = paste0(RdataDir, 
                   'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                   'cellCycleScoring_annot.v1_savedUMAP.subs.v2_savedTFactivities.decoupleR.rds'))

DefaultAssay(object = aa) <- "tfswmean"

ggs = c('Cdh1', 'Crabp2', 'Zeb2', 'Rbpj', 'Gli1', 'Shh', 'Tead1', 'Tead2', 'Tead3', "Tead4",
        "Gli2", "Esr1", "Kdm5b", 'Lef1', 'Otx2','Pou5f1', 'Sox1', 'Nanog', 'Nkx6-1', 
        'Nkx2-2', "Olig2", "Pax3", "Pax7", "Arx", "Sox10", "Pbx3", "Prdm1",  "Sox2", "Tcf7l1",
        "Tcf7l2", "Zeb1", "Xbp1", "Tfdp1","Pax6", "Foxa2", "Smad1", "Smad4", "Smad2", "Smad5", "Smad3", "Smad7", 
        "Smad9", "Bmp4", 'Bmpr2', 'Bmpr1b', 'Bmpr1a', "Acvr2a", "Acvr2b", "Smo", "Notch1", 
        "Lrp5", "Lrp6", "Fzd4", 'Fzd8', 'Wnt6', "Wnt5a", "Wnt5b", "Wnt4", "Wnt1", "Wnt3a", "Ptch1",
        "Ptch2", "Lrp2", "Gpc1", "Gas1", "Cdon", "Hhip", "Fgf8", "Cdh2")
mm = match(ggs, rownames(aa))
ggs[which(is.na(mm))]

ggs = ggs[which(!is.na(mm))]

for(gene in ggs)
{
  # gene = "Foxa2"
  cat("--", gene, "-- \n")
  #FeaturePlot(object = aa,  features = g)
  (FeaturePlot(aa, features = gene) & 
      scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) + 
    ggtitle(paste0(gene, ' activity'))
  ggsave(paste0("../results/plots_MondaySeminar/TF_activity", 
                '/umap_d2.to.d6_TFactivity_geneExamples_', gene, '_withLabel.pdf'),
         width=8, height = 6) 
  
  (FeaturePlot(aa, features = gene) & 
      scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) + 
    ggtitle(paste0(gene, ' activity')) + NoLegend()
  
  ggsave(paste0("../results/plots_MondaySeminar/TF_activity", 
                '/umap_d2.to.d6_TFactivity_geneExamples_', gene, '.pdf'),
         width=7, height = 6)
  
}

