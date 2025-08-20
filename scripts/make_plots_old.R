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
# all samples 
##########################################
outDir = paste0(resDir, '/UMAP_allSamples/')
if(!dir.exists(outDir)) dir.create(outDir)

aa =  readRDS(file = paste0('../results/Rdata/', 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', 
                            species, '_R13547_10x_mNT_20220813', '.rds'))
Idents(aa) = factor(aa$condition, levels = levels)
table(aa$condition)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)

ggsave(filename = paste0(resDir, '/UMAP_overview_initial_beforeFiltering.pdf'), 
       width = 12, height = 8)


Filter_weirdCluster_7_8 = FALSE
if(Filter_weirdCluster_7_8){
  cells_2filter = readRDS(paste0(RdataDir, 'cellNames_weirdClusters_7_8.rds'))
  mm = match(colnames(aa), cells_2filter)
  aa = subset(aa, cells = colnames(aa)[which(is.na(mm))])
  
}

source(paste0(functionDir, '/functions_scRNAseq.R'))

explore.umap.params.combination(sub.obj = aa, resDir = outDir, 
                                pdfname = 'UMAP_test_allSamples.pdf',
                                use.parallelization = FALSE,
                                group.by = 'condition',
                                cols = cols, 
                                weight.by.var = TRUE,
                                nfeatures.sampling = c(3000, 5000),
                                nb.pcs.sampling = c(30, 50), 
                                n.neighbors.sampling = c(30, 50, 100),
                                min.dist.sampling = c(0.1, 0.3)
)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
p1 + p2


##########################################
# no early time points
##########################################
outDir = paste0(resDir, '/UMAP_overview_allSamples_fromDay2/')
if(!dir.exists(outDir)) dir.create(outDir)

levels_sels = levels[which(levels != 'day0_beforeRA' & levels != 'day1_beforeRA')]
cols_sel = cols[match(levels_sels, names(cols))]

bb = subset(aa, idents = levels_sels)

source(paste0(functionDir, '/functions_scRNAseq.R'))

explore.umap.params.combination(sub.obj = bb, resDir = outDir, 
                                pdfname = 'UMAP_test_allSamples_no.day0.day1.pdf',
                                use.parallelization = FALSE,
                                group.by = 'condition',
                                cols = cols_sel, 
                                weight.by.var = TRUE,
                                nfeatures.sampling = c(3000, 5000),
                                nb.pcs.sampling = c(30, 50), 
                                n.neighbors.sampling = c(30, 50, 100),
                                min.dist.sampling = c(0.1, 0.3)
)






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

version.analysis = '_R13547_10x_mNT_20240522'
species = 'mNT_scRNAseq'

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_savedUMAP.v1_mNT_scRNAseq_R13547_10x_mNT_20220813.rds'))


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


subs = readRDS(file = paste0(RdataDir, 
              "seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_",
              "regressout.nCounts_cellCycleScoring_annot.v1_savedUMAP.subs.v2_mNT_scRNAseq_R13547_10x_mNT_20220813.rds"
))

FeaturePlot(object = subs,  features = c('Foxa2', 'Pax6'), blend = TRUE
            #cols =  c("lightgray","#00ff00",  "magenta")
            )

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

