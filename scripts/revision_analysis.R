##########################################################################
##########################################################################
# Project: RA competence for Hannah's project
# Script purpose: analysis for the first revision
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Nov 17 13:31:48 2026
##########################################################################
##########################################################################

rm(list = ls())

version.analysis = '_R13547_10x_mNT_20261118'
species = 'mNT_scRNAseq'

resDir = paste0("../results/revision_figures", version.analysis)

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
# Section I : compare the scRNA-seq embryo data with the mNTO
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam_clusters_4save.rds'))

levels_sels = c("day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
cols_sel = cols[match(levels_sels, names(cols))]


DimPlot(aa, label = FALSE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE) +
  NoLegend()


DimPlot(aa, label = TRUE, group.by =  'clusters', repel = TRUE, raster=FALSE, 
        cols = cols_cluster) + NoLegend()

aa$celltypes[which(aa$clusters == 2)] = 'Head'
aa$celltypes[which(aa$clusters == 3)] = 'pFP'
aa$celltypes[which(aa$clusters == 5)] = 'eNP'

DimPlot(aa, label = TRUE, group.by =  'celltypes', repel = TRUE, raster=FALSE) + NoLegend()

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


##########################################
# prepare the mouse gastrulation data  
##########################################
bb = readRDS(file = paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/scRNAseq_mNT/saved_seuratObj/',
                           'seuratObject_EmbryoAtlasData_all36sample.rds'))

celltype_colours = c("Epiblast" = "#635547",
                     "Primitive Streak" = "#DABE99",
                     "Caudal epiblast" = "#9e6762",
                     
                     "PGC" = "#FACB12",
                     
                     "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     
                     "Visceral endoderm" = "#F6BFCB",
                     "ExE endoderm" = "#7F6874",
                     "ExE ectoderm" = "#989898",
                     "Parietal endoderm" = "#1A1A1A"
                     
)

stage_colours = c("E6.5" = "#D53E4F",
                  "E6.75" = "#F46D43",
                  "E7.0" = "#FDAE61",
                  "E7.25" = "#FEE08B",
                  "E7.5" = "#FFFFBF",
                  "E7.75" = "#E6F598",
                  "E8.0" = "#ABDDA4",
                  "E8.25" = "#66C2A5",
                  "E8.5" = "#3288BD"
)

coexpcol = c("Not" = "grey", "CoExp" = "purple")

## subset cell types
Idents(bb) = as.factor(bb$celltype)

VlnPlot(bb, features = c('Foxa2', 'Pax6'), group.by = 'celltype', ncol = 1)

ggsave(filename = paste0(resDir, '/mouse_gastrulation_VlnPlot_Foxa2_Pax6.pdf'), 
       width = 12, height = 8)

## select the follow cell types specified by Hannah
'Spinal cord'
'Forebrain/Midbrain/Hindbrain'
'Rostral neurectoderm'
'Caudal neurectoderm'
"Forebrain/Midbrain/Hindbrain"

## select the cell types expressed either Pax6 or FoxA2
celltypes_sels = c('Anterior Primitive Streak',
                   
                   'Def. endoderm', 
                   'Epiblast', 
                   'Forebrain/Midbrain/Hindbrain',
                   'Gut', 
                   'NMP', 
                   'Notochord',
                   'Primitive Streak',
                   'Rostral neurectoderm',
                   'Spinal cord', 
                   'Surface ectoderm',
                   'Visceral endoderm',
                   
                   'Forebrain/Midbrain/Hindbrain',
                   'Rostral neurectoderm',
                   'Caudal neurectoderm')

celltypes_sels = unique(celltypes_sels)

sub.obj = subset(bb, idents = celltypes_sels)

VlnPlot(sub.obj, features = c('Foxa2', 'Pax6'), group.by = 'celltype', ncol = 1)


FeatureScatter(sub.obj, feature1 = 'Foxa2', feature2 = 'Pax6', 
               group.by = 'stage', cols = stage_colours, raster = FALSE, pt.size = 1.5) +
  theme(plot.title = element_blank()) + 
  labs( x = 'FoxA2', y = 'Pax6' ) +
  theme(axis.text.x = element_text(angle = 0, size = 12),
        axis.text.y = element_text(angle = 0, size = 12))

## or subset the FoxA2+ & Pax6+ cells
# sub.obj = subset(x = bb, subset = Foxa2 > 0 & Pax6 > 0)
# 
# counts = table(sub.obj$celltype)
# 
# pcts = counts
# for(n in 1:length(pcts))
# {
#   pcts[n] = pcts[n]/length(which(bb$celltype == names(pcts)[n]))
# }
# 
# pcts = pcts[order(-pcts)]

sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

sub.obj = ScaleData(sub.obj)

sub.obj <- RunPCA(sub.obj, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(sub.obj, ndims = 50)

#Idents(sub.obj) = sub.obj$condition

sub.obj <- RunUMAP(sub.obj, reduction = 'pca.corrected', dims = 1:30, n.neighbors = 50, min.dist = 0.1)

p1 = DimPlot(sub.obj, group.by = 'celltype', label = TRUE, repel = TRUE, raster = FALSE)

p2 = FeaturePlot(sub.obj, features = c('Foxa2', 'Pax6'))

p1 / p2


