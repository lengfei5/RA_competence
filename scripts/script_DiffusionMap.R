##########################################################################
##########################################################################
# Project: run diffusion map
# Script purpose: 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Sep 29 10:59:04 2022
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R13547_10x_mNT_20220813'

resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R13547_10x'
functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(DropletUtils)
library(edgeR)
library(future)
options(future.globals.maxSize = 80000 * 1024^2)
mem_used()

species = 'mNT_scRNAseq'

### Parameters for this project
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

##########################################
# main code:  test Diffusion Map reduction
##########################################
library(pheatmap)
library(RColorBrewer)
library(grid)
library(Seurat)
library(scater)
library(SingleCellExperiment)
library(scran)
library(destiny)

aa =  readRDS(file = paste0(RdataDir, 
                'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                'cellCycleScoring_annot.v1_', species, version.analysis, '.rds'))

Idents(aa) = factor(aa$condition, levels = levels)
aa$condition = factor(aa$condition, levels = levels)

# DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)
#aa = subset(aa, DF_out == 'Singlet')
levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", "day3.5_RA",   "day4_RA", 
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA")
cols_sel = cols[match(levels(aa$condition), names(cols))]


aa = subset(aa, idents = levels_sels)
aa$condition = droplevels(aa$condition)

####
#### not use the processing steps from Seurat
####
#aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
#aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
#ElbowPlot(aa, ndims = 50)

#ll.pca = Embeddings(object = aa, reduction = "pca")

#aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
#ElbowPlot(aa, ndims = 50)

#ll.pca.notWeighted = Embeddings(object = aa, reduction = "pca")

#Idents(aa) = aa$condition
#aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.2)
#DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
#DimPlot(aa, reduction = 'pca', label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
metadata = aa@meta.data
sce = as.SingleCellExperiment(aa)
#rm(aa)
dec <- modelGeneVar(sce)

#plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
#curve(metadata(dec)$trend(x), col="blue", add=TRUE)


for(nb_features in c(3000, 5000, 8000))
{
  for(n_neighbors in c(100, 200, 300, 500))
  {
    top.hvgs <- getTopHVGs(dec, n=nb_features)
    
    sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 100)
    # reducedDimNames(sce)
    ll.pca = reducedDim(sce, 'PCA')[, c(1:50)]
    
    tic()
    dm <- DiffusionMap(ll.pca, sigma = 'local', k = n_neighbors, 
                       n_eigs = 50, distance = 'euclidean')
    
    saveRDS(dm, file = paste0(RdataDir, 'diffusionMap_firstBifurcation_RA.vs.noRA_with.sce.PCA_euclidean_',
                              'localSignal_nfeatures.', nb_features, '_nbNeighbors.', n_neighbors,
                              '.rds'))
    
    toc()
    
    tic()
    dm <- DiffusionMap(ll.pca, sigma = 'global', k = n_neighbors, n_eigs = 50, distance = 'euclidean')
    
    saveRDS(dm, file = paste0(RdataDir, 'diffusionMap_firstBifurcation_RA.vs.noRA_with.sce.PCA_euclidean_',
                              'globalSignal_nfeatures.', nb_features, '_nbNeighbors.', n_neighbors,
                              '.rds'))
    toc()
    
  }
}

##########################################
# downstream analysis
##########################################
# # plot(dm)
# cells = names(dm$DC1)
# xx = data.frame(DC1 = dm$DC1, DC2 = dm$DC2, DC3 = dm$DC3, bc = cells, stringsAsFactors = FALSE)
# xx$condition = metadata$condition[match(cells, rownames(metadata))]
# xx$condition = factor(xx$condition)
# 
# names(cols) = levels
# 
# ggplot(aes(x = DC1, y = DC2, col = condition), data = xx[which(xx$DC2<0.005),]) +
#   geom_point(size = 0.1) + 
#   scale_color_manual(values=cols[match(levels(xx$condition), names(cols))])
# 
# ggplot(aes(x = DC1, y = DC3, col = condition), data = xx) +
#   geom_point()
# 
# ggplot(aes(x = DC2, y = DC3, col = condition), data = xx) +
#   geom_point(size = 0.1)
# 
# #install.packages("plot3D")
# #library("plot3D")
# library("car")
# 
# scatter3d(x = xx$DC1, 
#           y = xx$DC2, 
#           z = xx$DC3, 
#           #groups = factor(xx$condition),
#           col = cols[1:9],
#           surface=FALSE, ellipsoid = FALSE)
# 
# scatter3d(x = sep.l, y = pet.l, z = sep.w, groups = iris$Species, col = 
#             surface=FALSE, ellipsoid = TRUE)
# #library(destiny, quietly = TRUE)
# #dm <- DiffusionMap(sce)
# #rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)
# # plot(rd2, col = topo.colors(100), pch=16, asp = 1)
# #reducedDims(sim) <- SimpleList(PCA = rd1, DiffMap = rd2)
# # dm <- DiffusionMap(ll.pca, sigma = 'local', n_eigs = 5)
# 
# #plot(dm)
# #plot(dm$DC1, dm$DC2)
# dcs = as.matrix(cbind(dm$DC1, dm$DC2))
# ll.obj[["DP"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", assay = DefaultAssay(ll.obj))
# DimPlot(ll.obj, reduction = 'DP', group.by = 'manual.annot.ids')
# 
# 
