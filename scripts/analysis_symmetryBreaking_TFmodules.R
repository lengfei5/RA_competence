##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: search for genes with asymmetric expression at early time points
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Mar 30 16:04:20 2023
##########################################################################
##########################################################################
library(Seurat)
library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(future)

options(future.globals.maxSize = 120000 * 1024^2)

source('functions_utility.R')

outDir = paste0(resDir, '/RA_symetryBreaking/TF_modules/')
system(paste0('mkdir -p ', outDir))

levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA")
data_version = "_d2_d2.5_d3_d3.5_d4_d5"

names(cols) = levels
cols_sel = cols[match(levels_sels, names(cols))]


#aa = readRDS(file = paste0(RdataDir,
#                           'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
#                           'cellCycleScoring_annot.v1_savedUMAP.subs.v2_', species, version.analysis, '.rds'))
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d6_',
                           species, version.analysis, '.rds'))

##########################################
# import the all data for RA treatment
##########################################
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE, pt.size = 1)

ggsave(filename = paste0(outDir, '/UMAP_conditions.pdf'), 
       width = 10, height = 6)

Idents(aa) = aa$condition
p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)

p1 + p2

ggsave(filename = paste0(outDir, 'UMAP_RAtreatment_d2.to.d6_.pdf'), 
       width = 18, height = 8)

FeaturePlot(aa, features = c('Zfp42', 'Tcf15', 'Skil', 'Lef1'))
ggsave(filename = paste0(outDir, '/asymmetric_feature_expression.pdf'), 
       width = 14, height = 10)

FeaturePlot(aa, features = c('Pax6', 'Foxa2', 'Sox1'))
ggsave(filename = paste0(outDir, '/asymmetric_feature_expression_v2.pdf'), 
       width = 14, height = 10)

p1 = FeaturePlot(aa, features = c("Pax6", "Sox1"), blend = TRUE)
p2 = FeaturePlot(aa, features = c("Foxa2", "Sox1"), blend = TRUE)

p1 /p2

ggsave(filename = paste0(outDir, '/Overlapping_Foxa2_Pax6_Sox1.pdf'), 
       width = 14, height = 8)


FeatureScatter(aa, feature1 = "Foxa2", feature2 = "Pax6", group.by = 'condition')
FeatureScatter(aa, feature1 = "Foxa2", feature2 = "Sox1", group.by = 'condition')
FeatureScatter(aa, feature1 = "Pax6", feature2 = "Sox1", group.by = 'condition')


##########################################
# clean and subset the samples for scFates test 
##########################################
Clean_Subset_for_scFates = FALSE
if(Clean_Subset_for_scFates){
  ## remove the cluster 8 and 9 mainly mature neurons and also day6_RA
  aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '8' & aa$celltypes != '9' & 
                                               aa$condition != 'day6_RA')])
  
  aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '8' & aa$celltypes != '9')])
  
  aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '8' & aa$celltypes != '9' & 
                                               aa$condition != 'day2_beforeRA')])
  
  # rerun the umap 
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 100, min.dist = 0.1)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)
  
  ggsave(filename = paste0(outDir, 'UMAP_RAtreatment_d2.to.d6.noMatureNeurons.pdf'), 
         width = 10, height = 6)
  
  
  Discard_cellCycle.corrrelatedGenes = TRUE
  if(Discard_cellCycle.corrrelatedGenes){
    library(scater)
    Idents(aa) = aa$condition
    
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
    
    print(intersect(genes_discard, gene_examples))
    
    aa = subset(aa, features = setdiff(rownames(aa), genes_discard))
    
  }
  
  aa = FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  Idents(aa) = aa$condition
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 100, min.dist = 0.2)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)
  
  saveRDS(aa, file = paste0(outDir, 
                            'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                            'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d6.noNeurons.rds'))
  
  
  ##########################################
  # select the time point to test 
  ##########################################
  aa = readRDS(file = paste0(outDir, 
                             'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                             'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d5.noNeurons.rds'))
  
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)
  
  aa$dataset = 'afterRA'
  aa$dataset[which(aa$condition == 'day2.5_RA')] = 'RA'
  aa$dataset[which(aa$condition == 'day2_beforeRA')] = 'beforeRA'
  
  ## test first the day3, day3.5, day4 and day5
  Select_timePoints_for_scFates = FALSE
  if(Select_timePoints_for_scFates){
    
    Idents(aa) = as.factor(aa$condition)
    
    aa = subset(aa, idents = c('day2.5_RA', 'day3_RA.rep1', 'day3.5_RA', 'day4_RA', 'day5_RA', 'day6_RA'))
    
    aa = FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)
    
    aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
    ElbowPlot(aa, ndims = 50)
    
    Idents(aa) = aa$condition
    aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 50, min.dist = 0.1)
    DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)
    
    ggsave(filename = paste0(outDir, 
                             'UMAP_RAtreatment_d2.5.to.d6_no.mautreNeurons_filtered.cellCycleGenes.pdf'), 
           width = 10, height = 6)
    
    library(SeuratDisk)
    mnt = aa
    
    VariableFeatures(mnt) = NULL
    #mnt@assays$RNA@scale.data = NULL
    #mnt@assays$RNA@data = NULL
    
    DefaultAssay(mnt) = 'RNA'
    
    mnt = DietSeurat(mnt, 
                     counts = TRUE, 
                     data = TRUE,
                     scale.data = FALSE,
                     features = rownames(mnt), 
                     assays = c('RNA'), 
                     dimreducs = c('umap', 'pca'), graphs = NULL, 
                     misc = TRUE
    )
    
    
    DefaultAssay(mnt) = 'RNA'
    VariableFeatures(mnt)
    
    
    #saveRDS(mnt, file = paste0(outDir, '/SeuratObj_splice.unspliced_RA_day2_to_day6_all.rds'))
    Idents(mnt) = mnt$condition
    
    mnt = subset(mnt, downsample = 2000)
        
    saveFile = '/RNAmatrix_RA_d2.5_d6_all.h5Seurat'
    
    SaveH5Seurat(mnt, filename = paste0(outDir, saveFile), 
                 overwrite = TRUE)
    Convert(paste0(outDir, saveFile), 
            dest = "h5ad", overwrite = TRUE)
    
  }
  
  
}

########################################################
########################################################
# Section I : narrow down the TFs and SPs
# 
########################################################
########################################################
require(scater)
require(scran)

ggs = rownames(aa)
ggs = ggs[which(!is.na(match(ggs, c(tfs, sps))))]

sce <- as.SingleCellExperiment(subset(aa, features = ggs))

#rownames(sce) = toupper(rownames(sce))
ave.counts <- calculateAverage(sce, assay.type = "counts")

hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))

num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
              xlab=expression(Log[10]~"average count"))

# detected in >= 5 cells, ave.counts >=5 but not too high
genes.to.keep <- num.cells > 100 & ave.counts >= 10^-2
summary(genes.to.keep)

sce <- sce[genes.to.keep, ]

ggs = rownames(sce)

length(ggs[!is.na(match(ggs, tfs))])



########################################################
########################################################
# Section : subclustering each time points
# 
########################################################
########################################################
p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
#p2 = DimPlot(aa, group.by = 'clusters', label = TRUE, repel = TRUE)

conditions = unique(aa$condition)
print(conditions)
bb = aa;
rm(aa)

conditions = c('day2_beforeRA', 'day2.5_RA', 'day3_RA.rep1')
noisyGenes = readRDS(file = paste0(RdataDir, 'topGenes_localVaribility.gene.expression_VarID2.rds'))
gene_examples = unique(c(gene_examples, noisyGenes))

for(cc in conditions)
{
  #cc = c('day2_beforeRA')
  # cc = c('day2.5_RA')
  # cc = c('day3_RA.rep1')
  # cc = c('day3.5_RA')
  # cc = c('day4_RA')
  # cc = c('day5_RA')
  #cc = c('day3_RA.rep1', 'day3.5_RA')
  cc = c('day2_beforeRA', 'day2.5_RA', 'day3_RA.rep1')
  
  outDir_cc = paste0(outDir, '/', paste0(cc, collapse = "_"), '/')
  system(paste0('mkdir -p ', outDir_cc))
  
  aa = subset(bb, cells = colnames(bb)[which(!is.na(match(bb$condition, cc)))])
  
  Idents(aa) = aa$condition
  table(aa$condition)
 
  DimPlot(aa, group.by = 'condition', label = TRUE)
  FeaturePlot(aa, features = c('Pax6', 'Foxa2'))
  
  ##########################################
  # UMAP and clustering without sparse feature selection
  ##########################################
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000) # find subset-specific HVGs
  variableGenes = VariableFeatures(object = aa)
  select.method = 'HVGs2000_pca.weighted_'
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  #variableGenes = setdiff(VariableFeatures(object = aa), 'Lypd2')
  # variableGenes = intersect(variableGenes, c(tfs, sps))
  
  #variableGenes = VariableFeatures(object = aa)
  cat(length(variableGenes), ' variable genes using \n ')
  
  aa <- RunPCA(aa, features = variableGenes, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) + NoLegend()
  p11 = FeaturePlot(aa, features = 'nCount_RNA')
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(outDir_cc, '/UMAP_RA_condition_clustering_cellcyclePhase_', select.method, '.pdf'), 
         width = 12, height = 6)
  
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(outDir_cc, '/featureExamples_originalUMAP.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = aa, 
                           features = gene_examples[!is.na(match(gene_examples, rownames(aa)))])
  
  dev.off()
  
  ##########################################
  # regress out the cell cycle  
  ##########################################
  xx = aa
  xx$CC.Difference <- xx$S.Score - xx$G2M.Score
  #xx <- ScaleData(xx, vars.to.regress = c('nCount_RNA', "CC.Difference"), features = variableGenes)
  
  xx <- ScaleData(xx, vars.to.regress = c('nCount_RNA', "S.Score", "G2M.Score"), 
                  features = variableGenes)
  
  # cell cycle effects strongly mitigated in PCA
  xx <- RunPCA(xx, features = variableGenes, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(xx, ndims = 50)
  
  xx <- RunUMAP(xx, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
  DimPlot(xx, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  
  xx <- FindNeighbors(xx, dims = 1:10)
  xx <- FindClusters(xx, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  p1 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) + NoLegend()
  p11 = FeaturePlot(xx, features = 'nCount_RNA')
  p2 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(outDir_cc, '/UMAP_RA_condition_clustering_cellcyclePhaseRegression_', 
                           select.method, '.pdf'), 
         width = 12, height = 6)
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(outDir_cc, '/featureExamples_UMAP_cellcycleRegression.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = xx, 
                           features = gene_examples[!is.na(match(gene_examples, rownames(aa)))])
  
  dev.off()
  
  rm(xx)
  
  ##########################################
  # Removing cell cycle-related genes in the HVGs
  # original code :
  # https://bioconductor.org/books/3.12/OSCA/cell-cycle-assignment.html#removing-cell-cycle-effects
  ##########################################
  library(scater)
  #require(batchelor)
  
  Idents(aa) = aa$condition
  
  # Identifying the likely cell cycle genes between phases,
  # using an arbitrary threshold of 5%.
  scaledMatrix = GetAssayData(aa, slot = c("scale.data"))
  
  diff <- getVarianceExplained(scaledMatrix, data.frame(phase = aa$Phase))
  diff = data.frame(diff, gene = rownames(diff))
  diff = diff[order(-diff$phase), ]
  
  hist(diff$phase, breaks = 100); abline(v = c(1:5), col = 'red')
  
  genes_discard = diff$gene[which(diff$phase>5)]
  cat(length(genes_discard), 'genes to discard \n')
  
  hvgs = VariableFeatures(FindVariableFeatures(subset(aa, features = setdiff(rownames(aa), genes_discard)), 
                                               selection.method = "vst", nfeatures = 2000))
  # 
  # discard <- diff > diff_cutoff
  # summary(discard)
  # cat('cut off -- ', diff_cutoff, '--', sum(discard), ' genes related to cell cycle \n')
  # genes.discard = rownames(diff)[which(discard)]
  aa <- RunPCA(aa, 
               features = hvgs, 
               verbose = FALSE, 
               weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) + NoLegend()
  p11 = FeaturePlot(aa, features = 'nCount_RNA')
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(outDir_cc, '/UMAP_RA_condition_clustering_cellcycle.rmCellCycleGenes_', 
                           select.method, '.pdf'),  width = 12, height = 6)
  
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(outDir_cc, '/featureExamples_UMAP_cellcycle.rmCellCycleGenes.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = aa, 
                           features = gene_examples[!is.na(match(gene_examples, rownames(aa)))])
  
  dev.off()
  
  # Idents(aa) = aa$seurat_clusters
  # all.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
  # all.markers %>%
  #   group_by(cluster) %>%
  #   top_n(n = 10, wt = avg_log2FC) -> top10
  # 
  # DoHeatmap(aa, features = top10$gene) + NoLegend()
  # ggsave(filename = paste0(outDir_cc, '/Heatmap_clusterMarkers_', select.method, '.pdf'), 
  #        width = 14, height = 20)
  
  
}
