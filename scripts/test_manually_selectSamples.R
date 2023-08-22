##########################################################################
# Project: RA competence  
# Script purpose: manually select samples for the downstream analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Oct  7 14:09:17 2022
##########################################################################
##########################################################################

##########################################
# drop the day2.5 sample to focus on the symmetry breaking (day3, day3.5 and day4)
##########################################
Remove_day2.5_day6 = FALSE
if(Remove_day2.5_day6){
  aa = readRDS(file = paste0(RdataDir, 
                             'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                             'cellCycleScoring_annot.v2_newUMAP_clusters_noRAday6_',
                             species, version.analysis, '.rds'))
  aa$celltypes = aa$clusters
  aa$clusters = aa$seurat_clusters
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  p1 + p2
  
  
  Idents(aa) = aa$condition
  aa = subset(aa, cells = colnames(aa)[which(aa$condition != 'day2.5_RA')])
  
  # rerun the umap 
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
  #DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)
  #ggsave(filename = paste0(outDir, 'UMAP_RAsymmetryBreaking.onlyday3rep1_3000HVGs_noweighted.byvarPCA_',
  #                         '30pcs_30neighbors_minDist0.1.pdf'), 
  #       width = 10, height = 8)
  
  # quickly run clustering
  #ElbowPlot(aa, ndims = 50)
  aa <- FindNeighbors(aa, dims = 1:30)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p1 + p2
  
  ggsave(filename = 
           paste0(outDir, "noday6_noday2.5_downsample_v5/",
                  'UMAP_RAsymmetryBreaking.onlyday3rep1_timePoints_clustering.res0.7_noRAday6_noDay2.5.pdf'), 
         width = 20, height = 8)
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)
  
  aa$clusters = aa$seurat_clusters
  
  aa$celltypes = aa$clusters
  aa$celltypes = as.character(aa$celltypes)
  aa$clusters = as.character(aa$clusters)
  aa$clusters[which(aa$celltypes == "FP")] = "8"
  aa$clusters[which(aa$celltypes == "NP")] = "6"
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                            'cellCycleScoring_annot.v2_newUMAP_clusters_noRAday6_noRAday2.5_',
                            species, version.analysis, '.rds'))
  
  
}
