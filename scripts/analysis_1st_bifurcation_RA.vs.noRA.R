##########################################################################
##########################################################################
# Project:
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Oct  7 13:30:11 2022
##########################################################################
##########################################################################
names(cols) = levels
levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", "day3.5_RA",   "day4_RA", 
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA")

cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/RA.vs.noRA_firstBifurcation/')
system(paste0('mkdir -p ', outDir))

##########################################
# prepare and process the before RA, w/o RA day2 - day4 
##########################################
aa =  readRDS(file = paste0(RdataDir, 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', species, version.analysis, '.rds'))

Idents(aa) = factor(aa$condition, levels = levels)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)

ggsave(filename = paste0(resDir, '/UMAP_cellCyclePhase.pdf'), width = 10, height = 8)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)
ggsave(filename = paste0(resDir, '/UMAP_rmDoublet_rmRiboMT_regressed.nCounts_annot.v1.pdf'), 
       width = 10, height = 8)

aa = subset(aa, idents = levels_sels)


aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.2)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                          'cellCycleScoring_annot.v1_',
                          species, version.analysis, '.rds'))

##########################################
# test different approach of PCA calculation and diffusion map parameters
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_',
                           species, version.analysis, '.rds'))


# weighted PCA is not working well with DM
#dm = readRDS(file = paste0(RdataDir, 'diffusion_map_2save_allConditions.rds')) 
#dm = readRDS(file = paste0(RdataDir, 'diffusion_map_with.unweightedPCA_firstBifurcation_RA.vs.noRA.rds')) 
#dm = readRDS(file = paste0(RdataDir, 'diffusion_map_with.weightedPCA_firstBifurcation_RA.vs.noRA.rds'))  
#dm = readRDS(file = paste0(RdataDir, 'diffusionMap_firstBifurcation_RA.vs.noRA_with.sce.PCA_euclidean.rds'))
#dm = readRDS(file = paste0(RdataDir, 'diffusionMap_firstBifurcation_RA.vs.noRA_with.sce.PCA_cosine.rds'))

#dm = readRDS(file = paste0(RdataDir, 
#                           'diffusionMap_firstBifurcation_RA.vs.noRA_with.sce.PCA_euclidean_global.sigma.rds'))
require(destiny)

for(nb_features in c(3000, 5000, 8000))
{
  for(n_neighbors in c(100, 200, 300, 500))
  {
    # nb_features = 3000; n_neighbors = 200
    cat('nb_feature -- ', nb_features, '; n_neightors -- ', n_neighbors, '\n')
    dm = readRDS(file = paste0(RdataDir, 'diffusionMap_firstBifurcation_RA.vs.noRA_with.sce.PCA_euclidean_',
                               'globalSignal_nfeatures.', nb_features, '_nbNeighbors.', n_neighbors,
                               '.rds'))
    
    cells = names(dm$DC1)
    metadata = aa@meta.data
    dcs = data.frame(DC1 = dm$DC1, DC2 = dm$DC2, DC3 = dm$DC3, DC4 = dm$DC4, DC5 = dm$DC5, stringsAsFactors = FALSE)
    dcs = dcs[match(rownames(metadata), cells), ]
    
    dcs = as.matrix(dcs)
    aa[["DC"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", assay = DefaultAssay(aa))
    
    rm(metadata)
    
    p1 = DimPlot(aa, reduction = 'DC', dims = c(1, 2), cols = cols_sel)
    p2 = DimPlot(aa, reduction = 'DC', dims = c(1, 3), cols = cols_sel)
    p3 = DimPlot(aa, reduction = 'DC', dims = c(2, 3), cols = cols_sel)
    
    p4 = DimPlot(aa, reduction = 'DC', dims = c(1, 4), cols = cols_sel)
    
    p5=  DimPlot(aa, reduction = 'DC', dims = c(2, 4), cols = cols_sel)
    p6 = DimPlot(aa, reduction = 'DC', dims = c(3, 4), cols = cols_sel)
    
    (p1 + p2 + p3)/(p4 + p5 + p6)
    
    ggsave(filename = paste0(outDir, 'Diffusion_Map_components_test_globalSigma_nFeatures.', nb_features, 
                             '_nNeighbors.', n_neighbors,  '.pdf'), width = 16, height = 10)
    
    DimPlot(aa, reduction = 'DC', dims = c(1, 2), cols = cols_sel)
    ggsave(filename = paste0(outDir, 'Diffusion_Map_components_forTrajectory_globalSigma_nFeatures.', 
                             nb_features, 
                             '_nNeighbors.', n_neighbors,  '.pdf'), 
           width = 16, height = 8)
    
    # Diffusion pseudotime calculation (very slow)
    # Set index or tip of pseudotime calculation to be a zygotic cell (cell 268). 
    #dpt <- DPT(dm)
    #df <- data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2], 
    #                 dptval = dpt$dpt)
    #p1 <- ggplot(df) + geom_point(aes(x = DC1, y = DC2, color = dptval))
    
  }
}

saveRDS(aa, file = paste0(RdataDir, 
                           'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_reduction.DM.globalSignal.nfeature3000.neighbor200_',
                           species, version.analysis, '.rds'))

########################################################
########################################################
# Section : identify trojectory and pseudotime
# 
########################################################
########################################################
library(slingshot, quietly = FALSE)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(scater)
library(SingleCellExperiment)
library(scran)
library(RColorBrewer)

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_reduction.DM.globalSignal.nfeature3000.neighbor200_',
                           species, version.analysis, '.rds'))

## filter the weird clusters identified in make_plots_v0.R
cells_2filter = readRDS(file = paste0(RdataDir, 'subObj_clusters_to_filter.rds'))
mm = match(colnames(aa), colnames(cells_2filter))
aa = subset(aa, cells = colnames(aa)[which(is.na(mm))])

##########################################
# define clusters as slingshot and outlier detection because the principle curves are sensitive to the outliers
# manually correction of clusters 
##########################################
sce = as.SingleCellExperiment(aa)

rd1 = aa[['pca']]@cell.embeddings[, c(1:2)]
rd2 <- aa[['DC']]@cell.embeddings[, c(1:2)]
reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

plot(rd2, col = topo.colors(100), pch=16, cex = 0.1)

#cl1 <- Mclust(rd2)$classification
#colData(sce)$GMM <- cl1
#plot(rd2, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1, cex = 0.2)
set.seed(2022)
res_kmean = kmeans(rd2, centers = 14)
cl2 <- res_kmean$cluster

centers <- res_kmean$centers[res_kmean$cluster, ] 
distances <- sqrt(rowSums((rd2 - centers)^2))

hist(distances, breaks = 100)
abline(v = 0.003, col = 'red')
outliers <- which(distances>0.003)

cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[res_kmean$cluster]

pdfname = paste0(outDir, 'DMcomponents_clustering_outliersDetection_slingshot.pdf')
pdf(pdfname, width=16, height = 8)

plot(rd2, pch=19, col=cols, cex=0.1)
points(res_kmean$centers, col='darkred', pch=15, cex=1)
points(rd2[outliers, ], pch="+", col = 'darkgray', cex=0.8)

dev.off()

aa$dc_clusters_outliers = NA
aa$dc_clusters_outliers[outliers] = 1

colData(sce)$kmeans <- cl2

aa$dc_clusters = cl2
aa$dc_clusters[!is.na(aa$dc_clusters_outliers)] = NA


DimPlot(aa, reduction = 'DC', group.by = 'condition', cols = cols_sel)
ggsave(filename = paste0(outDir, 'DMcomponents_clustering_for_trajectory.pdf'), 
       width = 16, height = 8)

table(aa$condition, aa$dc_clusters)

### manually change the cluster labels if necessary
DimPlot(aa, reduction = 'DC', group.by = 'dc_clusters', label = TRUE)

cl2[which(cl2 == 9 | cl2 == 14)] = 9 # day4_noRA 
cl2[which(cl2== 2 | cl2 == 13 | cl2==11)] = 2 # day4_RA
colData(sce)$kmeans <- cl2

aa$dc_clusters = cl2
aa$dc_clusters[!is.na(aa$dc_clusters_outliers)] = NA
DimPlot(aa, reduction = 'DC', group.by = 'dc_clusters',  label = TRUE)

ggsave(filename = paste0(outDir, 
                         'DMcomponents_clustering.manaulCorrection_for_trajectory_slingshot_outliers_v2.pdf'), 
       width = 16, height = 8)

p1 = DimPlot(aa, reduction = 'DC', group.by = 'dc_clusters', label = TRUE)
p2 = DimPlot(aa, reduction = 'DC', group.by = 'condition', label = TRUE)
p1 + p2

ggsave(filename = paste0(outDir, 
                         'DMcomponents_clustering.manaulCorrection_for_trajectory_outliers',
                         '_clusters_condition_v2.pdf'), 
       width = 20, height = 8)



# default singshot 
#sce <- slingshot(sce, clusterLabels = 'kmeans', reducedDim = 'DC')
#summary(sce$slingPseudotime_1)
#colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
#plot(reducedDims(sce)$DC, col = colors[cut(sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
#lines(SlingshotDataSet(sce), lwd=2)
table(aa$condition, aa$dc_clusters)

if(trajectory_pseudotime_method == 'slingshot'){ #
  
  kk_clean = which(!is.na(aa$dc_clusters))
  
  rd3 = rd2[-outliers, ]
  cl3 = cl2[-outliers]
  
 
    
  plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
  lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')
  
  
  pdf(paste0(outDir, "Slingshot_getLineage_withStarting.cluster_v2.pdf"),
      height = 6, width =10, useDingbats = FALSE)
  
  # times = c('12' = 3, 
  #           '4' = 1, '6' = 1, 
  #           '8' = 2, '3' = 2, '13' = 2, '5' = 2,
  #           '7' = 3, '9' = 3, '1' = 3,
  #           '10' = 4) 
  #lin1 <- getLineages(rd3, cl3, start.clus = '12', end.clus = '1', times = times)
  lin1 <- getLineages(rd3, cl3, start.clus = '4', end.clus = c('2', '9'))
  lin1
  
  cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl3]
  plot(rd3, col = cols, asp = 1, pch = 16, cex = 0.2)
  lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')
  
  dev.off()
  
  crv1 <- getCurves(lin1, shrink = 1, stretch = 2,  smoother = "smooth.spline")
  #crv1 <- getCurves(lin1, shrink = TRUE, stretch = 2,  smoother = "loess")
  crv1
  
  plot(rd3, col = cols, asp = 1, pch = 16, cex = 0.2)
  lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')
  
  pdf(paste0(outDir, "Slingshot_getCurve_withLineage.pdf"),
      height = 6, width =10, useDingbats = FALSE)
  
  
  
  dev.off()
  
  
}


if(trajectory_pseudotime_method == 'principleCurve'){
  
  ##########################################
  # use the principle curve for each trajectory and 
  # merge those two later
  ##########################################
  library(destiny)
  library(princurve)
  
  pseudotime.scaling = function(X) {
    return((X - min(X))/diff(range(X)))
  }
  
  # save the pseudotime estimation 
  #pdfname = paste0(resDir, "/pseudotime_estimation_v1.pdf")
  #pdf(pdfname, width=12, height = 10)
  #par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  aa$pseudot_noRA = NA
  aa$pseudot_RA = NA 
  
  DimPlot(aa, reduction = 'DC', group.by = 'dc_clusters', label = TRUE)
  table(aa$condition, aa$dc_clusters)
  
  dcs_all = data.frame(aa[['DC']]@cell.embeddings[, c(1,2)])
  
  ## no RA conditions
  cluster_sels = c(4, 6, 12, 7, 5, 9)
  mm = which(!is.na(match(aa$dc_clusters, cluster_sels)))
  dcs = dcs_all[mm, ]
  
  cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[aa$dc_clusters[mm]]
  plot(dcs, col = cols, cex = 0.2)
  
  princurve = principal_curve(x = as.matrix(dcs[, c(1,2)]), 
                              #start = start,
                              smoother = 'smooth_spline', stretch = 2)
  
  saveRDS(princurve, file = paste0(outDir, 'principle_curve_noRA_v2.rds')) # save the principle curve found
  
  princurve = readRDS(file = paste0(outDir, 'principle_curve_noRA_v2.rds'))
  
  plot(dcs[, c(1:2)], cex = 0.1, col = cols)
  #points(start, cex = 0.4, col = 'green')
  lines(princurve$s[order(princurve$lambda),], lty=1,lwd=4,col="red",type = "l")
  
  pseudot = pseudotime.scaling(princurve$lambda)
  
  aa$pseudot_noRA[match(names(pseudot), colnames(aa))] = pseudot
  VlnPlot(aa, features = 'pseudot_noRA', group.by = 'condition')
  
  ## RA trajectory: special selection of cells from the border of day2_noRA and day2.5_noRA 
  VlnPlot(aa, features = 'pseudot_noRA', group.by = 'dc_clusters')
  range(aa$pseudot_noRA[which(aa$dc_clusters == 4)])
  range(aa$pseudot_noRA[which(aa$dc_clusters == 6)])
  
  ## select the cells around the cluster 4 centers
  p1 = DimPlot(aa, reduction = 'DC', group.by = 'dc_clusters', label = TRUE)
  p2 = DimPlot(aa, reduction = 'DC', group.by = 'condition', label = TRUE)
  p1 + p2
  
  res_kmean$centers
  head(dcs)  
  dd = apply(dcs, 1, function(x){return(sqrt((x[1] + 0.005982639)^2 + (x[2] - 0.0037537778)^2));})
  dd = dd[order(dd)]
  
  ntop = 1000
  #length(which(aa$pseudot_noRA > 0.18 & aa$pseudot_noRA < 0.219))
  index_startCells = match(names(dd)[1:ntop], colnames(aa))
  
  cluster_sels = c(8, 3, 10,  2,  1)
  #cluster_sels = cluster_sels[!is.na(cluster_sels)]
  mm = which(!is.na(match(aa$dc_clusters, cluster_sels)))
  mm = unique(c(index_startCells, mm))
  dcs = dcs_all[mm, ]
  
  plot(dcs, cex = 0.2)
  
  cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[aa$dc_clusters[mm]]
  plot(dcs, col = cols, cex = 0.2)
  
  princurve = principal_curve(x = as.matrix(dcs[, c(1,2)]), 
                              start = NULL,
                              smoother = 'smooth_spline', 
                              stretch = 2,
                              maxit = 50
                              )
  
  plot(dcs, cex = 0.1)
  #points(start, cex = 0.4, col = 'green')
  lines(princurve$s[order(princurve$lambda),], lty=1,lwd=4,col="red",type = "l")
  
  saveRDS(princurve, file = paste0(outDir, 'principle_curve_RA_borderCellsSelected_v3.rds'))
  
  princurve = readRDS(file = paste0(outDir, 'principle_curve_RA_borderCellsSelected_v3.rds'))
  
 
  pseudot = pseudotime.scaling(princurve$lambda)
  
  aa$pseudot_RA[match(names(pseudot), colnames(aa))] = pseudot
  aa$pseudot_RA[1] = aa$pseudot_noRA[1]
  VlnPlot(aa, features = 'pseudot_RA', group.by = 'condition') # seurat issue if the first element is NA
  
  saveRDS(aa, file = paste0(RdataDir, 
                             'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                             'cellCycleScoring_annot.v1_reduction.DM_princurves_',
                             species, version.analysis, '.rds'))
  
  ##########################################
  # try to merge two principle curves and 
  ##########################################
  pcurve_noRA = readRDS(file = paste0(outDir, 'principle_curve_noRA.rds'))
  pcurve_RA = readRDS(file = paste0(outDir, 'principle_curve_RA_borderCellsSelected.rds'))
  
  dcs_all = data.frame(aa[['DC']]@cell.embeddings[, c(1,2)])
  
  cluster_sels = unique(aa$dc_clusters[!is.na(aa$dc_clusters)])
  mm = which(!is.na(match(aa$dc_clusters, cluster_sels)))
  dcs = dcs_all[mm, ]
  
  
  #cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[aa$dc_clusters[mm]]
  mm = match(rownames(dcs), colnames(aa))
  cols = cols_sel[match(aa$condition[mm], names(cols_sel))]
  
  pdf(paste0(outDir, "Two_principleCurves_noRA_RA_v2.pdf"),
      height = 8, width =10, useDingbats = FALSE)
  
  plot(dcs, col = cols, cex = 0.1)
  
  xx = pcurve_noRA$s[order(pcurve_noRA$lambda),]
  pt = pseudotime.scaling(pcurve_noRA$lambda[order(pcurve_noRA$lambda)])
  jj = which(pt<0.16)
    
  lines(xx[jj, ], lty=1,lwd=4,col="black",type = "l")
  lines(xx[-jj, ], lty=1,lwd=4,col="black",type = "l")
  lines(pcurve_RA$s[order(pcurve_RA$lambda),], lty=1,lwd=4,col="black",type = "l")
  
  dev.off()
  
  
  ##########################################
  # cell weights and pseudo time correction/scaling
  # pseudo time of no_RA trajectory is not touched
  # pseudo time of RA trajectory will be scaled using no_RA as reference
  ##########################################
  aa = readRDS(file = paste0(RdataDir, 
                            'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_reduction.DM_princurves_',
                            species, version.analysis, '.rds'))
  
  pcurve_noRA = readRDS(file = paste0(outDir, 'principle_curve_noRA.rds'))
  pcurve_RA = readRDS(file = paste0(outDir, 'principle_curve_RA_borderCellsSelected.rds'))
  
  # assign cell weights to two trajectories
  index_cells = which(!is.na(aa$pseudot_noRA) | !is.na(aa$pseudot_RA))
  cellWeights = matrix(NA, ncol = 2, nrow = length(index_cells))
  rownames(cellWeights) = colnames(aa)[index_cells]
  colnames(cellWeights) = c('curve_noRA', 'curve_RA')
  
  cell_beforeRA = colnames(aa)[which(!is.na(aa$pseudot_noRA) & aa$condition == 'day2_beforeRA')]
  cell_noRA = colnames(aa)[which(!is.na(aa$pseudot_noRA) & grepl('_noRA', aa$condition))]
  cell_RA = colnames(aa)[which(!is.na(aa$pseudot_RA) & grepl('_RA', aa$condition))]
  
  kk = match(cell_beforeRA, rownames(cellWeights))
  cellWeights[kk, 1] = 1; cellWeights[kk, 2] = 1
  kk1 = match(cell_noRA, rownames(cellWeights))
  cellWeights[kk1, 1] = 1; cellWeights[kk1, 2] = 0
  kk2 = match(cell_RA, rownames(cellWeights))
  cellWeights[kk2, 1] = 0; cellWeights[kk2, 2] = 1
  
  rm(list =(c('kk', 'kk1', 'kk2')))
  
  index_keep = which(!is.na(cellWeights[,1]) & !is.na(cellWeights[,2])) 
  cellWeights = cellWeights[index_keep, ]
  
  ### assign pseudotime for two trajectories and scale the RA ones
  pseudotime = matrix(NA, ncol = 2, nrow = nrow(cellWeights))
  rownames(pseudotime) = rownames(cellWeights)
  colnames(pseudotime) = colnames(cellWeights)
  
  kk = match(cell_beforeRA, rownames(pseudotime)) # cell before RA shared by two trajectories
  cat(length(which(is.na(kk))), 'missing cells \n')
  pseudotime[kk, 1] = aa$pseudot_noRA[match(cell_beforeRA, colnames(aa))]
  pseudotime[kk, 2] = aa$pseudot_noRA[match(cell_beforeRA, colnames(aa))]
  max_pst_shared = max(aa$pseudot_noRA[match(cell_beforeRA, colnames(aa))])
  
  kk1 = match(cell_noRA, rownames(pseudotime)) # cell in noRA trajectory
  cat(length(which(is.na(kk1))), 'missing cells \n')
  pseudotime[kk1, 1] = aa$pseudot_noRA[match(cell_noRA, colnames(aa))] # keep the pseudo computed before
  
  dcs = data.frame(aa[['DC']]@cell.embeddings[, c(1,2)])
  mm = match(cell_noRA, rownames(dcs))
  dcs = as.matrix(dcs[mm, ])
  ptc = project_to_curve(dcs, pcurve_RA$s, stretch = 2)
  
  plot(dcs)
  lines(pcurve_RA$s)
  segments(dcs[, 1], dcs[, 2], ptc$s[, 1], ptc$s[, 2])
  
  ptc_scaled = scales::rescale(ptc$lambda, to = c(0.4, 1))  
  pseudotime[kk1, 2] = ptc_scaled 
  
  kk2 = match(cell_RA, rownames(pseudotime)) # cell in RA trajectory
  cat(length(which(is.na(kk2))), 'missing cells \n')
  ptx = aa$pseudot_RA[match(cell_RA, colnames(aa))]
  pseudotime[kk2, 1] = ptx # keep the pseudo computed before for RA
  
  dcs = data.frame(aa[['DC']]@cell.embeddings[, c(1,2)])
  mm = match(cell_RA, rownames(dcs))
  dcs = as.matrix(dcs[mm, ])
  ptc = project_to_curve(dcs, pcurve_noRA$s, stretch = 2)
  
  plot(dcs)
  lines(pcurve_RA$s)
  segments(dcs[, 1], dcs[, 2], ptc$s[, 1], ptc$s[, 2])
  
  ptc_scaled = scales::rescale(ptc$lambda, to = c(0.4, 1))  
  pseudotime[kk2, 2] = ptc_scaled 
  
  
  rm(list =(c('kk', 'kk1', 'kk2')))
  
  counts = GetAssayData(aa, slot = 'counts')
  counts = counts[, match(rownames(pseudotime), colnames(counts))]
  save(counts, pseudotime, cellWeights, 
       file = paste0(outDir, '/counts_pseudotime_cellWeights_for_tradeSeq.Rdata'))
  
  
}


########################################################
########################################################
# Section : DE gene analysis  
# either by pairwise comparison of samples
# or test trajectory/pseudotime-based DE (not working)
########################################################
########################################################

##########################################
# Method 1) DE candidates with pairwise comparisons between samples
##########################################
require(Seurat)
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_reduction.DM_princurves_',
                           species, version.analysis, '.rds'))


aa$condition[which(aa$condition == 'day3_RA.rep1')] = 'day3_RA'
Idents(aa) = factor(aa$condition)

candidates = c()
days = c('day2.5', 'day3', 'day3.5', 'day4')

for(n in 1:length(days))
{
  # n = 1
  markers = FindMarkers(aa, 
                        ident.1 = paste0(days[n], '_RA'), ident.2 = paste0(days[n], '_noRA'),
                        only.pos = FALSE, 
                        test.use = "wilcox",
                        min.pct = 0.1, 
                        logfc.threshold = 0.25)
  
  markers = data.frame(gene = rownames(markers), markers, 
                       comparison = rep(paste0(days[n], '_RA_vs_', days[n], '_noRA'), length.out = nrow(markers)),
                       stringsAsFactors = FALSE)
  
  candidates = rbind(candidates, markers)
  cat(nrow(markers), ' DE gened found -- ', length(unique(candidates$gene)), ' DE total \n')
  
}

colnames(candidates)[ncol(candidates)] = 'comparison' 

cc = levels(Idents(aa))
cc = cc[grep('day2_beforeRA|day3_RA.rep2', cc, invert = TRUE)]

for(n in c(1:2))
{
  # n = 2
  cat(n, '--', cc[n], '\n')
  markers = FindMarkers(aa, 
                        ident.1 = cc[n], ident.2 = 'day2_beforeRA',
                        only.pos = FALSE, 
                        test.use = "wilcox",
                        min.pct = 0.1, 
                        logfc.threshold = 0.25)
  
  markers = data.frame(gene = rownames(markers), markers, 
                       comparison = rep(paste0(cc[n], '_vs_day2_beforeRA'), length.out = nrow(markers)),
                       stringsAsFactors = FALSE)
  
  candidates = rbind(candidates, markers)
  cat(nrow(markers), ' DE gened found -- ', length(unique(candidates$gene)), ' DE total \n')
  
}

saveRDS(candidates, file = paste0(outDir, 'DElist_1961genes_pairwiseComaprison.rds'))

##########################################
# Method 2) DE with tradeSeq using pseudotime
##########################################
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)

load(file = paste0(outDir, '/counts_pseudotime_cellWeights_for_tradeSeq.Rdata'))

#RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
#data(countMatrix, package = "tradeSeq")
counts <- as.matrix(counts)
#rm(countMatrix)
#data(crv, package = "tradeSeq")
#data(celltype, package = "tradeSeq")


### Downstream of any trajectory inference method using pseudotime and cell weights
# slow, but still ok. with 60K cells, it takes 10-15 mins for each knot.
# tic()
# set.seed(7)
# index_sub = sample(c(1:ncol(counts)), 10000)
# set.seed(7)
# icMat <- evaluateK(counts = counts[, index_sub], 
#                     pseudotime = pseudotime[index_sub, ], 
#                     cellWeights = cellWeights[index_sub, ],
#                     k=3:10, 
#                     nGenes = 300, 
#                     verbose = TRUE, 
#                     plot = TRUE,
#                     parallel=TRUE, 
#                    BPPARAM = BPPARAM)
# 
# saveRDS(icMat, file = paste0(outDir, 'tradeSeqDE_evaluateK_res_v2.rds'))
# 
# toc()


pdfname = paste0(outDir, "Knot_number_estimation_genes300.pdf")
pdf(pdfname, width=12, height = 10)
par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

#icMat = readRDS(file = paste0(outDir, 'tradeSeqDE_evaluateK_res.rds'))
#icMat = readRDS(file = paste0(outDir, 'tradeSeqDE_evaluateK_res_k3.12_100genes_v7.rds'))
#icMat = readRDS(file = paste0(outDir, 'tradeSeqDE_evaluateK_res_k3.12_200genes_v6.rds'))
icMat = readRDS(file = paste0(outDir, 'tradeSeqDE_evaluateK_res_k3.12_300genes_v5.rds'))
icMat = data.frame(icMat)
colnames(icMat) = gsub('k..', '', colnames(icMat))

icDev = icMat - apply(icMat, 1, mean)
boxplot(icDev, ylim = c(-50, 50))

plot(as.numeric(as.character(colnames(icMat))), apply(icMat, 2, mean, na.rm = TRUE), type = 'b')

icRel =  icMat/apply(icMat, 1, max, na.rm = TRUE)
plot(as.numeric(as.character(colnames(icMat))), apply(icRel, 2, mean, na.rm = TRUE), type = 'b')

dev.off()

# downsample the cells and also the genes of interest to test fitGAM
# library(BiocParallel)
# BPPARAM <- BiocParallel::bpparam()
# BPPARAM # lists current options
# BPPARAM$workers <- 16 # use 2 cores
# 
# 
# set.seed(7)
# nb_cells = 10000
# set.seed(7)
# index_sub = sample(c(1:ncol(counts)), nb_cells)
# 
# candidates = readRDS(file = paste0(outDir, 'DElist_3512genes_pairwiseComaprison.rds'))
# genes.sel = match(candidates, rownames(counts))
# genes.sel = genes.sel[which(!is.na(genes.sel))]
# length(genes.sel)
# 
# nb.knots = 6;
# 
# 
# tic()
# set.seed(7)
# sce <- fitGAM(counts = counts[, index_sub], 
#               pseudotime = pseudotime[index_sub, ], 
#               cellWeights = cellWeights[index_sub, ],
#               genes = genes.sel,
#               nknots = nb.knots, 
#               verbose = TRUE, 
#               parallel=TRUE, 
#               BPPARAM = BPPARAM
#               )
# toc()
#save(sce, file = paste0(RdataDir, 'fitGAM_output_tradeSeq_v3.Rdata'))

## reload the fitGam results
load(paste0(outDir, 'tradeSeqDE_fitGAM_output_cellnb.10000.Rdata'))
table(rowData(sce)$tradeSeq$converged)
candidates = readRDS(file = paste0(outDir, 'DElist_3512genes_pairwiseComaprison.rds'))
load(file = paste0(outDir, '/counts_pseudotime_cellWeights_for_tradeSeq.Rdata'))

counts = counts[match(candidates, rownames(counts)), match(colnames(sce), colnames(counts))] 
pseudotime = pseudotime[match(colnames(sce), colnames(counts)), ] 
cellWeights = cellWeights[match(colnames(sce), colnames(counts)), ] 

plotGeneCount(curve = pseudotime, 
              counts = counts,
              #clusters = apply(slingClusterLabels(crv), 1, which.max),
              models = sce)
earlyDERes <- earlyDETest(sce, knots = c(1, 2))

oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
earlyDERes = earlyDERes[oEarly, ]
head(rownames(earlyDERes))

plotSmoothers(sce, counts, gene = rownames(earlyDERes)[oEarly][1])

plotSmoothers(sce, counts, gene = 'Zfp703', nPoints = 1000)
plotSmoothers(sce, counts, gene = 'Cyp26a1', nPoints = 1000)

assoRes <- associationTest(sce)
head(assoRes)

startRes <- startVsEndTest(sce)

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[3]]
plotSmoothers(sce, counts[, index_sub], gene = sigGeneStart)

endRes <- diffEndTest(sce)
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[1]]
plotSmoothers(sce, counts[, index_sub], sigGene)

plotGeneCount(crv, counts[, index_sub], gene = sigGene)

patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat])

plotSmoothers(sce, counts[, index_sub], gene = rownames(patternRes)[oPat][4])

########################################################
########################################################
# Section : visualize the gene candidates with/without intersecting the RAR genomic targets
# 
########################################################
########################################################
require(Seurat)
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_reduction.DM_princurves_',
                           species, version.analysis, '.rds'))
#aa$condition[which(aa$condition == 'day3_RA.rep1')] = 'day3_RA'
Idents(aa) = factor(aa$condition)


candidates = readRDS(file = paste0(outDir, 'DElist_1961genes_pairwiseComaprison.rds'))

### scale the pseudotime of RA trajectory
aa$pseudot = NA

cell_beforeRA = colnames(aa)[which(!is.na(aa$pseudot_noRA) & aa$condition == 'day2_beforeRA')]
cell_noRA = colnames(aa)[which(!is.na(aa$pseudot_noRA) & grepl('_noRA', aa$condition))]
cell_RA = colnames(aa)[which(!is.na(aa$pseudot_RA) & grepl('_RA', aa$condition))]

kk = match(cell_beforeRA, colnames(aa)) # cell before RA shared by two trajectories
cat(length(which(is.na(kk))), 'missing cells \n')
aa$pseudot[kk] = aa$pseudot_noRA[kk]
max_pst_shared = max(aa$pseudot_noRA[kk])

kk1 = match(cell_noRA, colnames(aa)) # cell in noRA trajectory
cat(length(which(is.na(kk1))), 'missing cells \n')
aa$pseudot[kk1] = aa$pseudot_noRA[kk1] # keep the pseudo computed before

kk2 = match(cell_RA, colnames(aa)) # cell in RA trajectory
cat(length(which(is.na(kk2))), 'missing cells \n')
ptx = aa$pseudot_RA[kk2]
aa$pseudot[kk2] = scales::rescale(ptx, to = c(max_pst_shared, 1))

VlnPlot(aa, features = 'pseudot', group.by = 'condition', cols = cols_sel, pt.size = 0.0)
ggsave(filename = paste0(outDir, 'pseudotime_vs_realTime.pdf'), width = 10, height = 8)

saveRDS(aa, file = paste0(RdataDir, 
                           'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_reduction.DM_princurves_pseudotime_',
                           species, version.analysis, '.rds'))

##########################################
# plot heatmap  
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_reduction.DM_princurves_pseudotime_',
                           species, version.analysis, '.rds'))

candidates = readRDS(file = paste0(outDir, 'DElist_1961genes_pairwiseComaprison.rds'))
candidates = unique(candidates$gene)

DimPlot(aa, cols = cols_sel, group.by = 'condition', reduction = 'DC')
FeaturePlot(aa, features = candidates$gene[1])

source('plotting_utility.R')
plot_genes_branched_heatmap(seuratObj = aa, 
                            gene_subset = candidates,
                            nbCell_condition = 50)


##########################################
# plot gene examples 
##########################################
p0 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = FALSE, cols = cols_sel)
p1 = FeaturePlot(aa, features = 'Adgra2')

p0 + p1

ggsave(filename = paste0(outDir, 'geneExamples_Adgra2.pdf'), width = 20, height = 8)

