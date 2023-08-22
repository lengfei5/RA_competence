##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: to identify early branching genes
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Apr 12 15:48:45 2023
##########################################################################
##########################################################################
names(cols) = levels

levels_sels = c("day2_beforeRA", 
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/RA_symetryBreaking/branching_genes_BGP/')
system(paste0('mkdir -p ', outDir))

version.analysis = paste0('_R13547_10x_mNT_20220813', '_ES.beforeRA.and.RA')

library(slingshot, quietly = FALSE)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(scater)
library(SingleCellExperiment)
library(scran)
library(RColorBrewer)

##########################################
# prepare the dataset from seurat and d2.before.RA and RA samples
##########################################
Reselect_samples_to_include = FALSE
if(Reselect_samples_to_include){
    
  aa = readRDS(file = paste0(RdataDir, 
                             'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                             'cellCycleScoring_annot.v2_newUMAP_clusters_time_',
                             species, version.analysis, '.rds'))
  
  Idents(aa) = aa$condition
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)
  
  aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '12' &
                                               aa$celltypes != '4' &
                                               aa$celltypes != '2' &
                                               aa$celltypes != '13' &
                                               aa$celltypes != '15')])
  
  aa = subset(aa, cells = colnames(aa)[which(aa$condition != 'day0_beforeRA' &
                                                aa$condition != 'day1_beforeRA')])
  
  aa = subset(aa, cells = colnames(aa)[which(aa$condition != 'day6_RA')])
  aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != 'Neurons')])
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.1)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  #DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)
  
  aa <- FindNeighbors(aa, dims = 1:30)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  p1 + p2
  
  ggsave(filename = paste0(outDir, 
                           'UMAP_RAsymmetryBreaking.onlyday3rep1_timePoints_selected_forDM',
                           version.analysis, '.pdf'), 
         width = 18, height = 8)
  
  saveRDS(aa, file = paste0(outDir, 'RA_d2.beforeRA_no.day6_noNeurons.rds'))
  
  
}else{
  aa = readRDS(file = paste0(outDir, 'd2.beforeRA_RA_noNeurons_downsample_10.5k.rds'))
  
  aa = DietSeurat(aa, counts = TRUE, data = TRUE,
                  scale.data = FALSE,
                  features = NULL,
                  assays = 'RNA',
                  dimreducs = 'umap',
                  graphs = NULL,
                  misc = TRUE
  )
  
  saveRDS(aa, file = paste0(outDir, 'd2.beforeRA_RA_noNeurons_downsample_10.5k_singleAssay.RNA.rds'))
  
}


########################################################
########################################################
# Section I : monocle2 
# original code from 
# http://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle
########################################################
########################################################
Test_monocle2_BEAMtest = FALSE
if(Test_monocle2_BEAMtest){
  library(monocle) # monocle2, BEAM was discarded in monocle3
  
  aa = readRDS(file = paste0(outDir, 'd2.beforeRA_RA_noNeurons_downsample_10.5k_singleAssay.RNA.rds'))
  #monocle::importCDS(aa)
  Idents(aa) = aa$condition
  aa = subset(aa, downsample = 500)
  
  pd <- new("AnnotatedDataFrame", data = data.frame(aa@meta.data))
  fd <- new("AnnotatedDataFrame", data = data.frame(row.names = rownames(aa), genes = rownames(aa), 
                                                    gene_short_name = rownames(aa)))
  count_matrix = as.matrix(aa@assays$RNA@counts)
  HSMM <- newCellDataSet(aa@assays$RNA@counts,
                         phenoData = pd,
                         featureData = fd, 
                         expressionFamily=negbinomial.size())
  
  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)
  
  # Filtering low-quality cells
  HSMM <- detectGenes(HSMM, min_expr = 0.1)
  print(head(fData(HSMM)))
  
  expressed_genes <- row.names(subset(fData(HSMM),
                                      num_cells_expressed >= 10))
  
  print(head(pData(HSMM)))
  #HSMM <- detectGenes(HSMM, min_expr = 0.1)
  
  # Log-transform each value in the expression matrix.
  # L <- log(exprs(HSMM[expressed_genes,]))
  
  # Standardize each gene, so that they are all on the same scale,
  # Then melt the data with plyr so we can plot it easily
  # melted_dens_df <- reshape2::melt(Matrix::t(scale(Matrix::t(L))))
  
  # Plot the distribution of the standardized gene expression values.
  # qplot(value, geom = "density", data = melted_dens_df) +
  #   stat_function(fun = dnorm, size = 0.5, color = 'red') +
  #   xlab("Standardized log(FPKM)") +
  #   ylab("Density")
  
  ##########################################
  # Constructing Single Cell Trajectories 
  ##########################################
  HSMM_myo = HSMM
  rm(HSMM)
  require(tictoc)
  
  tic()
  diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                        fullModelFormulaStr = "~time", cores = 1)
  toc()
  
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  
  HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
  plot_ordering_genes(HSMM_myo)
  
  tic()
  HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,
                              #method = 'ICA'
                              method = 'DDRTree'
  )
  toc()
  
  
  HSMM_myo <- orderCells(HSMM_myo, num_paths = 2)
  
  plot_cell_trajectory(HSMM_myo, color_by = "condition")
  
  plot_cell_trajectory(HSMM_myo, color_by = "celltypes")
  
  
  saveRDS(HSMM_myo, file = paste0(outDir, 'd2.beforeRA_RA_noNeurons_downsample_2k_monocle2_v1.rds'))
  
  
  
  GM_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
      T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
      return(as.numeric(names(T0_counts)[which
                                         (T0_counts == max(T0_counts))]))
    } else {
      return (1)
    }
  }
  HSMM_myo <- orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))
  plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")
  
}


########################################################
########################################################
# Section II : test slingshot to get the trajectory and pseudotime 
# original code from 
# http://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle
########################################################
########################################################
aa = readRDS(file = paste0(outDir, 'RA_d2.beforeRA_no.day6_noNeurons.rds'))
Idents(aa) = aa$condition

aa = subset(aa, downsample = 1000)

##########################################
# Diffusion map for dimension reduction 
##########################################
#metadata = aa@meta.data
sce = as.SingleCellExperiment(aa)
#rm(aa)
dec <- modelGeneVar(sce)

#plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
#curve(metadata(dec)$trend(x), col="blue", add=TRUE)

for(nb_features in c(3000, 5000, 8000))
{
  for(n_neighbors in c(100, 200, 300, 500))
  {
    # nb_features = 5000;n_neighbors = 200
    top.hvgs <- getTopHVGs(dec, n=nb_features)
    
    sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 100)
    # reducedDimNames(sce)
    ll.pca = reducedDim(sce, 'PCA')[, c(1:50)]
    
    # tic()
    # dm <- DiffusionMap(ll.pca, sigma = 'local', k = n_neighbors, 
    #                    n_eigs = 50, distance = 'euclidean')
    # 
    # saveRDS(dm, file = paste0(RdataDir, 'diffusionMap_firstBifurcation_RA.vs.noRA_with.sce.PCA_euclidean_',
    #                           'localSignal_nfeatures.', nb_features, '_nbNeighbors.', n_neighbors,
    #                           '.rds'))
    # 
    #toc()
    
    tic()
    dm <- DiffusionMap(ll.pca, sigma = 'global', k = n_neighbors, n_eigs = 50, distance = 'euclidean')
    
    toc()
    
    cells = names(dm$DC1)
    metadata = aa@meta.data
    dcs = data.frame(DC1 = dm$DC1, DC2 = dm$DC2, DC3 = dm$DC3, 
                     DC4 = dm$DC4, DC5 = dm$DC5, stringsAsFactors = FALSE)
    dcs = dcs[match(rownames(metadata), cells), ]
    
    dcs = as.matrix(dcs)
    aa[["DC"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", assay = DefaultAssay(aa))
    
    rm(metadata)
    
    p1 = DimPlot(aa, reduction = 'DC', dims = c(1, 2),cols = cols_sel)
    p2 = DimPlot(aa, reduction = 'DC', dims = c(1, 3), cols = cols_sel)
    p3 = DimPlot(aa, reduction = 'DC', dims = c(1, 4), cols = cols_sel)
    p1 / p2 / p3
    
    DimPlot(aa, reduction = 'DC', dims = c(1, 4), group.by = 'celltypes')
    
    # library("plot3D")
    # library(rgl)
    # library(magick)
    # colors <- c("royalblue1", "darkcyan", "oldlace")
    # iris$color <- colors[ as.numeric( as.factor(iris$Species) ) ]
    # 
    # # Static chart
    # plot3d( iris[,1], iris[,2], iris[,3], col = iris$color, type = "s", radius = .2 )
    # play3d( spin3d( axis = c(0, 0, 1), rpm = 20), duration = 10 )
    # scatter3D(dcs[,1], dcs[,2], dcs[,3], theta = 150, phi = 20)
    
    saveRDS(aa, file = paste0(outDir, 'd2.beforeRA_RA_noNeurons_downsample_6k_DM.rds'))
    
  }
}


##########################################
# define clusters as slingshot and outlier detection because the principle curves are sensitive to the outliers
# manually correction of clusters 
##########################################
library(plotly)

aa = readRDS(file = paste0(outDir, 'd2.beforeRA_RA_noNeurons_downsample_6k_DM.rds'))

dcs = data.frame(aa[['DC']]@cell.embeddings[, c(1:5)])
dcs$condition = aa$condition

plot_ly(data.frame(dcs), x = ~DC_1, y = ~DC_2, z = ~DC_4, size = 3) %>%
  add_markers(color = ~ condition)


#data("slingshotExample")
#rd <- slingshotExample$rd
#cl <- slingshotExample$cl
#aa$celltypes = droplevels(aa$celltypes)

p1 = DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = FALSE)
p2 = DimPlot(aa, reduction = 'DC', dims = c(1, 4), group.by = 'celltypes', label = TRUE, repel = TRUE)
p1 + p2

aa$celltypes = as.character(aa$celltypes)
#aa$celltypes[which(aa$celltypes == '12'|aa$celltypes == '11'|aa$celltypes == '10')] = '0'
#DimPlot(aa, reduction = 'DC', dims = c(1, 4), group.by = 'celltypes', label = TRUE, repel = TRUE)

#DimPlot(aa, group.by = 'celltypes', label = TRUE, repel = TRUE)

rd =  aa[['DC']]@cell.embeddings[, c(1, 4)]
cl = aa$celltypes

# reclustering 
set.seed(2022)
res_kmean = kmeans(rd, centers = 10)
cl2 <- res_kmean$cluster

centers <- res_kmean$centers[res_kmean$cluster, ] 
distances <- sqrt(rowSums((rd - centers)^2))

hist(distances, breaks = 100)
abline(v = 0.005, col = 'red')
outliers <- which(distances>0.005)
cat(length(outliers), ' outliers found \n')


cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[res_kmean$cluster]

plot(rd, pch=19, col=cols, cex=0.2)
points(res_kmean$centers, col='darkred', pch=15, cex=1)
points(rd[outliers, ], pch="+", col = 'darkgray', cex=0.8)


# Process.data.for.slingshot = FALSE
# if(Process.data.for.slingshot){
#   sce = as.SingleCellExperiment(aa)
#   pca <- prcomp(t(assays(sce)$logcounts), scale. = FALSE)
#   rd1 <- pca$x[,1:2]
#   
#   plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
#   
#   rd1 = aa[['pca']]@cell.embeddings[, c(1:2)]
#   rd2 <- aa[['DC']]@cell.embeddings[, c(1:2)]
#   reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)
#   
# }


colData(sce)$kmeans <- cl2
# 
aa$dc_clusters = cl2
aa$dc_clusters[which(aa$celltypes == '14'| aa$celltypes == '11')] = NA
# 
# 
DimPlot(aa, reduction = 'DC', dims = c(1, 4), group.by = 'dc_clusters', label = TRUE, repel = TRUE)
# ggsave(filename = paste0(outDir, 'DMcomponents_clustering_for_trajectory_slingshot_outliers.pdf'), 
#        width = 16, height = 8)

table(aa$condition, aa$dc_clusters)

### manually change the cluster labels
cl2 = aa$dc_clusters
cl2[which(cl2== 9 | cl2 == 8 | cl2==5)] = 5 # day4_RA
cl2[which(cl2 == 10 | cl2 == 2)] = 10 # day4_noRA 
colData(sce)$kmeans <- cl2

aa$dc_clusters = cl2
#aa$dc_clusters[!is.na(aa$dc_clusters_outliers)] = NA
DimPlot(aa, reduction = 'DC', dims = c(1, 4), group.by = 'dc_clusters', label = TRUE, repel = TRUE)

ggsave(filename = paste0(outDir, 
                         'DMcomponents_clustering.manaulCorrection_for_trajectory_slingshot_outliers.pdf'), 
       width = 16, height = 8)


# default singshot 
#sce <- slingshot(sce, clusterLabels = 'kmeans', reducedDim = 'DC')
#summary(sce$slingPseudotime_1)
#colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
#plot(reducedDims(sce)$DC, col = colors[cut(sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
#lines(SlingshotDataSet(sce), lwd=2)

table(aa$condition, aa$dc_clusters)

if(trajectory_pseudotime_method == 'slingshot'){ #
  
  #kk_clean = which(!is.na(aa$dc_clusters))
  #rd3 = rd2[-outliers, ]
  #cl3 = cl2[-outliers]
  
  times = c('0' = 4,
            '1' = 3,
            '2' = 5, '10' = 5,
            '4' = 1, 
            '3' = 2, 
            'NP' = 6, 'FP' = 6
            ) 
  
  rd =  aa[['DC']]@cell.embeddings[, c(1,4)]
  cl2 = aa$dc_clusters
  
  index_keep = which(!is.na(cl2))
  rd = rd[index_keep, ]
  cl2 = cl2[index_keep]
  
  lin1 <- getLineages(rd, cl2, start.clus = '5', end.clus = c('1', '6')
                      )
  
  #cl.uniq = unique(cl)
  cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl2]
  #names(cols) = cl.uniq
  
  pdf(paste0(outDir, "Slingshot_getLineage_withStarting.cluster.pdf"),
      height = 6, width =10, useDingbats = FALSE)
  
  plot(rd, col = cols, asp = 1, pch = 1, cex = 0.5)
  lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')
  
  dev.off()
  
  crv1 <- getCurves(lin1, shrink = 1, stretch = 2,  smoother = "smooth.spline", thresh = 0.1, maxit = 50, 
                    shrink.method = 'density')
  # crv1 <- getCurves(lin1, shrink = TRUE, stretch = 2,  smoother = "loess",
  #                   thresh = 0.001, maxit = 50, 
  #                   shrink.method = 'density')
  
  
  pdf(paste0(outDir, "Slingshot_getCurve_withLineage.pdf"),
      height = 6, width =10, useDingbats = FALSE)
  
  
  #crv1
  plot(rd, col = cols, asp = 1, pch = 16, cex = 0.2)
  lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')
  
  dev.off()
  
  scrv = slingCurves(crv1)
  pst = as.data.frame(slingPseudotime(crv1))
  
  pst$weight.lineage1 = scrv$Lineage1$w
  pst$weight.lineage2 = scrv$Lineage2$w
  
  save(aa, sce, lin1, crv1, 
       file = paste0(outDir, 
                     'd2.beforeRA_RA_noNeurons_downsample.6k_DM_slingshot_lineage_pseudottime.Rdata'))
  
  ##########################################
  # prepare the input files for BGP 
  ##########################################
  sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
  sps = unique(sps$gene)
  
  tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
  tfs = unique(tfs$`HGNC symbol`)
  tfs = as.character(unlist(sapply(tfs, firstup)))
  
  pseudotime.scaling = function(X) {
    return((X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE)))
  }
  
  load(file = paste0(outDir, 
                     'd2.beforeRA_RA_noNeurons_downsample.6k_DM_slingshot_lineage_pseudottime.Rdata'))
  
  scrv = slingCurves(crv1)
  pst = data.frame(Lineage1 = scrv$Lineage1$lambda, Lineage2 = scrv$Lineage2$lambda,
                   weight.lineage1 = scrv$Lineage1$w,
                   weight.lineage2 = scrv$Lineage2$w)
  
  #pst$Lineage1 = pseudotime.scaling(pst$Lineage1)
  #pst$Lineage2 = pseudotime.scaling(pst$Lineage2)
  dpt = data.frame(aa@reductions$DC@cell.embeddings[, c(1,4)])
  dpt = dpt[match(rownames(pst), rownames(dpt)), ]
  dpt$state = 0
  dpt$pseudotime = 0
  dpt = dpt[, c(3,4,1,2)]
  
  dpt = data.frame(dpt, pst)
  
  dpt$dc_clusters = aa$dc_clusters[match(rownames(dpt), colnames(aa))] 
  dpt$state[which(dpt$weight.lineage1 == dpt$weight.lineage2)] = 1
  dpt$state[which(dpt$weight.lineage1 > dpt$weight.lineage2)] = 2
  dpt$state[which(dpt$weight.lineage1 < dpt$weight.lineage2)] = 3
  dpt$state[which(dpt$DC_1 <0.0011)] = 1
  
  ggplot(data.frame(dpt), aes(DC_1, DC_4, color= state))  + 
    geom_point(size=4) 
  #scale_shape_manual(values=1:nlevels(pca2save$time)) +
  #geom_text(hjust = 0.5, nudge_y = 0.5, size=2.5)
  
  jj = which(dpt$state == 1)
  dpt$pseudotime[jj] = apply(dpt[jj, c(5,6)], 1, mean) 
  
  jj = which(dpt$state == 2)
  dpt$pseudotime[jj] = dpt[jj, c(5)]
  
  jj = which(dpt$state == 3)
  dpt$pseudotime[jj] = dpt[jj, c(6)]
  
  dpt$pseudotime = pseudotime.scaling(dpt$pseudotime)
  
  aa$dpt = dpt$pseudotime[match(colnames(aa), rownames(dpt))]
  
  DimPlot(aa, reduction = 'DC', dims = c(1, 4), group.by = 'dc_clusters', label = TRUE, repel = TRUE)
  FeaturePlot(aa, features = 'dpt')
  
  
  sbdata = t(as.matrix(aa@assays$RNA@data))
  sbdata = sbdata[match(rownames(dpt), rownames(sbdata)), ]
  
  ss = apply(sbdata, 2, mean)
  ss2 = apply(sbdata, 2, function(x) length(which(x>0)))
  
  plot(log10(ss), log10(ss2))
  hist(ss2/ncol(sbdata), breaks = 100)
  abline(v = 0.025)
  
  min_pct = 0.025
  cel_pct = ss2/ncol(sbdata)
  geneIndex_sels = which(cel_pct > min_pct)
  cat(length(geneIndex_sels), ' genes selected \n')
  
  sbdata = sbdata[, geneIndex_sels]
  
  ## save the cell annotation
  write.csv(dpt[, c(1:4)], file = paste0(outDir, 'symmetry_breaking_dpt.csv'), 
            row.names = TRUE, quote = FALSE)
  
  # only TFs 
  mm = match(colnames(sbdata), c(tfs))
  sbdata = sbdata[, which(!is.na(mm))]
  
  ss1 = apply(sbdata, 2, function(x){length(which(x>0))})
  ss2 = apply(sbdata, 2, mean)
  
  sbdata = sbdata[, order(-ss1)]
  
  write.csv(sbdata, file = paste0(outDir, 'symmetry_breaking_scRNAseq_data_tfs.csv'), 
            row.names = TRUE, quote = FALSE)
  
  # only sps
  mm = match(colnames(sbdata), c(sps))
  sbdata = sbdata[, which(!is.na(mm))]
  write.csv(sbdata, file = paste0(outDir, 'symmetry_breaking_scRNAseq_data_sps.csv'), 
            row.names = TRUE, quote = FALSE)
  
  
  
}

########################################################
########################################################
# Section : summary the BGP results
# 
########################################################
########################################################
## TFs
file_list = list.files(path = paste0(outDir, 'BGP_out_tfs'), pattern =  '*.txt', full.names = TRUE)
file_name = basename(file_list)

## signaling molecules
file_list2 = list.files(path = paste0(outDir, 'BGP_out_sps'), pattern =  '*.txt', full.names = TRUE)
file_name2 = basename(file_list2)
mm = match(file_name2, file_name)
file_list2 = file_list2[which(is.na(mm))]
file_name2 = file_name2[which(is.na(mm))]

file_list = c(file_list, file_list2)
file_name = basename(file_list)

res = data.frame(gene = rep(NA, length(file_list)), 
                 branching.time = rep(NA, length(file_list)),
                 bf = rep(NA, length(file_list)))

for(n in 1:length(file_list))
{
  # n = 1
  cat(n, '--', basename(file_list[n]), '\n')
  test = read.table(file_list[n], sep = '\t', header = FALSE)
  res$gene[n] = test[1, 1]
  res$branching.time[n] = test[1, 2]
  res$bf[n] = test[1, 3]
  
}

res = res[order(res$branching.time), ]

saveRDS(res, file = paste0(RdataDir, 'symmetry_breaking_early.tfs.sps_BGP_output.rds'))

##########################################
# visualize the branching time
##########################################
library(ggplot2)
require("ggrepel")

#res = readRDS(file = paste0(RdataDir, 'symmetry_breaking_early.tfs.sps_BGP_output.rds'))
res = readRDS(file = paste0(RdataDir, 'symmetry_breaking_early.tfs.sps_BGP_output.rds'))
res = res[which(res$bf>0), ]


#res$bf = log10(res$bf)
res %>%
  as_tibble() %>%
  ggplot(aes(x = branching.time, y = bf, label = gene)) +
  geom_point(size=1.5) + 
  #geom_text(size = 3) +
  theme_classic() + 
  geom_text_repel(label.size = 0.2, box.padding = unit(0.2, "lines"), point.padding = unit(0.2, "lines"))

ggsave(filename = paste0(outDir, 'TFs_SPs_branchingtime_bayesianFactor',
                         version.analysis, '.pdf'), 
       width = 8, height =6)

library(forcats)
# Reorder following the value of another column:
res %>%
  mutate(name = fct_reorder(gene, dplyr::desc(branching.time))) %>%
  top_n(n = 300) %>%
  ggplot( aes(x=name, y=branching.time)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  theme_bw()

ggsave(filename = paste0(outDir, 'TF_branchingtime_top30',
                         version.analysis, '.pdf'), 
       width = 10, height = 6)

##########################################
# test the umap with branching genes 
##########################################
genes_bgp = res$gene[which(res$bf > 5)]

sub_obj = subset(aa, features = genes_bgp)
sub_obj <- RunPCA(sub_obj, features = rownames(sub_obj), verbose = FALSE, weight.by.var = FALSE, 
                  npcs = 30)
ElbowPlot(sub_obj, ndims = 30)

sub_obj <- RunUMAP(sub_obj, dims = 1:10, n.neighbors = 30, min.dist = 0.1)
DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

FeaturePlot(sub_obj, features = genes_bgp[c(1:20)])

p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) +
  NoLegend()
p1 =  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)+
  NoLegend()

p0 + p1

