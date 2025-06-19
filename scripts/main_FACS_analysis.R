##########################################################################
##########################################################################
# Project: RA competence project
# Script purpose: analyze FACS 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Sat Apr 23 23:09:17 2022
##########################################################################
##########################################################################
rm(list = ls())

resDir = '../results/FACS_analysis_clusteringWT'
RdataDir = paste0(resDir, '/Rdata')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

#dataDir = '/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/FACS_6h_secondrun'
dataDir = '/groups/tanaka/People/current/jiwang/projects/RA_competence/FACS_timecourse/'

library(CATALYST)
library(dplyr)
library(flowCore)
library(flowWorkspace)
library(ggcyto)
library(ggplot2)
library(mvtnorm)
library(openCyto)
library(scDataviz)
library(cowplot)
library(SingleCellExperiment)


########################################################
########################################################
# Section I: data collectio nand metadata
# 
########################################################
########################################################
files = list.files(path = dataDir, pattern = '*.csv', full.names = TRUE)

for(n in 1:length(files))
{
  # n = 1
  
  test = read.csv(file = files[n], header = TRUE)
  names = basename(files[n])
  
  cat(n, '--', names, '\n')
  
  names = gsub('  #', '_', names)
  names = unlist(strsplit(as.character(names), '_'))
  if(n == 1){
    metadata = data.frame(condition = rep(names[2], nrow(test)), 
                          time = rep(paste0(names[4], 'h'), nrow(test)), 
                          rep = rep(names[3], nrow(test)))
    aa = test
  }else{
    metadata = rbind(metadata, 
                     data.frame(condition = rep(names[2], nrow(test)), 
                                time = rep(paste0(names[4], 'h')), 
                                rep = rep(names[3], nrow(test))))
    aa = rbind(aa, test)
  }
  
}

save(aa, metadata, file = paste0(RdataDir, '/facs_wt_data_meta.Rdata'))

########################################################
########################################################
# Section II: data processing, QCs and dimension reduction
# 
########################################################
########################################################
# import our own data
load(paste0(RdataDir, '/facs_wt_data_meta.Rdata')) # aa and metadata

colnames(aa) = gsub('FJComp.', '', colnames(aa))
colnames(aa) = gsub('.H', '', colnames(aa))
mat = aa
rm(aa)

saveRDS(metadata, file = paste0(RdataDir, '/metadata_cells_CyTOF.rds'))

save(mat, metadata, file = paste0(RdataDir, '/cytof_mat_metadata.Rdata'))

##########################################
# test Robinson's analysis 
# https://www.bioconductor.org/packages/release/workflows/vignettes/cytofWorkflow/inst/doc/cytofWorkflow.html
# #data-organization
# https://f1000research.com/articles/6-748/v3
##########################################


##########################################
# normalize data
# original code from https://www.biostars.org/p/393974/
# or https://github.com/kevinblighe/scDataviz
##########################################
Use.scDataviz.analysis = FALSE
if(Use.scDataviz.analysis){
  # After reading in the data, you need to normalise it and eliminate junk cells. 
  # As you know, CyTOF data is typically normalised by hyperbolic arc-sine with a factor of 5.
  
  # exmaple code from https://github.com/kevinblighe/scDataviz 
  
  load(file = paste0(RdataDir, '/cytof_mat_metadata.Rdata'))  
  
  ##########################################
  # normalization and transformation of fact data matrix
  ##########################################
  # Set background noise threshold - values below this are set to 0
  BackgroundNoiseThreshold <- 1
  
  # Euclidean norm threshold - this is the square root of the sum of all the squares
  EuclideanNormThreshold <- 1
  
  # Choose a transformation function (any mathematical function)
  transFun <- function (x) asinh(x)
  
  # Set hyperbolic arc-sine factor (NB - asinh(x/5) is recommended for CyTOF and FACS data)
  asinhFactor <- 5
  
  x = as.matrix(mat)
  
  x <- x[apply(x, 1, FUN=function(x) sqrt(sum(x^2)))>EuclideanNormThreshold,]
  NoiseCorrected <- x
  NoiseCorrected[NoiseCorrected<BackgroundNoiseThreshold] <- 0
  x <- transFun(NoiseCorrected/asinhFactor)
  
  ##########################################
  # prepare the metadata 
  ##########################################
  metadata = data.frame(cells = paste0('cell_', c(1:nrow(metadata))), metadata, stringsAsFactors = FALSE) 
  rownames(metadata) = metadata$cells
  rownames(x) = rownames(metadata)
  rownames(mat) = rownames(metadata)
  
  metadata$treatment = metadata$condition
  metadata$time2 = metadata$time
  metadata$time[which(metadata$time == '48h')] = 'd2'
  metadata$time[which(metadata$time == '54h')] = 'd2.6h'
  metadata$time[which(metadata$time == '60h')] = 'd2.12h'
  metadata$time[which(metadata$time == '66h')] = 'd2.18h'
  metadata$time[which(metadata$time == '72h')] = 'd3'
  metadata$time[which(metadata$time == '78h')] = 'd3.6h'
  metadata$time[which(metadata$time == '84h')] = 'd3.12h'
  metadata$time[which(metadata$time == '90h')] = 'd3.18h'
  metadata$time[which(metadata$time == '96h')] = 'd4'
  metadata$condition = paste0(metadata$treatment, '_', metadata$time)
  #metadata$condition = gsub('noRA_d2', "beforeRA_d2", metadata$condition)
  
  cc.levels = c("noRA_d2", 
                "noRA_d2.6h", "noRA_d2.12h", "noRA_d2.18h", "noRA_d3",
                "noRA_d3.6h", "noRA_d3.12h", "noRA_d3.18h", "noRA_d4", 
                "RA_d2.6h", "RA_d2.12h", "RA_d2.18h", "RA_d3",
                "RA_d3.6h", "RA_d3.12h", "RA_d3.18h", "RA_d4")
  
  metadata$condition = factor(metadata$condition, levels = cc.levels)
  
  save(metadata, x, mat, file = paste0(RdataDir, '/cytof_mat_transformedData_metadata.Rdata'))  
  
  ###################
  # make singleCellExperiment object
  load(file = paste0(RdataDir, '/cytof_mat_transformedData_metadata.Rdata'))
  #subsample = sample(c(1:nrow(x)), size = 15000, replace = FALSE)
  subsample = c(1:nrow(x))
  
  sce <- importData(x[subsample, ],
                    assayname = 'normcounts',
                    metadata = metadata[subsample, ])
  sce
  
  ## PCA analysis
  library(PCAtools)
  assayNames(sce)
  
  p <- pca(assay(sce, 'normcounts'), metadata = metadata(sce))
  
  p1 = biplot(p,
              x = 'PC1', y = 'PC2',
              lab = NULL,
              xlim = c(min(p$rotated[,'PC1'])-1, max(p$rotated[,'PC1'])+1),
              ylim = c(min(p$rotated[,'PC2'])-1, max(p$rotated[,'PC2'])+1),
              pointSize = 1.0,
              colby = 'condition',
              #shape = 'time',
              legendPosition = 'right',
              title = 'PCA applied to CyTOF data',
              caption = paste0('10000 cells randomly selected after ',
                               'having filtered for low variance'))
  
  p2 = biplot(p,
              x = 'PC3', y = 'PC4',
              lab = NULL,
              xlim = c(min(p$rotated[,'PC2'])-1, max(p$rotated[,'PC2'])+1),
              ylim = c(min(p$rotated[,'PC3'])-1, max(p$rotated[,'PC3'])+1),
              pointSize = 1.0,
              colby = 'condition',
              #shape = 'time',
              legendPosition = 'right',
              title = 'PCA applied to CyTOF data',
              caption = paste0('10000 cells randomly selected after ',
                               'having filtered for low variance'))
  
  p1 + p2
  
  reducedDim(sce, 'PCA') <- p$rotated
  
  config <- umap::umap.defaults
  config$min_dist <- 0.01
  config$n_neighbors = 10
  config$metric = "euclidean"
  
  library(tictoc)
  tic()
  sce = performUMAP(sce, assay = 'normcounts', config = config)
  toc()
  
  saveRDS(sce, file = paste0(RdataDir, '/cytof_mat_transformedData_metadata_scDataviz_umap.rds'))
  
  # sce <- performUMAP(sce, reducedDim = 'PCA', dims = c(1:5))
  
  # metadataPlot(sce,
  #              colby = 'condition',
  #              #colkey = c(Healthy = 'royalblue', Disease = 'red2'),
  #              title = 'Disease status',
  #              subtitle = 'UMAP performed on expression values',
  #              legendLabSize = 16,
  #              axisLabSize = 20,
  #              titleLabSize = 20,
  #              subtitleLabSize = 16,
  #              captionLabSize = 16)
  
  p1 = metadataPlot(sce,
                    colby = 'condition',
                    #colkey = c(Healthy = 'royalblue', Disease = 'red2'),
                    title = 'Disease status',
                    subtitle = 'UMAP performed on expression values',
                    legendLabSize = 16,
                    axisLabSize = 20,
                    titleLabSize = 20,
                    subtitleLabSize = 16,
                    captionLabSize = 16)
  
  p2 = metadataPlot(sce,
                    colby = 'time',
                    title = 'Treatment type',
                    subtitle = 'UMAP performed on expression values',
                    legendLabSize = 16,
                    axisLabSize = 20,
                    titleLabSize = 20,
                    subtitleLabSize = 16,
                    captionLabSize = 16)
  
  p1 + p2
  
  ggout1 <- contourPlot(sce,
                        reducedDim = 'UMAP',
                        bins = 150,
                        subtitle = 'UMAP performed on expression values',
                        legendLabSize = 18,
                        axisLabSize = 22,
                        titleLabSize = 22,
                        subtitleLabSize = 18,
                        captionLabSize = 18)
  
  
}


##########################################
# Test again the Robinson's workflow 
# original code from:
# https://www.bioconductor.org/packages/release/workflows/vignettes/cytofWorkflow/inst/doc/
# cytofWorkflow.html#data-transformation
##########################################
Use.Robinson.workflow = FALSE
if(Use.Robinson.workflow){
  if(Test_example){
    library(readxl)
    url <- "https://zenodo.org/records/10039274/files"
    md <- "PBMC8_metadata.xlsx"
    download.file(file.path(url, md), destfile = md, mode = "wb")
    md <- read_excel(md)
    head(data.frame(md))
    
    library(HDCytoData)
    fs <- Bodenmiller_BCR_XL_flowSet()
    
    panel <- "PBMC8_panel_v3.xlsx"
    download.file(file.path(url, panel), destfile = panel, mode = "wb")
    panel <- read_excel(panel)
    head(data.frame(panel))
    
    all(panel$fcs_colname %in% colnames(fs))
    
    # specify levels for conditions & sample IDs to assure desired ordering
    # specify levels for conditions & sample IDs to assure desired ordering
    md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))
    md$sample_id <- factor(md$sample_id, 
                           levels = md$sample_id[order(md$condition)])
    
    # construct SingleCellExperiment
    sce <- prepData(fs, panel, md, features = panel$fcs_colname)
    
    p <- plotExprs(sce, color_by = "condition")
    p$facet$params$ncol <- 6
    p
    
    n_cells(sce) 
    
    plotCounts(sce, group_by = "sample_id", color_by = "condition")
    
    pbMDS(sce, color_by = "condition", label_by = "sample_id")
    
    plotExprHeatmap(sce, scale = "last",
                    hm_pal = rev(hcl.colors(10, "YlGnBu")))
    
    # run t-SNE/UMAP on at most 500/1000 cells per sample
    set.seed(1234)
    sce <- runDR(sce, "TSNE", cells = 500, features = "type")
    sce <- runDR(sce, "UMAP", cells = 1e3, features = "type")
    
    plotDR(sce, "UMAP", color_by = "sample_id")
    
  }
  
  ############
  ## start with the same table as Meri
  ############
  mat = read.csv2(file = '../data/data_metadata_clusterIDs_perCondition_nbClusters_7.csv', 
                 header = TRUE, row.names = c(1))
  
  kk = which(mat$condition == 'noRA_d2')
  mat$treatment[kk] = 'beforeRA'
  mat$condition[kk] = 'beforeRA_d2'
  
  kk = which(mat$treatment != 'noRA')
  mat = mat[kk, ]
  
  counts = as.matrix(t(mat[, c(1:5)]))
  metadata = mat[, -c(1:5)]
  metadata$marker_class = 'type'
  
  #load(file = paste0(RdataDir, '/cytof_mat_metadata.Rdata'))
  
  
  #counts = as.matrix(t(mat))
  #colnames(counts) = paste0('cell_', c(1:ncol(counts)))
  #rownames(metadata) = colnames(counts)
  metadata = data.frame(metadata, stringsAsFactors = FALSE)
  
  sce <- SingleCellExperiment(assays=list(counts=counts),
                              colData=metadata, 
                              metadata = metadata)
  
  # here, x = input SCE with raw data
  # stored in assay "counts"
  cf <- 5 
  y <- assay(sce, "counts")
  y <- asinh(sweep(y, 1, cf, "/"))
  assay(sce, "exprs", FALSE) <- y
  
  # res <- normCytof(sce, beads = "dvs", k = 50, 
  #                  assays = c("counts", "exprs"), overwrite = FALSE)
  # # check number & percentage of bead / removed events
  # n <- ncol(sce); ns <- c(ncol(res$beads), ncol(res$removed))
  # data.frame(
  #   check.names = FALSE, 
  #   "#" = c(ns[1], ns[2]), 
  #   "%" = 100*c(ns[1]/n, ns[2]/n),
  #   row.names = c("beads", "removed"))
  # 
  # sce <- res$data
  # assayNames(sce)
  
  # load(file = paste0(RdataDir, '/cytof_mat_metadata.Rdata'))
  # metadata$time = gsub('_','.', metadata$time)
  # metadata$condition = paste0(metadata$treatment, '_', metadata$time)
  # #metadata$condition = gsub('beforeRA', 'RA', metadata$condition)
  # #metadata$treatment = gsub('beforeRA', 'RA', metadata$treatment)
  # 
  # # Set background noise threshold - values below this are set to 0
  # BackgroundNoiseThreshold <- 1
  # 
  # # Euclidean norm threshold - this is the square root of the sum of all the squares
  # EuclideanNormThreshold <- 1
  # 
  # # Choose a transformation function (any mathematical function)
  # transFun <- function (x) asinh(x)
  # 
  # Set hyperbolic arc-sine factor (NB - asinh(x/5) is recommended for CyTOF and FACS data)
  # asinhFactor <- 5
  # 
  # x = as.matrix(mat)
  # 
  # x <- x[apply(x, 1, FUN=function(x) sqrt(sum(x^2)))>EuclideanNormThreshold,]
  # NoiseCorrected <- x
  # NoiseCorrected[NoiseCorrected<BackgroundNoiseThreshold] <- 0
  # x <- transFun(NoiseCorrected/asinhFactor)
  # 
  # sce <- importData(x,
  #                   assayname = 'exprs',
  #                   metadata = metadata)
  
  sce$sample_id = as.character(sce$condition)
  #sce$condition = gsub('noRA_d2', "beforeRA_d2", sce$condition)
  
  cc.levels = c("beforeRA_d2",
                "RA_d2.6h", "RA_d2.12h", "RA_d2.18h", "RA_d3",
                "RA_d3.6h", "RA_d3.12h", "RA_d3.18h", "RA_d4")
  
  sce$condition = factor(sce$condition, levels = cc.levels)
  sce@metadata$condition = sce$condition
  
  rowData(sce)$marker_name = rownames(sce)
  rowData(sce)$channel_name = NULL
  rowData(sce)$marker_class = 'type'
  
  
  p <- plotExprs(sce, color_by = "condition")
  p$facet$params$ncol <- 2
  p
  
  n_cells(sce) 
  
  plotCounts(sce, group_by = "sample_id", color_by = "condition")
  
  pbMDS(sce, color_by = "condition", label_by = "sample_id")
  
  plotExprHeatmap(sce, scale = "last",
                  hm_pal = rev(hcl.colors(10, "YlGnBu")))
  
  # run t-SNE/UMAP on at most 500/1000 cells per sample
  set.seed(1234)
  #sce <- runDR(sce, "TSNE", cells = 500, features = "type")
  
  sce <- runDR(sce, "PCA", ncomponents = 4,
               cells = NULL, features = "type")
  
  p1 = plotDR(sce, "PCA", color_by = "condition")
  p1
  
  pcs = reducedDim(sce, 'PCA')
  
  ### filter outliers cells
  sels = which(pcs[ ,1] < 2 & pcs[ ,2] <2)
  sce = sce[, sels]
  
  p1 = plotDR(sce, "PCA", color_by = "condition")
  p1
  
  cf <- 5 
  y <- assay(sce, "counts")
  #y <- asinh(sweep(y, 1, cf, "/"))
  assay(sce, "exprs", FALSE) <- y
  
  
  #p2 = plotDR(sce, "PCA", color_by = "FoxA2")
  #p3 = plotDR(sce, 'PCA', color_by = 'Pax6')
  #p2 + p3
  
  
  sce <- runDR(sce, "DiffusionMap", cells = 2000,
               features = "type")
  
  p1 = plotDR(sce, "DiffusionMap", color_by = "condition")
  p1
    
  p2 = plotDR(sce, "DiffusionMap", color_by = "FoxA2")
  p3 = plotDR(sce, 'DiffusionMap', color_by = 'Pax6')
  
  p2 + p3
  
  set.seed(1234)
  sce <- runDR(sce, "UMAP", cells = 3000, 
               features = "type",
               n_neighbors = 50, scale = TRUE,
               min_dist = 0.05, metric = "cosine"
               )
  
  figureDir = "../results/figures_tables_R13547_10x_mNT_20240522/"
  
  p1 = plotDR(sce, "UMAP", color_by = "condition")
  p1
  
  ggsave(paste0(figureDir, 'FACS_RA_umap.pdf'), width=8, height = 6) 
  
  p2 = plotDR(sce, "UMAP", color_by = "FoxA2", scale = TRUE)
  p3 = plotDR(sce, 'UMAP', color_by = 'Pax6', scale = TRUE)
  p4 = plotDR(sce, 'UMAP', color_by = 'Oct4', scale = TRUE)
  p5 = plotDR(sce, 'UMAP', color_by = 'Sox1', scale = TRUE)
  
  (p2 + p3)/(p4 + p5)
  
  ggsave(paste0(figureDir, 'FACS_RA_umap_genes.pdf'), width=16, height = 12) 
  
  
  saveRDS(sce, file = paste0(RdataDir, '/sce_FACS_RA_umap_saved.rds'))
  
    
}


########################################################
########################################################
# Section III: fitting the Gaussian Mixture Model for cluster 
# 
########################################################
########################################################
library(mclust)
load(file = paste0(RdataDir, '/cytof_mat_transformedData_metadata.Rdata'))

outDir = paste0(resDir, '/pooling_treatment_time_inclSOX2_v2.8')
if(!dir.exists(outDir)) dir.create(outDir)

##########################################
# pooling all treatment and time
##########################################
#subsample = sample(c(1:nrow(mat)), size = 10000, replace = FALSE)
subsample = c(1:nrow(mat))
#cat('time point -- ', t, '\n')
#subsample = which(metadata$time == t)

# not considering Sox2, not informative
mat = mat[subsample, c(5, 2, 4, 3, 1)]
metadata = metadata[subsample, ]

print(dim(mat))
print(dim(metadata))

# clPairs(mat, metadata$condition)
# BIC <- mclustBIC(mat)
# plot(BIC)

for(nb_clusters in c(3:8))
{
  # nb_clusters = 7
  cat('nb of cluster -- ', nb_clusters, '\n')
  
  Search_for_optimal_initiation = FALSE
  if(Search_for_optimal_initiation){
    logliks = c()
    
    for(n in 0:25)
    {
      
      set.seed(2000)
      mb = Mclust(mat, G = nb_clusters)
      
      # optimal selected model
      #mb$modelName
      
      # optimal number of cluster
      #mb$G
      cat(n, "--", mb$loglik, "\n")
      logliks = c(logliks, mb$loglik)
      
    }
    
    xx = data.frame(seeds = c(0:11, 0:25), loglik = logliks)
    saveRDS(xx, file = paste0(outDir, '/seed_loglikelihood_saved_v2.rds'))
    
    xx = readRDS(file = paste0("../results/FACS_analysis_clusteringWT/",
                               "pooling_treatment_time_inclSOX2_v2.6_testInitiation/",
                               "seed_loglikelihood_saved_v2.rds"))
    
  }
  
  set.seed(1000)
  mb = Mclust(mat, G = nb_clusters, control = emControl(itmax=500, tol = 1.e-6))
  
  cat('loglike --', mb$loglik, "\n")
  
  # probality for an observation to be in a given cluster
  #head(mb$z)
  
  # get probabilities, means, variances
  #summary(mb, parameters = TRUE)
  
  saveRDS(mb, file = paste0(outDir, '/res_mclust_nbClusters.', nb_clusters, '.rds'))
  
  clusters = mb$classification
  clusters = clusters[match(rownames(mat), names(clusters))]  
  
  keep = table(metadata$condition, mb$classification)
  
  manual_modify_clusterIndex = FALSE
  
  if(manual_modify_clusterIndex){
    
    # 3 > 7
    #index_map = c(4, 5, 2, 6, 7, 1, 3)
    index_map = c(1, 3, 4, 7, 6, 2, 5)
    
    xx = keep[, index_map]
    colnames(xx) = c(1:7)
    
    newclusters = clusters
    for(m in 1:length(index_map))
    {
      newclusters[which(clusters == index_map[m])] = m
    }
    
    keep = xx
    clusters = newclusters
    
  }
  
  
  #Compare amount of the data within each cluster
  write.csv2(keep, file = paste0(outDir, '/cellNumbers_perCluster_perCondition_nbClusters_', 
                                 nb_clusters, '.csv'))
  
  for(n in 1:nrow(keep)){
    keep[n, ] = keep[n, ]/sum(keep[n,])
  }
  
  write.csv2(keep, file = paste0(outDir, '/cellProportions_perCluster_perCondition_nbCluste_', 
                                 nb_clusters, '.csv'))
  
  res = data.frame(mat, metadata[match(rownames(mat), rownames(metadata)), ], stringsAsFactors = FALSE)
  
  res = data.frame(res, clusters, stringsAsFactors = FALSE)
  
  table(res$condition, res$clusters)
  write.csv2(res, file = paste0(outDir, '/data_metadata_clusterIDs_perCondition_nbClusters_', 
                                nb_clusters, '.csv'))
  
  #metadata$cluster = mb$classification
  cc = unique(clusters)
  cc = cc[order(cc)]
  
  pdf(paste0(outDir, "/markerIntensity_incl.SOX2_nbClusters_", nb_clusters, ".pdf"), 
      height = 3*nb_clusters, width =16, useDingbats = FALSE)
  
  attach(mtcars)
  par(mfrow=c(length(cc), ncol(mat))) 
  for(n in 1:length(cc))
  {
    c = cc[n];
    for(m in 1:ncol(mat))
    {
      # c = 1; m = 1;
      hist(mat[which(clusters == c), m], breaks = 50, xlim = range(mat),
           xlab = '', ylab = paste0('cluster_', c), main = colnames(mat)[m],
           col = n);
      
    }
  }
  
  dev.off()
  
  
  pdf(paste0(outDir, "/clusterProjection_incl.SOX2_nbClusters_", nb_clusters, ".pdf"), 
      height = 12, width =16, useDingbats = FALSE)
  
  #After the data is fit into the model, we plot the model based on clustering results.
  # plot(mb, "density")
  source('functions_plotMclust.R')
  plot.Mclust_cutomized(mb, what=c("classification"), cex = 0.01, 
                        addEllipses = TRUE, cex_clusterlabels = 2.0)
  
  
  #plot.surface_customized(mb)
  dev.off()
  
  
}

##########################################
# split treatment or time
##########################################
#subsample = which(metadata$treatment == 'RA')
times = unique(metadata$time)
for(t in times)
{
  load(file = paste0(RdataDir, '/cytof_mat_transformedData_metadata.Rdata'))
  
  cat('time point -- ', t, '\n')
  subsample = which(metadata$time == t)
  outDir = paste0(resDir, '/split_timepoints/time_', t, '_poolingTreatment')
  if(!dir.exists(outDir)) dir.create(outDir)
  
  # not considering Sox2, not informative
  mat = mat[subsample, c(2, 4, 3, 1)]
  metadata = metadata[subsample, ]
  
  print(dim(mat))
  print(dim(metadata))
  
  # clPairs(mat, metadata$condition)
  # BIC <- mclustBIC(mat)
  # plot(BIC)
  
  for(nb_clusters in c(3:8))
  {
    # nb_clusters = 3
    cat('nb of cluster -- ', nb_clusters, '\n')
    mb = Mclust(mat, G = nb_clusters)
    
    # optimal selected model
    mb$modelName
    
    # optimal number of cluster
    mb$G
    
    # probality for an observation to be in a given cluster
    head(mb$z)
    
    # get probabilities, means, variances
    summary(mb, parameters = TRUE)
    
    #Compare amount of the data within each cluster
    keep = table(metadata$condition, mb$classification)
    write.csv2(keep, file = paste0(outDir, '/cellNumbers_perCluster_perCondition_nbClusters_', nb_clusters, '.csv'))
    
    metadata$cluster = mb$classification
    cc = unique(metadata$cluster)
    cc = cc[order(cc)]
    
    pdf(paste0(outDir, "/markerIntensity_noSOX2_nbClusters_", nb_clusters, ".pdf"), 
        height = 3*nb_clusters, width =16, useDingbats = FALSE)
    
    attach(mtcars)
    par(mfrow=c(length(cc), ncol(mat))) 
    for(n in 1:length(cc))
    {
      c = cc[n];
      for(m in 1:ncol(mat))
      {
        # c = 1; m = 1;
        hist(mat[which(metadata$cluster == c), m], breaks = 50, xlim = range(mat),
             xlab = '', ylab = paste0('cluster_', c), main = colnames(mat)[m],
             col = n);
        
      }
    }
    
    dev.off()
    
    pdf(paste0(outDir, "/clusterProjection_noSOX2_nbClusters_", nb_clusters, ".pdf"), 
        height = 12, width =16, useDingbats = FALSE)
    
    #After the data is fit into the model, we plot the model based on clustering results.
    # plot(mb, "density")
    source('functions_plotMclust.R')
    plot.Mclust_cutomized(mb, what=c("classification"), cex = 0.01, 
                          addEllipses = TRUE, cex_clusterlabels = 2.0)
    
    
    #plot.surface_customized(mb)
    dev.off()
    
    
  }
}

