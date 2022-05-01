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
resDir = '../results/FACS_analysis'
RdataDir = paste0(resDir, '/Rdata')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/groups/tanaka/People/current/jiwang/projects/RA_competence/FACS_timecourse/'

########################################################
########################################################
# Section : metadata and data processing 
# 
########################################################
########################################################
metadata = read.delim(file = paste0(dataDir, 'metadata.txt'), sep = '\t', header = FALSE)

colnames(metadata) = 'fileName'
metadata$fileName = gsub('.fcs', '', metadata$fileName)
metadata$condition = sapply(metadata$fileName, function(x) unlist(strsplit(as.character(x), '[#]'))[1])
metadata$sample = sapply(metadata$fileName, function(x) unlist(strsplit(as.character(x), '[#]'))[2])
metadata$sampleID = sapply(metadata$sample, function(x) unlist(strsplit(as.character(x), '_'))[4])

metadata$time = metadata$condition
metadata = metadata[, c(4, 2, 5, 1, 3)]

kk = grep('_noRA', metadata$condition)
metadata$condition[kk] = 'noRA'

kk = grep('_RA', metadata$condition)
metadata$condition[kk] = 'RA'

kk = grep('d2', metadata$condition)
metadata$condition[kk] = 'beforeRA'

metadata$sampleID = gsub('Sample', '', metadata$sampleID)

metadata$time = gsub('noRA_', '', metadata$time)
metadata$time = gsub('RA_', '', metadata$time)

metadata$time = gsub('_noRA', '', metadata$time)
metadata$time = gsub('_RA', '', metadata$time)

metadata = metadata[order(as.integer(metadata$sampleID)), ]

saveRDS(metadata, file = paste0(RdataDir, '/metadata.rds'))

##########################################
# data processing 
##########################################
aa = read.csv(file = paste0(dataDir, 'AIgood_AllSingleLive_AllParam_1_ScaleVal.csv'), header = TRUE)

jj = grep('Comp.*.H', colnames(aa))
aa = aa[, c(1, jj, 2:5, 31:32)]

saveRDS(aa, file = paste0(RdataDir, '/facs_data.rds'))

########################################################
########################################################
# Section : data processing and QCs
# 
########################################################
########################################################
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

# import our own data
aa =  readRDS(file = paste0(RdataDir, '/facs_data.rds'))
rownames(aa) = paste0('cell_', c(1:nrow(aa)))

mat = aa[, c(2:6)]
metadata = aa[, -c(2:6)]
colnames(mat) = gsub('.H', '', gsub('Comp.','', colnames(mat)))

meta2 = readRDS(file = paste0(RdataDir, '/metadata.rds'))
metadata = data.frame(metadata, meta2[match(metadata$SampleID, meta2$sampleID), ], stringsAsFactors = FALSE)
metadata$treatment = metadata$condition
metadata$condition = paste0(metadata$condition, '_', metadata$time)

rm(meta2)

metadata = metadata[, c(1, 10, 14, 11, 2:9, 12:13)]

saveRDS(metadata, file = paste0(RdataDir, '/metadata_cells_CyTOF.rds'))

save(mat, metadata, file = paste0(RdataDir, '/cytof_mat_metadata.Rdata'))

##########################################
# test Robinson's analysis 
# https://www.bioconductor.org/packages/release/workflows/vignettes/cytofWorkflow/inst/doc/cytofWorkflow.html#data-organization
# https://f1000research.com/articles/6-748/v3
##########################################
Use.Robinson.workflow = FALSE
if(Use.Robinson.workflow){
  library(readxl)
  url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
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
  md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))
  md$sample_id <- factor(md$sample_id, 
                         levels = md$sample_id[order(md$condition)])
  
  # construct SingleCellExperiment
  sce0 <- prepData(fs, panel, md, features = panel$fcs_colname)
  
  
  ## load our own data
  #load(file = paste0(RdataDir, '/cytof_mat_metadata.Rdata'))
  #fs <- as(mat, "flowSet") 
  
  load(file = paste0(RdataDir, '/cytof_mat_metadata.Rdata'))
  metadata$time = gsub('_','.', metadata$time)
  metadata$condition = paste0(metadata$treatment, '_', metadata$time)
  #metadata$condition = gsub('beforeRA', 'RA', metadata$condition)
  #metadata$treatment = gsub('beforeRA', 'RA', metadata$treatment)
  
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
  
  sce <- importData(x,
                    assayname = 'exprs',
                    metadata = metadata)
  
  sce$sample_id = as.character(sce@metadata$SampleID)
  cc.levels = c("beforeRA_d2", 
                "noRA_d2.6h", "noRA_d2.12h", "noRA_d2.18h", "noRA_d3",
                "noRA_d3.6h", "noRA_d3.12h", "noRA_d3.18h", "noRA_d4", 
                "RA_d2.6h", "RA_d2.12h", "RA_d2.18h", "RA_d3",
                "RA_d3.6h", "RA_d3.12h", "RA_d3.18h", "RA_d4")
  
  sce$condition = factor(sce@metadata$condition, levels = cc.levels)
  sce$treatment = sce@metadata$treatment
  sce$time = factor(sce@metadata$time, levels = unique(sapply(cc.levels, function(x) 
    {x = unlist(strsplit(as.character(x), '_')); x[length(x)]})))
  
  rowData(sce)$marker_name = rownames(sce)
  rowData(sce)$channel_name = NULL
  rowData(sce)$marker_class = NULL
  # # make singleCellExperiment object
  # subsample = sample(c(1:nrow(x)), size = 15000, replace = FALSE)
  # sce <- importData(x[subsample, ],
  #                   assayname = 'normcounts',
  #                   metadata = metadata[subsample, ])
  # sce
  
  p <- plotExprs(sce, color_by = "condition")
  p$facet$params$ncol <- 2
  p
  
  ggsave(filename = paste0(RdataDir, '/QCs_transformed_intensities_distribution.pdf'), width = 12, height = 8)
 
  n_cells(sce) 
  
  plotCounts(sce, group_by = "sample_id", color_by = "condition")
  ggsave(filename = paste0(RdataDir, '/QCs_cellNumbers_perCondition.pdf'), width = 12, height = 8)
  
  
  # MDS plot
  pbMDS(sce, color_by = "condition", label_by = "sample_id")
  
  plotExprHeatmap(sce, scale = "last",
                  hm_pal = rev(hcl.colors(20, "YlGnBu")))
  
  plotNRS(sce, features = NULL, color_by = "condition")
  
  
  set.seed(1234)
  sce <- CATALYST::cluster(sce, features = NULL,
                 xdim = 10, ydim = 10, maxK = 20, seed = 1234)
  
  set.seed(1234)
  #sce <- runDR(sce, "TSNE", cells = 500, features = "type")
  sce <- runDR(sce, "MDS", cells = 100, features = NULL)
  plotDR(sce, "MDS", color_by = "time", facet_by = 'treatment')
  
  set.seed(1234)
  sce <- runDR(sce, "UMAP", cells = 500, features = NULL, scale = FALSE, n_neighbors = 10, n_threads = 16, min_dist = 0.01)
  plotDR(sce, "UMAP", color_by = "condition", facet_by = 'treatment')
  
  plotDR(sce, "UMAP", color_by = "treatment", facet_by = 'time')
  #plotDR(sce, "UMAP", color_by = "condition")
  
  plotDR(sce, "UMAP", color_by = "FoxA2")
  plotDR(sce, 'UMAP', color_by = 'Pax6')
  
  
  
  
}

##########################################
# normalize data
# original code from https://www.biostars.org/p/393974/
# or https://github.com/kevinblighe/scDataviz
##########################################
Use.scDataviz.analysis = FALSE
if(Use.scDataviz.analysis){
  # After reading in the data, you need to normalise it and eliminate junk cells. 
  # As you know, CyTOF data is typically normalised by hyperbolic arc-sine with a factor of 5.
  
  # exmaple from https://github.com/kevinblighe/scDataviz 
  # mat <- jitter(matrix(
  #   MASS::rnegbin(rexp(50000, rate=.1), theta = 4.5),
  #   ncol = 20))
  # colnames(mat) <- paste0('CD', 1:ncol(mat))
  # rownames(mat) <- paste0('cell', 1:nrow(mat))
  # 
  # metadata <- data.frame(
  #   group = rep('A', nrow(mat)),
  #   row.names = rownames(mat),
  #   stringsAsFactors = FALSE)
  # head(metadata)
  
  
  ##########################################
  # start the normalization and transformation
  ##########################################
  load(file = paste0(RdataDir, '/cytof_mat_metadata.Rdata'))
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
  
  # make singleCellExperiment object
  subsample = sample(c(1:nrow(x)), size = 15000, replace = FALSE)
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
              colby = 'time',
              shape = 'treatment',
              legendPosition = 'right',
              title = 'PCA applied to CyTOF data',
              caption = paste0('10000 cells randomly selected after ',
                               'having filtered for low variance'))
  
  p2 = biplot(p,
              x = 'PC2', y = 'PC3',
              lab = NULL,
              xlim = c(min(p$rotated[,'PC2'])-1, max(p$rotated[,'PC2'])+1),
              ylim = c(min(p$rotated[,'PC3'])-1, max(p$rotated[,'PC3'])+1),
              pointSize = 1.0,
              colby = 'time',
              shape = 'treatment',
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
                    colby = 'treatment',
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

