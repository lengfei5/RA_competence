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

##########################################
# normalize data
# original code from https://www.biostars.org/p/393974/
# or https://github.com/kevinblighe/scDataviz
##########################################
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



