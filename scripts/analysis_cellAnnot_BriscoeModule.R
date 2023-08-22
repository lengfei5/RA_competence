##########################################################################
##########################################################################
# Project: RA competence  
# Script purpose: annotate the terminal cells with James's group knowlegde-based module
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Sep 30 16:56:41 2022
##########################################################################
##########################################################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene_regressed.nCounts_clusterd0.7_manualAnnot_v1_', 
                           species, version.analysis, '.rds'))


source('../data/MouseSpinalCordAtlas/R_files/MouseSpinalCordAtlas_tools.R')
partition_input = '../data/MouseSpinalCordAtlas/input_files/partitioning_table.csv'
refs = read.csv(partition_input, sep = '\t')
refs = refs[which(refs$Type == 'Progenitor'|refs$Type == 'Neuron'), ]


cell_partition = doCellPartition(known_template_file="./input_files/partitioning_table.csv", 
                                 readcounts=m$getReadcounts(data_status='Raw'))

pop_colors = getPopulationColors(known_template_file="./input_files/partitioning_table.csv")

for(md in c("Type_step1", "Type_step2", "Type_step2_unique", "DV")){pData(m$expressionSet)[[md]] <- cell_partition[[md]]}

##########################################
# subset the neurons and NP cells in cluster 11, 12 and 13
##########################################
cluster_sels = c('11', '12', '13')
Idents(aa) = factor(aa$clusters)
sub.obj = subset(x = aa, idents = cluster_sels)

# sub.obj = NormalizeData(sub.obj, normalization.method = "LogNormalize", scale.factor = 10000)
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
sub.obj = ScaleData(sub.obj, vars.to.regress = 'nCount_RNA')

sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 30)

nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
n.neighbors = 30; min.dist = 0.1;

sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                   dims = 1:nb.pcs, 
                   n.neighbors = n.neighbors,
                   min.dist = min.dist)

sub.obj <- FindNeighbors(sub.obj, dims = 1:nb.pcs)
sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 1.0)
