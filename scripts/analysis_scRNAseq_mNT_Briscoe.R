##########################################################################
##########################################################################
# Project: 
# Script purpose: process scRNA-seq data of mouse NT in animal from 
# https://github.com/juliendelile/MouseSpinalCordAtlas#data-availability
# prepare the data as reference for cell type annotation
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep 14 11:39:53 2022
##########################################################################
##########################################################################
# devtools::install_github("juliendelile/Antler", ref = "Development2019")
# 
# library(Antler)


inputDir = '../data/scRNA_NT_mouse/'
source(paste0(inputDir, 'MouseSpinalCordAtlas_tools.R'))

require(data.table)
count.data = fread(input = paste0(inputDir, 'UMI_count.tsv'))
ggs = as.character(count.data$V1)

count.data = as.matrix(count.data[, -1])

## change the gene names
topdir = paste0(dataDir, '/cellRanger_outs/',  design$sampleID[1], '/outs/raw_feature_bc_matrix/')
g = read.csv(paste0(topdir, "/features.tsv.gz"), header = F, stringsAsFactors = F, sep = '\t')
g$name = g$V2
gg.counts = table(g$V2)
gg.dup = names(gg.counts)[which(gg.counts>1)]
index.dup = which(!is.na(match(g$V2, gg.dup)))
g$name[index.dup] = paste0(g$V2[index.dup], '_', g$V1[index.dup])

names = g$name[match(ggs, g$V1)]

jj = which(!is.na(names))

count.data = count.data[jj, ]
rownames(count.data) = names[jj]


meta = read.csv(file = paste0(inputDir, 'phenoData_annotated.csv'), sep = '\t', row.names = 1, header = TRUE)
mm = match(colnames(count.data), rownames(meta))
meta = meta[mm, ]

mnt =  CreateSeuratObject(counts = count.data,
                               meta.data = meta, 
                               min.cells = 20, min.features = 100)

kk = which(mnt$Type_step1 == 'Neuron'|mnt$Type_step1 == 'Progenitor')

mnt = mnt[, kk]


saveRDS(mnt, file = paste0(RdataDir, 'mouse_NT_scRNAseq_NP.Neurons.rds'))

########################################################
########################################################
# Section : utility functions for the function reference.based.cluster.annotation.scmap
#  
########################################################
########################################################
library(SingleCellExperiment)
library(scmap)

cell.sels = which(aa$celltypes == 'NP_RA'|aa$celltypes == 'NP_noRA'|aa$celltypes == 'Neurons')
seurat.obj = subset(x = aa, cells = cell.sels)

# process aleks data for scmap
sce = Seurat::as.SingleCellExperiment(seurat.obj)
sce <- sce[!duplicated(rownames(sce)), ]

rowData(sce)$feature_symbol <- rownames(sce)
counts(sce) = as.matrix(counts(sce)) # sce object converted from seurat object was using spare matrix
logcounts(sce) = as.matrix(logcounts(sce))

ee = readRDS(file = paste0(RdataDir, 'mouse_NT_scRNAseq_NP.Neurons.rds'))
ee = Seurat::as.SingleCellExperiment(ee)
counts(ee) = as.matrix(counts(ee))
logcounts(ee) = as.matrix(logcounts(ee))
rowData(ee)$feature_symbol <- rownames(ee)
ee$cell_type1 = ee$Type_step2

markers = read.csv2(file = '../data/scRNA_NT_mouse/TableS1_binary_knowledge_matrix.csv')
markers = markers$Gene

##########################################
# run scmap-cluster
##########################################
keep = data.frame(colnames(sce), stringsAsFactors = FALSE)
colnames(keep) = 'cell'

for(nb.features.scmap in c(20, 30, 40, 50, 100, 200, 500))
{
  ## feature selection for scmap
  # nb.features.scmap = 3000
  #threshold.scmap = 0.5
  cat('nb of features selected : ', nb.features.scmap, '\n')
  
  #ee_ref <- selectFeatures(ee, suppress_plot = FALSE, n_features = nb.features.scmap)
  ee_ref = setFeatures(ee, features = intersect(markers, rownames(ee)))
  #table(rowData(ee)$scmap_features)
  ee_ref = indexCluster(ee_ref, cluster_col = 'Type_step2')
  
  #head(metadata(ee_ref)$scmap_cluster_index)
  #heatmap(as.matrix(metadata(ee_ref)$scmap_cluster_index))
  
  scmapCluster_results <- scmapCluster(
    projection = sce, 
    index_list = list(
      ref = metadata(ee_ref)$scmap_cluster_index
    ),
    threshold = 0
  )
  
  #seurat.obj = AddMetaData(seurat.obj, as.factor(scmapCluster_results$scmap_cluster_labs), 
  #                         col.name = paste0('scmap.pred.id.features.', nb.features.scmap))
  
  keep = data.frame(keep, as.character(scmapCluster_results$scmap_cluster_labs), 
                    as.numeric(scmapCluster_results$scmap_cluster_siml), 
                    stringsAsFactors = FALSE)
  
  #length(scmapCluster_results$scmap_cluster_labs)
  #length(scmapCluster_results$combined_labs)
  ident.ref = unique(ee$cell_type1)
  ident.projection = unique(scmapCluster_results$scmap_cluster_labs)
  ident.missed = ident.ref[which(is.na(match(ident.ref, ident.projection)))]
  cat('cell identities missed : \n')
  print(ident.missed)
  #head(scmapCluster_results$scmap_cluster_labs)
  #head(scmapCluster_results$scmap_cluster_siml)
  
  hist(scmapCluster_results$scmap_cluster_siml, breaks = 100)
  #abline(v = threshold.scmap, col = 'red')
  head(scmapCluster_results$combined_labs)
  
  predicted.id = scmapCluster_results$scmap_cluster_labs
  counts.pred.ids = table(predicted.id)
  counts.pred.ids = counts.pred.ids[order(-counts.pred.ids)]
  # print(counts.pred.ids)
  
  predicted.id[which(predicted.id == 'unassigned')] = NA
  
  cat('nb of assigned cells :',  length(predicted.id[!is.na(predicted.id)]), '\n')
  cat('percent of assigned cells: ', length(predicted.id[!is.na(predicted.id)])/length(predicted.id), '\n')
  
}

colnames(keep)[-1] = paste0(rep(c('scmap.pred.id.', 'scmap.corr.'), 7), 
                            rep(c(20,30,40, 50, 100, 200, 500), each =2))

keep = keep[, c(1, 4:5)]
colnames(keep)[-1] = c('scmap.pred.id', 'scmap.corr')
rownames(keep) = colnames(seurat.obj)
keep = keep[, c(1, grep('.30|.50|.100|.500', colnames(keep)))]
keep = as.data.frame(keep)
rownames(keep) = keep$cell # metadata rownames must be the cell barcodes
seurat.obj= AddMetaData(seurat.obj, metadata = keep[match(colnames(seurat.obj), keep$cell), -1])


sub.obj = seurat.obj
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

DimPlot(sub.obj, group.by = 'scmap.pred.id')

saveRDS(seurat.obj, file = paste0(RdataDir, 
                                  'annotation_ref.scmap.rds'))

##########################################
# run scmap-cell 
##########################################
run.scamp.cell = FALSE
if(run.scmap.cell){
  set.seed(1)
  
  ## feature selection for scmap
  ee <- selectFeatures(ee, suppress_plot = FALSE, n_features = 3000)
  table(rowData(ee)$scmap_features)
  #as.character(unique(ee$cell_type1))
  
  ee <- indexCell(ee)
  
  names(metadata(ee)$scmap_cell_index)
  
  length(metadata(ee)$scmap_cell_index$subcentroids)
  
  dim(metadata(ee)$scmap_cell_index$subcentroids[[1]])
  
  metadata(ee)$scmap_cell_index$subcentroids[[1]][,1:5]
  
  scmapCell_results <- scmapCell(
    sce, 
    list(
      ref = metadata(ee)$scmap_cell_index
    ),
    w = 10
  )
  
  scmapCell_clusters <- scmapCell2Cluster(
    scmapCell_results, 
    list(
      as.character(colData(ee)$cell_type1)
    ),
    w = 3,
    threshold = 0.
  )
  
  head(scmapCell_clusters$scmap_cluster_labs)
  table(scmapCell_clusters$scmap_cluster_labs)
  
  head(scmapCell_clusters$scmap_cluster_siml)
  
}

