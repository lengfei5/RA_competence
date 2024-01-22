##########################################################################
##########################################################################
# Project: RA competence project
# Script purpose: mapping mNT cells into mouse gastrulation datsets 
# check the Pax6 and Fox2 co-expression in mouse data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Mar 29 11:02:25 2023
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_MouseGastrulationData/'

#resDir = paste0("../results/dataset_scRNAseq", version.analysis)
resDir = paste0('../results/scRNAseq_R13547_10x_mNT_20220813/mapping_to', version.analysis)
RdataDir = paste0(resDir, 'Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source(paste0(functionDir, '/functions_scRNAseq.R'))
source(paste0(functionDir, '/functions_Visium.R'))
source(paste0(functionDir, '/functions_dataIntegration.R'))

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(DropletUtils)
library(edgeR)
library(future)
library(data.table)
library(tidyverse)

options(future.globals.maxSize = 160000 * 1024^2)
mem_used()

species = '_mouseGastrulation'

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

## subset our scRNA-seq data 
levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA",
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA")

#levels_sels = c("day2_beforeRA",  "day2.5_RA", "day3_RA.rep1", "day3.5_RA",
#                "day4_RA", "day5_RA", "day6_RA")

cols_sel = cols[match(levels_sels, names(cols))]


col_mouse = readRDS(file = paste0(RdataDir, 'cols_mouse_gastrulation_celltypes.rds'))

sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))


##########################################
# import mNT scRNA-seq data
##########################################
aa =  readRDS(file = paste0('../results/Rdata/',  
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', 'mNT_scRNAseq',
                            '_R13547_10x_mNT_20220813', '.rds'))


Idents(aa) = factor(aa$condition)
aa = subset(aa, idents = levels_sels)

# downsample for each condition
aa = subset(x = aa, downsample = 1000)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition.rds'))

########################################################
########################################################
# Section I: Explore the mapping our mNT scRNA-seq data to mouse gastrulation atlas
# Reference choices
# 1) full atlas from Marioni2019 or Chan2019 or integrated Marioni.Chan
# 2) selected relevant cell types from those full atlas 
# Methods choices: seurat_CCA, seurat_RPCA, harmony, fastMNN, scanorama, scVI
########################################################
########################################################

##########################################
# import the mNT data and reference
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition_reclustered.rds'))

ref = readRDS(file = paste0(RdataDir,  
                             'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v2.rds'))

cols_mouse = sapply(ref$colour, function(x) {paste0('#', x, collapse = '')})
names(cols_mouse) =ref$celltype
cols_mouse = cols_mouse[match(unique(names(cols_mouse)), names(cols_mouse))]


mapping_method = "seurat_rpca"
data_version = 'mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes'

outDir = paste0(resDir, mapping_method, '/',  data_version, '/')
system(paste0('mkdir -p ', outDir))

##########################################
# test Seurat data integration
##########################################
features.common = intersect(rownames(aa), rownames(ref))

aa = subset(aa, features = features.common)
ref = subset(ref, features = features.common)

aa$dataset = 'mNT'
aa$stage = aa$condition
aa$sequencing.batch = 'mNT'
ref$dataset = 'ref'

aa$celltype = paste0('mNT_', aa$condition)

##########################################
# calculate similarity before data integration
##########################################
Test_cal_similarity = FALSE
if(Test_cal_similarity){
  # aa = FindVariableFeatures(aa, selection.method = "vst")
  # aa <- ScaleData(aa, verbose = FALSE)
  # aa <- RunPCA(aa, npcs = 50, verbose = FALSE)
  # ElbowPlot(aa, ndims = 50)
  # 
  # aa <- RunUMAP(aa, reduction = "pca", dims = 1:30, n.neighbors = 50, 
  #               min.dist = 0.2) 
  # DimPlot(aa, group.by = 'condition')
  # 
  # aa <- FindNeighbors(aa, reduction = "pca", dims = 1:20)
  # aa <- FindClusters(aa, resolution = 0.7)
  # aa$clusters = aa$seurat_clusters
  #saveRDS(aa, file = paste0(RdataDir, 
  #                          'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition_reclustered.rds'))
  
  p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
  p2 = DimPlot(aa, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)
  p1 + p2
  
  ggsave(paste0(outDir, '/umap_query_conditions_clusters.pdf'), 
         width = 16, height = 6)
  
  FeaturePlot(aa, features = 'Foxa2')
  
  source(paste0(functionDir, '/functions_dataIntegration.R'))
  
  ref = FindVariableFeatures(ref, selection.method = 'vst', nfeatures = 1000)
  
  cc = c(3, 6, 5, 1, 7, 4, 11, 10)
  
  for(n in 1:length(cc))
  {
    # n = 4
    cat(n, ' -- ', cc[n], '\n')
    
    subs = subset(aa, cells = colnames(aa)[which(aa$clusters == cc[n])]);
    px = calculate_similarity_query_ref(query = subs, 
                                        ref = ref, 
                                        nHVGs = 1000, 
                                        method = c("spearman"),
                                        group.by = 'celltype')
    
    pdfname = paste0(outDir, '/spearman_similarity_withRefCelltypes_', cc[n], '.pdf')
    
    pdf(pdfname, width=16, height = 8)
    plot(px)
    
    dev.off()
    
  }
    
}

##########################################
# test integration method
##########################################
refs.merged = merge(aa, y = ref, add.cell.ids = c("mNT", "mouseGastrulation"), project = "RA_competence")

ref.list <- SplitObject(refs.merged, split.by = "dataset")

rm(list = c('refs.merged')) # remove big seurat objects to clear memory

# normalize and identify variable features for each dataset independently
ref.list <- lapply(X = ref.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize")
  x <- FindVariableFeatures(x, selection.method = "vst")
  
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ref.list)

ref.list <- lapply(X = ref.list, FUN = function(x) {
  x <- ScaleData(x, features = features.common, verbose = TRUE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  
})

ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                      anchor.features = features, 
                                      #reference = c(2),
                                      #reduction = "cca", 
                                      reduction = 'rpca',
                                      k.anchor = 5,
                                      dims = 1:50)

rm(ref.list)

# this command creates an 'integrated' data assay
ref.combined <- IntegrateData(anchorset = ref.anchors, features.to.integrate = features.common, 
                              dims = 1:50) ## take ~100G memory

rm(ref.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(ref.combined) <- "integrated"

ref.combined <- ScaleData(ref.combined, verbose = FALSE)
ref.combined <- RunPCA(ref.combined, npcs = 50, verbose = FALSE)

ElbowPlot(ref.combined, ndims = 50)

kk = which(ref.combined$dataset == 'mNT') 
ref.combined$celltype[kk] = paste0('mNT_', ref.combined$condition[kk])
names(cols_sel) = paste0('mNT_', names(cols_sel))

cols_used = c(cols_mouse, cols_sel)
#ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
#ref.combined <- FindClusters(ref.combined, resolution = 0.2)

ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:50, n.neighbors = 50, 
                        min.dist = 0.2) 

saveRDS(cols_used, file = paste0(outDir, '/integrated_mNT_mouseGastrulation_colorsUsed.rds'))

# Visualization
DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
        repel = TRUE, raster=FALSE, cols = c(cols_mouse, cols_sel)) 

ggsave(paste0(outDir, '/Integration_mNT_celltypes_noref.batchcorrection.pdf'), 
       width = 16, height = 8)

#DimPlot(ref.combined, reduction = "umap")
saveRDS(ref.combined, file = paste0(outDir, '/integrated_mNT_mouseGastrulation_SeuratRPCA.rds'))


ref.combined = readRDS(file = paste0(outDir, '/integrated_mNT_mouseGastrulation_SeuratRPCA.rds'))
cols_used = readRDS(file = paste0(outDir, '/integrated_mNT_mouseGastrulation_colorsUsed.rds'))


# Visualization
DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
        repel = TRUE, raster=FALSE, cols = cols_used) 

ggsave(paste0(outDir, '/Integration_mNT_celltypes_noref.batchcorrection.pdf'), 
       width = 16, height = 8)

DimPlot(ref.combined, reduction = "umap", group.by = "dataset", raster=FALSE)
ggsave(paste0(outDir, '/Integration_dataset.pdf'), 
       width = 16, height = 8)

ggs = c('Pax6', 'Foxa2', 'Pou5f1', 'Sox17', 'Sox1', 'Sox2', 
        "Ifitm1", "T",
"Krt8", 'Tuba1a', "Eno1", 'Krt18','Foxj1', 'Sp5', 'Noto', 'Pou5f1','Foxa2','Cdx2', 'Gsta4', 'Sox9','Fst',
'Nog', 'Shh', 'Slc2a1', 'Cxx1b','Igfbp2','Epcam', 'Lhx1','Ptch1',
'Flrt3','Foxp1', 'Hoxb1', 'Dnmt3b', 'Nr6a1','Cdh2','Rspo3','Zfp503','Fgf5')

ggs = ggs[which(!is.na(match(ggs, rownames(ref.combined))))]

pdf(paste0(outDir, '/FeaturePlot_Markers.pdf'),
    width =10, height = 8, useDingbats = FALSE)
for(n in 1:length(ggs))
{
  cat(n, '--', ggs[n], '\n')
  p1 = FeaturePlot(ref.combined, features = ggs[n], min.cutoff = 'q5')
  #FeaturePlot(ref.combined, features = 'Foxa2', min.cutoff = 'q5')
  #FeaturePlot(ref.combined, features = 'Sox17', min.cutoff = 'q5')
  plot(p1)
  
}

dev.off()

DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE, split.by = 'dataset', 
        cols = cols_used,
        repel = TRUE, raster=FALSE) + NoLegend()

ggsave(paste0(outDir, '/Integration_celltypes_split.dataset.pdf'), 
       width = 24, height = 8)

DimPlot(ref.combined, reduction = "umap", group.by = "stage", label = TRUE,
        repel = TRUE, raster=FALSE)

ggsave(paste0(outDir, '/Integration_stage.pdf'), 
       width = 16, height = 8)


##########################################
# test seurat query-ref-mapping 
# test projection method (slight different from integration)
# https://satijalab.org/seurat/articles/integration_mapping.html (original code)
# did not work well and don't know the reason
##########################################
Test_Seurat_query_ref_mapping = FALSE 
if(Test_Seurat_query_ref_mapping){
  mapping_method = "seurat_query_ref_mapping"
  
  aa = readRDS(file = paste0(RdataDir, 
                             'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition_reclustered.rds'))
  ref = readRDS(file = paste0(RdataDir,  
                              'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v3.rds'))
  
  data_version = 'mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes_test2'
  
  features.common = intersect(rownames(aa), rownames(ref))
  
  aa = subset(aa, features = features.common)
  ref = subset(ref, features = features.common)
  
  aa$dataset = 'mNT'
  aa$stage = aa$condition
  aa$sequencing.batch = 'mNT'
  ref$dataset = 'ref'
  
  aa$celltype = paste0('mNT_', aa$condition)
  
  outDir = paste0(resDir,  mapping_method, '/', data_version, '/')
  system(paste0('mkdir -p ', outDir))
  
  
  ElbowPlot(ref, ndims = 50, reduction = 'pca')
  ref = RunUMAP(ref, reduction = "pca", dims = 1:30, n.neighbors = 30, 
                min.dist = 0.1, return.model = TRUE) 
  
  DimPlot(ref, reduction = "umap", 
          group.by = "celltype", label = TRUE,
          repel = TRUE, raster=FALSE, cols = cols_mouse) 
  
  ref$labels = factor(paste0(ref$celltype, '_', ref$stage))
  
  # In data transfer, Seurat has an option (set by default) to project the PCA structure of a reference 
  # onto the query, instead of learning a joint structure with CCA. 
  # We generally suggest using this option when projecting data between scRNA-seq datasets.
  anchors <- FindTransferAnchors(reference = ref, 
                                 query = aa, 
                                 dims = 1:50,
                                 normalization.method = "LogNormalize",
                                 #reference.reduction = "pca.corrected",
                                 max.features = 200,
                                 k.anchor = 10,
                                 reduction = "cca" 
  )
  
  #pancreas.anchors <- FindTransferAnchors(reference = pancreas.ref, query = pancreas.query, dims = 1:30,
  #                                        reference.reduction = "pca")
  predictions <- TransferData(anchorset = anchors, refdata = ref$labels, 
                              reference = ref,
                              query = aa,
                              weight.reduction = 'rpca.ref',
                              dims = 1:30)
  
  saveRDS(anchors, file = paste0(outDir, '/cca_anchors.rds'))
  query = aa;
  query$predicted.id = predictions$predicted.id 
  query$predicted.score = predictions$predicted.id.score
  
  
  p1 = DimPlot(query, reduction = "umap", 
          group.by = "predicted.id", label = TRUE,
          repel = TRUE, raster=FALSE) 
  p2 = FeaturePlot(query, features = 'predicted.score')
  
  p1
  ggsave(paste0(outDir, '/transfer_learning_seurat_cca.ref_v3.pdf'), 
         width = 16, height = 8)
  
  
  
  ## visualize the query cells alongside our reference and didn't work well
  ## (MapQuery in Seurat didn't work well)
  # query <- Seurat::MapQuery(anchorset = anchors, 
  #                   reference = ref, 
  #                   query = aa,
  #                   refdata = list(labels = "labels"),
  #                   #refdata = list(celltype = "celltype"), 
  #                   transferdata.args = list(weight.reduction = 'rpca.ref'), 
  #                   reference.reduction = "pca", 
  #                   reduction.model = "umap")
  # 
  # p1 <- DimPlot(ref, reduction = "umap", group.by = "celltype", 
  #               label = TRUE, label.size = 3,
  #               repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
  # p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.labels", 
  #               label = TRUE,
  #               label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
  # p1 + p2
  
  
}

##########################################
# test the symphony 
##########################################
Test_reference_mapping_Symphony = FALSE
if(Test_reference_mapping_Symphony){
  library(symphony)
  library(singlecellmethods)
  source('utils_symphony.R')
  
  fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
  }
  
  mapping_method = 'symphony_mapping'
  data_version = "mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes"
  
  outDir = paste0(resDir,  mapping_method, '/', data_version, '/')
  system(paste0('mkdir -p ', outDir))
  
  ## import and prepare the ref and query
  aa = readRDS(file = paste0(RdataDir, 
                             'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition_reclustered.rds'))
  ref = readRDS(file = paste0(RdataDir,  
                              'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v3.rds'))
  
  features.common = intersect(rownames(aa), rownames(ref))
  aa = subset(aa, features = features.common)
  ref = subset(ref, features = features.common)
  
  aa$dataset = 'mNT'
  aa$stage = aa$condition
  aa$sequencing.batch = 'mNT'
  ref$dataset = 'ref'
  aa$celltype = paste0('mNT_', aa$condition)
  
  ref$labels = paste0(ref$celltype, '_', ref$stage)
  
  p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
  p2 = DimPlot(aa, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)
  p1 + p2
  
  ggsave(paste0(outDir, '/umap_query_conditions_clusters.pdf'), 
         width = 16, height = 6)
  
  
  
  # Read in normalized expression and metadata
  #exprs_norm = readRDS('../data/data_symphony/exprs_norm_all.rds')
  #metadata = read.csv('../data/data_symphony/meta_data_subtypes.csv', row.names = 1)
  #dim(exprs_norm)
  #dim(metadata)
  
  #idx_query = which(metadata$donor == "5'") # use 5' dataset as the query
  ref_exp_full = ref@assays$RNA@data
  ref_metadata = ref@meta.data
  query_exp = aa@assays$RNA@data
  query_metadata = aa@meta.data
  
  # Sparse matrix with the normalized genes x cells matrix
  ref_exp_full[1:5, 1:2]
  
  # Select variable genes and subset reference expression by variable genes
  var_genes = vargenes_vst(ref_exp_full, groups = as.character(ref_metadata[['sample']]), topn = 500)
  ref_exp = ref_exp_full[var_genes, ]
  dim(ref_exp)
  
  # Build reference
  reference = symphony::buildReference(
    ref_exp,                   # reference expression (genes by cells)
    ref_metadata,              # reference metadata (cells x attributes)
    vars = NULL,         # variable(s) to integrate over
    K = 100,                   # number of Harmony soft clusters
    verbose = TRUE,            # display verbose output
    do_umap = TRUE,            # run UMAP and save UMAP model to file
    do_normalize = FALSE,      # perform log(CP10k) normalization on reference expression
    vargenes_method = 'vst',   # variable gene selection method: 'vst' or 'mvp'
    vargenes_groups = 'sample', # metadata column specifying groups for variable gene selection within each group
    topn = 500,               # number of variable genes (per group)
    theta = 2,                 # Harmony parameter(s) for diversity term
    d = 20,                    # number of dimensions for PCA
    save_uwot_path = '/groups/tanaka/People/current/jiwang/projects/RA_competence/results/scRNAseq_R13547_10x_mNT_20220813/mapping_to_MouseGastrulationData/symphony_mapping/mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes/model_test'
    , # file path to save uwot UMAP model
    additional_genes = NULL    # vector of any additional genes to force include
  )
  
  #Run Harmony integration
  Run_Harmony_integration = FALSE
  if(Run_Harmony_integration){
    # Calculate and save the mean and standard deviations for each gene
    vargenes_means_sds = tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
    vargenes_means_sds$stddev = singlecellmethods::rowSDs(ref_exp,  row_means = vargenes_means_sds$mean)
    head(vargenes_means_sds)
    
    #Scale data using calculated gene means and standard deviations
    ref_exp_scaled = singlecellmethods::scaleDataWithStats(ref_exp, vargenes_means_sds$mean, 
                                                           vargenes_means_sds$stddev, 1)
    
    #Run SVD, save gene loadings (s$u)
    set.seed(0)
    s = irlba::irlba(ref_exp_scaled, nv = 20)
    Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
    loadings = s$u
    
    set.seed(0)
    ref_harmObj = harmony::HarmonyMatrix(
      data_mat = t(Z_pca_ref),  ## PCA embedding matrix of cells
      meta_data = ref_metadata, ## dataframe with cell labels
      theta = c(2),             ## cluster diversity enforcement
      vars_use = NULL,    ## variable to integrate out
      nclust = NULL,             ## number of clusters in Harmony model
      max.iter.harmony = 20,
      return_object = TRUE,     ## return the full Harmony model object
      do_pca = FALSE            ## don't recompute PCs
    )
    
    # To run the next function buildReferenceFromHarmonyObj(), 
    # you need to input the saved gene loadings (loadings) and vargenes_means_sds.
    # Compress a Harmony object into a Symphony reference
    reference = symphony::buildReferenceFromHarmonyObj(
      ref_harmObj,            # output object from HarmonyMatrix()
      ref_metadata,           # reference cell metadata
      vargenes_means_sds,     # gene names, means, and std devs for scaling
      loadings,               # genes x PCs matrix
      verbose = TRUE,         # verbose output
      do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
      save_uwot_path = '/groups/tanaka/People/current/jiwang/projects/RA_competence/results/scRNAseq_R13547_10x_mNT_20220813/mapping_to_MouseGastrulationData/symphony_mapping/mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes/model_test')
    
  }
  
  # Optionally, you can specify which normalization method was
  # used to build the reference as a custom slot inside the Symphony object to 
  # help record this information for future query users
  #reference$normalization_method = 'log(CP10k+1)'
  saveRDS(reference, paste0(outDir, 'testing_reference_mouseGastrulation_1.rds'))
  
  str(reference)
  
  # The harmonized embedding is located in the Z_corr slot of the reference object.
  dim(reference$Z_corr)
  reference$Z_corr[1:5, 1:5]
  
  # reference = readRDS(paste0(outDir, 'testing_reference1.rds'))
  umap_labels = cbind(ref_metadata, reference$umap$embedding)
  
  fig.size(3, 5)
  plotBasic(umap_labels, title = 'Reference', color.by = 'celltype')
  
  # In order to map a new query dataset onto the reference, 
  # you will need a reference object saved from the steps above, 
  # as well as query cell expression and metadata.
  # Map query
  query = mapQuery(query_exp,             # query gene expression (genes x cells)
                   query_metadata,        # query metadata (cells x attributes)
                   reference,             # Symphony reference object
                   do_normalize = FALSE,  # perform log(CP10k+1) normalization on query
                   do_umap = TRUE)        # project query cells into reference UMAP
  
  ## Symphony assumes that the query is normalized in the same manner as the reference. 
  # Our implementation currently uses log(CP10k) normalization.
  
  str(query)
  
  #Let's take a look at what the query object contains:

   # Z: query cells in reference Harmonized embedding
   #Zq_pca: query cells in pre-Harmony reference PC embedding (prior to correction)
   # R: query cell soft cluster assignments
   #  Xq: query cell design matrix for correction step
   #  umap: query cells projected into reference UMAP coordinates (using uwot)
   #  meta_data: metadata
  # Predict query cell types using k-NN
  query = knnPredict(query, reference, reference$meta_data$celltype, k = 5)
  
  ## Query cell type predictions are now in the cell_type_pred_knn column. 
  ## The cell_type_pred_knn_prob column reports the proportion of nearest neighbors with the winning vote 
  # (can help identify query cells that fall "on the border" between 2 reference cell types).
  
  head(query$meta_data)
  
  # Add the UMAP coordinates to the metadata
  reference$meta_data$cell_type_pred_knn = NA
  reference$meta_data$cell_type_pred_knn_prob = NA
  reference$meta_data$ref_query = 'reference'
  query$meta_data$ref_query = 'query'
  
  # Add the UMAP coordinates to the metadata
  umap_combined = rbind(query$umap, reference$umap$embedding)
  xx = query$meta_data[, c(24, 22)]
  colnames(xx)[2] = 'celltype' 
  meta_data_combined = rbind(xx, reference$meta_data[, c(27, 18)])
  
  umap_combined_labels = cbind(meta_data_combined, umap_combined)
  
  fig.size(6, 14)
  plotBasic(umap_combined_labels, title = 'Reference and query cells', 
            color.by = 'celltype', facet.by = 'ref_query')
  
  ggsave(paste0(outDir, '/mapping_symphony_ref_v1_celltype_coembedding.pdf'), 
         width = 16, height = 8)
  
  aa$predicted.id = query$meta_data$cell_type_pred_knn
  aa$predicted.score = query$meta_data$cell_type_pred_knn_prob
      
  
  p2 = FeaturePlot(aa, features = 'predicted.score')
  
  p1
  p1 = DimPlot(aa, reduction = "umap", 
          group.by = "predicted.id", label = TRUE,
          repel = TRUE, raster=FALSE) 
  plot(p1)
  ggsave(paste0(outDir, '/mapping_symphony_ref_celltype_v1.pdf'), 
         width = 16, height = 8)
  
  
}

##########################################
# test scArches and bipartite graph
##########################################



##########################################
# test the similarity calculation between mNT cells and cell types in the reference
# after data integration
# 
##########################################
DefaultAssay(ref.combined) = 'integrated'
#ref.combined = FindVariableFeatures(ref.combined, selection.method = 'vst', nfeatures = 1000)

refs_subs = subset(ref.combined, cells = colnames(ref.combined)[grep('mouseGastrulation', 
                                                                    colnames(ref.combined))]);
cc = c(3, 6, 5, 1, 7, 4, 11, 10)
source(paste0(functionDir, '/functions_dataIntegration.R'))

refs_subs = FindVariableFeatures(refs_subs, selection.method = 'vst', nfeatures = 3000, assay = 'integrated')
ggs = VariableFeatures(refs_subs)
ggs = intersect(ggs, c(sps, tfs))

for(n in 1:length(cc))
{
  # n = 1
  cat('cluster -- ', cc[n], '\n')
  
  cells = colnames(aa)[which(aa$clusters == cc[n])]
  mm1 = match(colnames(ref.combined), paste0('mNT_', cells))
  subs = subset(ref.combined, cells = colnames(ref.combined)[!is.na(mm1)]);
  #subs = FindVariableFeatures(subs, selection.method = 'vst', nfeatures = 2000)
  
  px = calculate_similarity_query_ref(query = subs, 
                                      ref = refs_subs, 
                                      assay_use = 'integrated',
                                      find_hvg = FALSE,
                                      features_use = ggs,
                                      method = c("pearson"),
                                      group.by = 'celltype')
  
  pdfname = paste0(outDir, '/spearman_similarity_withRefCelltypes_dataIntegration_', cc[n], '_test_4.pdf')
  
  pdf(pdfname, width=16, height = 8)
  plot(px)
  
  dev.off()
  
}

##########################################
# clustering the combined data
##########################################
ref.combined = readRDS(file = paste0(outDir, '/integrated_mNT_mouseGastrulation_SeuratRPCA.rds'))

ElbowPlot(ref.combined, ndims = 50)
ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
ref.combined <- FindClusters(ref.combined, resolution = 1.0)

DimPlot(ref.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
        repel = TRUE, raster=FALSE)

ggsave(paste0(outDir, '/integrated_ref_mNT_27clusters_resolution1.0.pdf'), 
       width = 16, height = 10)

cluster19.markers <- FindMarkers(ref.combined, ident.1 = 19)
head(cluster19.markers, n = 10)

cluster.markers <- FindMarkers(ref.combined, ident.1 = 19, ident.2 = c(16, 24))
head(cluster.markers, n = 10)

markers <- FindAllMarkers(ref.combined, only.pos = TRUE)

save(cluster.markers, cluster19.markers, markers, 
     file = paste0(outDir, '/marker_Genes.Rdata'))

markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 30) %>%
  ungroup() -> top10

DoHeatmap(ref.combined, features = top10$gene) + NoLegend()
ggsave(paste0(outDir, '/markerGenes_all.clusters.for.cluster19_top30.pdf'), 
       width = 16, height = 30)


########################################################
########################################################
# Section II: first check Pax6-FoxA2 positive cells
# 
########################################################
########################################################
srat = readRDS(file = paste0(RdataDir,  
                             'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v2.rds'))

##########################################
# check the Pax6 and FoxA2 co-expression
##########################################
p1 = VlnPlot(srat, features = c('Pax6'), group.by = 'celltype') + NoLegend()
p2 = VlnPlot(srat, features = 'Foxa2', group.by = 'celltype') + NoLegend()

p1 / p2

ggsave(filename = paste0(resDir, '/Vlnplot_Pax6_Foxa2_celltypes.pdf'), width = 14, height = 10)


FeaturePlot(srat, features = c('Pax6', 'Foxa2'), blend = TRUE, raster = TRUE, order = TRUE)

ggsave(filename = paste0(resDir, '/Featureplots_Pax6_Foxa2_blended_ordered.pdf'), width = 20, height = 8)


## subset potential clusters
Idents(srat) = as.factor(srat$celltype)

celltypes_sels = c('Surface ectoderm', 'Spinal cord', 'Rostral neurectoderm',
                   'NMP', 'Forebrain/Midbrain/Hindbrain', 
                   'Caudal neurectoderm', 'Caudal epiblast')

sub.obj = subset(srat, idents = celltypes_sels)

Idents(sub.obj) = as.factor(sub.obj$celltype)

DimPlot(sub.obj, reduction = 'umap', group.by = 'celltype', 
        label = TRUE, repel = TRUE,
        raster=FALSE)


sub.obj <- RunUMAP(sub.obj, reduction = "pca.corrected", dims = 1:50, 
                   n.neighbors = 30, min.dist = 0.1)

p1 = DimPlot(sub.obj, reduction = 'umap', group.by = 'celltype', 
        label = TRUE, repel = TRUE,
        raster=FALSE)
p2 = DimPlot(sub.obj, reduction = 'umap', group.by = 'stage', 
             label = TRUE, repel = TRUE,
             raster=FALSE)

p1 + p2

ggsave(filename = paste0(resDir, '/umap_subsetting_celltypes_stages.pdf'), width = 25, height = 8)

FeaturePlot(sub.obj, features = c('Pax6', 'Foxa2'), blend = TRUE, raster = TRUE, order = TRUE)
ggsave(filename = paste0(resDir, '/Featureplots_Pax6_Foxa2_blended_ordered_subsetting.pdf'), 
       width = 20, height = 8)


