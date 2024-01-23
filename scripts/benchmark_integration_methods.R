##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: to compare and benchmark the integration and projection (mapping methods)
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Nov 23 11:59:22 2023
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
library('pals')

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
names(cols_sel) = paste0('mNT_', names(cols_sel))

# xx = readRDS(paste0(RdataDir, 'seuratObject_EmbryoAtlasData_all36sample_RNAassay.rds'))
# cols_mouse = sapply(xx$colour, function(x) {paste0('#', x, collapse = '')})
# cols_mouse = unique(c(cols_mouse, cols25(n = 25)))
# saveRDS(cols_mouse, file = paste0(RdataDir, 'cols_mouse_gastrulation_celltypes.rds'))

col_mouse = readRDS(file = paste0(RdataDir, 'cols_mouse_gastrulation_celltypes.rds'))

sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

########################################################
########################################################
# Section : test seurat_RPCA
# 
########################################################
########################################################
Test_seurat_RPCA_Chan2019 = FALSE
if(Test_seurat_RPCA_Chan2019){
  ##########################################
  # import the mNT data and reference
  ##########################################
  mapping_method = "seurat_rpca"
  aa = readRDS(file = paste0(RdataDir, 'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition.rds'))
  
  ## Marioni2019 subset data
  ref = readRDS(file = paste0(RdataDir,  
                              'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v2.rds'))
  
  ## Chan2019 full data
  ref = readRDS(file = paste0('../results/Rdata/',  
                              'seuratObject_mm10_mouse_gastrulation_Chan.et.al_',
                              'lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltypes_fastmnn.rds'))
  
  data_version = 'mapping_mNT.noRA.RA.d2_d5_Chan2019_allCelltypes_TFs.SPs.only'
  
  outDir = paste0(resDir, mapping_method, '/',  data_version, '/')
  system(paste0('mkdir -p ', outDir))
  
  ##########################################
  # test Seurat data integration
  ##########################################
  features.common = intersect(rownames(aa), rownames(ref))
  features.common = intersect(features.common, c(tfs, sps))
  
  aa = subset(aa, features = features.common)
  ref = subset(ref, features = features.common)
  
  aa$dataset = 'mNT'
  aa$stage = aa$condition
  aa$sequencing.batch = 'mNT'
  ref$dataset = 'ref'
  
  aa$celltype = paste0('mNT_', aa$condition)
  
  celltypes = unique(ref$celltype)
  cols_mouse = cols_mouse[1:length(celltypes)]
  names(cols_mouse) = celltypes
  rm(celltypes)
  
  ##########################################
  # test integration method
  ##########################################
  refs.merged = merge(aa, y = ref, add.cell.ids = c("mNT", "mouseGastrulation"), project = "RA_competence")
  
  source(paste0(functionDir, '/functions_dataIntegration.R'))
  
  ref.combined = IntegrateData_Seurat_RPCA(seuratObj = refs.merged, group.by = 'dataset', nfeatures = 1000)
  rm(refs.merged)
  
  jj = which(ref.combined$dataset == 'mNT')
  ref.combined$celltype[jj] = paste0('mNT_', ref.combined$condition[jj])
  
  # Visualization
  DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
          repel = TRUE, raster=FALSE, cols = c(cols_sel, cols_mouse)) + NoLegend()
  
  ggsave(paste0(outDir, '/Integration_mNT_celltypes_noref.batchcorrection_1000HVGs.pdf'), 
         width = 16, height = 8)
  
  
  ## highlight the RA cells in the umap
  Idents(ref.combined) = factor(ref.combined$condition)
  for(cc in c('day3_RA.rep1', 'day3.5_RA', 'day4_RA', 'day5_RA'))
  {
    # cc = 'day5_RA'
    cat(cc, '\n')
    cells = which(ref.combined$condition == cc)
    DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
            repel = TRUE, raster=FALSE, 
            #cols = c(cols_sel, cols_mouse),
            #order = c('day3_RA.rep1'),
            cells.highlight = cells,
            #shuffle = TRUE,
            cols.highlight = 'red'
    ) + NoLegend()
    
    ggsave(paste0(outDir, '/Integration_mNT_celltypes_noref.batchcorrection_highlight.', cc, '_2000HVGs.pdf'), 
           width = 16, height = 8)
    
  }
  
  #DimPlot(ref.combined, reduction = "umap")
  saveRDS(ref.combined, file = paste0(outDir, '/integrated_mNT_mouseGastrulation_SeuratRPCA.rds'))
  
  ref.combined = readRDS(file = paste0(outDir, '/integrated_mNT_mouseGastrulation_SeuratRPCA.rds'))
  
  ref.combined[['mnn.reconstructed']] = NULL

  ref.combined$celltype[which(is.na(ref.combined$celltype))] = 'unknown'
  
  DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE, split.by = 'dataset',
          repel = TRUE, raster=FALSE, cols = c(cols_sel, cols_mouse)) + NoLegend()
  
  ggsave(paste0(outDir, '/Integration_celltypes_split.dataset.pdf'), 
         width = 24, height = 8)
  
  DimPlot(ref.combined, reduction = "umap", group.by = "stage", label = TRUE,
          repel = TRUE, raster=FALSE)
  
  ggsave(paste0(outDir, '/Integration_stage.pdf'), 
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
  
  
  
  ##########################################
  # calculate the similarity distribution between mNT cells and cell types in the reference 
  ##########################################
  aa = FindVariableFeatures(aa, selection.method = "vst")
  aa <- ScaleData(aa, verbose = FALSE)
  aa <- RunPCA(aa, npcs = 50, verbose = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, reduction = "pca", dims = 1:30, n.neighbors = 50, 
                min.dist = 0.2) 
  DimPlot(aa, group.by = 'condition')
  
  aa <- FindNeighbors(aa, reduction = "pca", dims = 1:20)
  aa <- FindClusters(aa, resolution = 1.0)
  aa$clusters = aa$seurat_clusters
  
  DimPlot(aa, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)
  ggsave(paste0(outDir, '/umap_query_clusters.pdf'), 
         width = 12, height = 8)
  
  FeaturePlot(aa, features = 'Foxa2')
  
  source(paste0(functionDir, '/functions_dataIntegration.R'))
  
  ref = FindVariableFeatures(ref, selection.method = 'vst', nfeatures = 1000)
  
  cc = c(2, 1, 10, 15, 0, 4, 3, 12, 13)
  
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
# test projection method (slight different from integration)
# https://satijalab.org/seurat/articles/integration_mapping.html (original code)
# did not work well and don't know the reason
##########################################
Test_Seurat_projection = FALSE 
if(Test_Seurat_projection){
  mapping_method = "seurat_projection"
  
  aa = readRDS(file = paste0(RdataDir, 'seuratObject_mNT_selectedCondition_downsampled.1k.perCondition.rds'))
  
  ref = readRDS(file = paste0(RdataDir,  
                              'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v3.rds'))
  
  cols_mouse = sapply(ref$colour, function(x) {paste0('#', x, collapse = '')})
  names(cols_mouse) =ref$celltype
  cols_mouse = cols_mouse[match(unique(names(cols_mouse)), names(cols_mouse))]
  
  data_version = 'mapping_mNT.noRA.RA.d2_d5_Marioni2019_selectedCelltypes'
  
  
  features.common = intersect(rownames(aa), rownames(ref))
  
  aa = subset(aa, features = features.common)
  ref = subset(ref, features = features.common)
  
  aa$dataset = 'mNT'
  aa$stage = aa$condition
  aa$sequencing.batch = 'mNT'
  ref$dataset = 'ref'
  
  aa$celltype = paste0('mNT_', aa$condition)
  
  outDir = paste0(resDir,  data_version, '/', mapping_method, '/')
  system(paste0('mkdir -p ', outDir))
  
  
  ElbowPlot(ref, ndims = 50, reduction = 'pca')
  ref = RunUMAP(ref, reduction = "pca", dims = 1:30, n.neighbors = 30, 
                min.dist = 0.1, return.model = TRUE) 
  
  DimPlot(ref, reduction = "umap", 
          group.by = "celltype", label = TRUE,
          repel = TRUE, raster=FALSE, cols = cols_mouse) 
  
  anchors <- FindTransferAnchors(reference = ref, 
                                 query = aa, 
                                 dims = 1:50,
                                 normalization.method = "LogNormalize",
                                 reference.reduction = "pca",
                                 reduction = "rpca"
  )
  
  query <- MapQuery(anchorset = anchors, 
                    reference = ref, 
                    query = aa,
                    #refdata = list(celltype = "celltype"), 
                    reference.reduction = "pca", 
                    reduction.model = "umap")
  
  p1 <- DimPlot(ref, reduction = "umap", group.by = "celltype", 
                label = TRUE, label.size = 3,
                repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
  p2 <- DimPlot(query, reduction = "ref.umap", group.by = "condition", label = TRUE,
                label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
  p1 + p2
  
  
}

##########################################
# test seurat rpca but trying to use the orignal mnn reudction 
##########################################
Test_Seurat_rpca_bc = FALSE
if(Test_Seurat_rpca_bc){
  mapping_method = "seurat_rpca_ref.mnn"
  
  outDir = paste0(resDir,  data_version, '/', mapping_method, '/')
  system(paste0('mkdir -p ', outDir))
  
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
  
  
}
