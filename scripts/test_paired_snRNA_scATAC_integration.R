##########################################################################
##########################################################################
# Project: RA project 
# Script purpose: test integration methods for paired scRNA and scATAC
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon May 18 15:36:11 2026
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R16597_mNT_10xmultiome_reseq_20240517'

resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/'

source(paste0(functionDir,  'functions_scATAC.R'))
source(paste0(functionDir, 'functions_scRNAseq.R'))
source(paste0(functionDir, 'functions_Visium.R'))

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
require(SeuratObject)
library(data.table)
library(DropletUtils)
library(edgeR)
library(future)
library(tictoc)

library(RColorBrewer)
require(Matrix)
library(RaceID)
require(tictoc)
require(gridExtra)
require(cowplot)

options(future.globals.maxSize = 200 * 1024^3)
set.seed(1234)
mem_used()

species = 'mNT_multiome'

### Parameters for this project
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


## subset the mutliome conditions
levels_sels = c("day2_beforeRA",
                "day2.5_RA", "day3_RA", "day3.5_RA",  "day4_RA", "day5_RA", 
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA")

cc = names(cols)
cc[which(cc == 'day3_RA.rep1')] = 'day3_RA'
names(cols) = cc

cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/mNTs_multiome_RA_scRNAscATAC_integraton_test')
if(!dir.exists(outDir)) dir.create(outDir)


##########################################
# features  
##########################################
sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

chromfactors = readRDS(file = '../data/annotations/human_chromatin_remodelers_Epifactors.database.rds')
chromfactors = as.character(unlist(sapply(chromfactors, firstup)))

features_1 = unique(c('Sox2', 'Sox1', 'Tubb3', 'Elavl3', 
                      'Irx3', 'Irx5', 'Pax3', 'Pax7',
                      'Pax6', 'Olig2', 'Nkx2-9', 'Nkx2-2', 
                      'Nkx6-1', 'Foxa2', 'Arx', 'Shh'
)) # DV overview

features_2 = c('Dhrs3', 'Rarg', 'Cyp26a1',
               'Pou3f1', 'Hoxa1', 'Gas1', 'Spry4', 'Sox11', 
               'Cdh1', 'Cdh2', 'Shh', 'Rfx4', 'Zfp42', 'Tcf15', 'Prrx2', 'Gdf3',
               'Etv5', 'Fgf4', 'Otx2', 'Zscan10', 'Apoe', 'Peg10', 'Klf9', 'Tshz1', 'Skil', 'Zfp703')

features = unique(c(c('Pax6', 'Foxa2', 'Sox1', 'Sox2', 'Tubb3', 'Shh', 'Arx',
                      'Zfp703', 'Lef1', 'Irx5', 'Pou5f1', 'Otx2', 'Adgra2', 'Hoxb4', 
                      'Nkx2-2', 'Nkx2-9', 'Nkx6-1', 'Olig2', 'Pax3', 'Pax7', 'Cyp26a1', 'Dhrs3'), 
                    features_1,
                    features_2)) # marker wanted by Hannah

gene_examples = unique(c('Foxa2', 'Pax6', c('Zfp42', 'Tcf15', 'Skil', 'Lef1',
                                            'Sox2', 'Pou5f1', 'Sall4', 'Tdgf1', # pluripotency markers
                                            'Nanog', 'Nr5a2', #'Prdm14', 
                                            'Klf4', 'Fgf4', 'Esrrb', 'Tcf3', 'Tbx3'), # naive pluripotency
                         c('Zfp42', 'Tcf15', 'Skil',
                           'Fgf5', 'Otx2', 'Pou3f1', 'Lef1', 'Dnmt3b', 'Dnmt3a',	
                           'Foxd3', 'Utf1', 'Tcf15', 'Zic3', 'Rhox5', 'Etv5', 'Etv4',	
                           'Lin28b', 'Sox4', 'Sox3', 'Sox11'
                         ),
                         c('Lhx1','Eomes', 'Sox2', 'Hoxb4', 'Hoxb5', 'Hoxb6','Zfp703'),
                         c('Zfp42', 'Tcf15', 'Skil', 'Lef1', 'Dhrs3', 'Rarg', 'Cyp26a1'),
                         features,
                         features_1,
                         features_2
                         
))


# sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))

########################################################
########################################################
# Section I : check gene examples of Pluripotent, RA response, AP identity genes
# test subclustering for each time points
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir,
                           'seuratObj_RAmultiome_snRNA.normalized.umap_scATAC.normalized',
                           'noMatureNeurons.smallClusters_wnn.umap.clusters.rds'))


##########################################
# check again the data QC 
##########################################
check.QC.plots = FALSE
if(check.QC.plots.QC.plots){
  
  #TSSPlot(aa, group.by = 'high.tss') + NoLegend()
  
  VlnPlot(aa, features = "nCount_RNA", ncol = 1, y.max = 12000, group.by = 'condition', pt.size = 0., 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(1000, 2000, 5000))
  
  ggsave(filename = paste0(outDir, '/QCs_nCount_RNA_cellRanger.pdf'), height =8, width = 12)
  
  VlnPlot(aa, features = "nFeature_RNA", ncol = 1, y.max = 5000, group.by = 'condition', pt.size = 0., 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(1000, 1500, 2000))
  
  ggsave(filename = paste0(outDir, '/QCs_nFeature_RNA_cellRanger.pdf'), height = 8, width = 12)
  
  VlnPlot(aa, features = "nCount_ATAC", ncol = 1, y.max = 60000, group.by = 'condition', pt.size = 0., 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(1500, 10000))
  
  ggsave(filename = paste0(outDir, '/QCs_nCount_ATAC_cellRanger_toCompareFernando.pdf'), height =8, width = 12)
  
  VlnPlot(aa, features = "nCount_ATAC", ncol = 1, y.max = 100000, group.by = 'condition', pt.size = 0., 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(1000, 5000, 10000))
  
  ggsave(filename = paste0(outDir, '/QCs_nCount_ATAC_cellRanger_to1stBatch.pdf'), height =8, width = 12)
  
  
  VlnPlot(aa, features = "nFeature_ATAC", ncol = 1, y.max = 50000, group.by = 'condition', pt.size = 0., 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(1000, 5000, 10000))
  
  ggsave(filename = paste0(outDir, '/QCs_nFeature_ATAC_cellRanger.pdf'), height =8, width = 12)
  
  
  VlnPlot(aa, features = c("pct_reads_in_peaks"), raster=FALSE, group.by = 'condition',  pt.size = 0.1)
  ggsave(filename = paste0(outDir, '/QCs_pct_readsWithinPeaks_cellRangerPeaks.pdf'), 
         height =8, width = 12 )
  
  VlnPlot(aa, features = "percent.mt", ncol = 1, y.max = 100, group.by = 'condition', 
          pt.size = 0.0, 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(10, 20, 30))
  
  ggsave(filename = paste0(outDir, '/QCs_percent.mt_RNA_cellRanger.pdf'), height =8, width = 12)
  
  aa$atac_percent.mt = aa$atac_mitochondrial_reads/aa$nCount_ATAC*100
  VlnPlot(aa, features = "atac_percent.mt", ncol = 1, group.by = 'condition', 
          pt.size = 0.0, 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(10, 20, 30))
  
  
  VlnPlot(object = aa, features = c("TSS.enrichment"), pt.size = 0, y.max = 20, raster=FALSE) +
    geom_hline(yintercept = c(1, 2))
  
  ggsave(filename = paste0(outDir, '/QCs_TSS.enrichment_cellRangerPeaks.pdf'), height =8, width = 12)
  
  VlnPlot(object = aa, features = c("nucleosome_signal"), pt.size = 0, raster=FALSE) + 
    geom_hline(yintercept = c(2))
  
  ggsave(filename = paste0(outDir, '/QCs_nucleosome_signal_cellRangerPeaks.pdf'), height =8, width = 12)
  
  
  
  DefaultAssay(aa) <- "ATAC"
  FindTopFeatures(aa, min.cutoff = 100)
  srat_cr <- FindTopFeatures(aa, min.cutoff = 'q50')
  
  
  # quick filtering 
  aa <- subset(
    x = aa,
    subset = nFeature_RNA < 7500 &
      nFeature_RNA > 500 &
      percent.mt < 30 &
      nFeature_ATAC < 30000 &
      nFeature_ATAC > 1000 &
      TSS.enrichment > 1 &
      nucleosome_signal < 2
  )
  
  
}

DensityScatter(aa, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

VlnPlot(
  object = aa, 
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4, group.by = 'condition',
  pt.size = 0
)



table(aa$condition)

DimPlot(aa, reduction = 'umap_rna', group.by = 'condition', cols = cols_sel, label = TRUE, repel = TRUE)


##########################################
# rerun the Seurat WNN
##########################################
DefaultAssay(aa) <- "RNA"
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 100, min.dist = 0.1, 
              reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)


DefaultAssay(aa) <- "ATAC"

#aa <- FindTopFeatures(aa, min.cutoff = 100)
#aa <- FindTopFeatures(aa, min.cutoff = 'q10')
aa <- FindTopFeatures(aa, min.cutoff = 500) ## stringently select top peaks > 500 cells
aa <- RunTFIDF(aa, method = 1)
aa <- RunSVD(aa)

p1 <- ElbowPlot(aa, ndims = 30, reduction="lsi")
p2 <- DepthCor(aa, n = 30)
p1 | p2

aa <- RunUMAP(object = aa, reduction = 'lsi', 
              dims = 2:30, 
              n.neighbors = 30, 
              min.dist = 0.1, 
              reduction.name = "umap_atac",
              reduction.key = "UMAPATAC_"
)

DimPlot(object = aa, label = TRUE, reduction = 'umap_atac', group.by = 'condition') 

aa <- FindMultiModalNeighbors(aa, 
                              reduction.list = list("pca", "lsi"), 
                              dims.list = list(1:30, 2:30)
                              )

aa <- RunUMAP(aa, nn.name = "weighted.nn", 
              reduction.name = "wnn.umap", 
              reduction.key = "wnnUMAP_", min.dist = 0.1)

p1 <- DimPlot(aa, reduction = "umap_rna", group.by = "condition", 
              label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(aa, reduction = "umap_atac", group.by = "condition", label = TRUE, 
              label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(aa, reduction = "wnn.umap", group.by = "condition", label = TRUE, 
              label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3  & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_RAsamples_v2.pdf'), 
       width = 16, height = 6)


saveRDS(aa, file = paste0(RdataDir,
                          'seuratObj_RAmultiome_snRNA.normalized.umap_scATAC.normalized',
                          'noMatureNeurons.smallClusters_wnn.umap.v2.rds'))



##########################################
# add gene activities matrix 
##########################################
aa = readRDS(file = paste0(RdataDir,
                           'seuratObj_RAmultiome_snRNA.normalized.umap_scATAC.normalized',
                           'noMatureNeurons.smallClusters_wnn.umap.v2.rds'))

DefaultAssay(aa) <- "ATAC"

calculate_geneActivityMatrix = FALSE
if(calculate_geneActivityMatrix) {
  gene.activities <- GeneActivity(aa, 
                                  extend.upstream = 2000,
                                  extend.downstream = 0,
                                  biotypes = "protein_coding")
  
  saveRDS(gene.activities, file = paste0(RdataDir,
                                         'seuratObj_RAmultiome_snRNA.normalized.umap_scATAC.normalized',
                                         'noMatureNeurons.smallClusters_geneActivityMatrix.rds'))
  
}else{
  gene.activities = readRDS(file = paste0(RdataDir,
                                          'seuratObj_RAmultiome_snRNA.normalized.umap_scATAC.normalized',
                                          'noMatureNeurons.smallClusters_geneActivityMatrix.rds'))
  
}


# add the gene activity matrix to the Seurat object as a new assay
aa[['geneActivity']] <- CreateAssayObject(counts = gene.activities)

aa <- NormalizeData(
  object = aa,
  assay = 'geneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(aa$nCount_geneActivity)
  
)

DefaultAssay(aa) <- 'geneActivity'
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
aa = ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.1, 
              reduction.name = "umap_geneActivity", reduction.key = "UMAPGeneActivity_")

FeaturePlot(aa, features = c('nCount_geneActivity', 'nFeature_geneActivity'))

VlnPlot(object = aa, features = c("nCount_geneActivity"), pt.size = 0, y.max = 70000, raster=FALSE) +
  geom_hline(yintercept = c(1000, 2000))

ggsave(filename = paste0(outDir, '/QCs_nCount_geneActivityMatrix.pdf'), height =8, width = 12)

VlnPlot(object = aa, features = c("nFeature_geneActivity"), group.by = 'condition',
        pt.size = 0, y.max = 15000, raster=FALSE) +
  geom_hline(yintercept = c(1000, 2000, 3000))

ggsave(filename = paste0(outDir, '/QCs_nFeature_geneActivityMatrix.pdf'), height =8, width = 12)

FeaturePlot(aa, features = c('nCount_geneActivity', 'nFeature_geneActivity'))

# quick filtering 
aa <- subset(
  x = aa,
  subset = nFeature_geneActivity > 3000 &
    nFeature_geneActivity < 16000
)


aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
aa = ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.1, 
              reduction.name = "umap_geneActivity", reduction.key = "UMAPGeneActivity_")

DimPlot(aa, group.by = 'condition', label = TRUE, 
        label.size = 2.5, repel = TRUE)

FeaturePlot(aa, features = c('nCount_geneActivity', 'nFeature_geneActivity'))


aa = FindNeighbors(object = aa, reduction = 'pca', dims = 1:20, 
                        force.recalc = T, graph.name = "geneactivitygraph")

aa = FindClusters(object = aa, verbose = FALSE, algorithm = 3, 
                       graph.name = "geneactivitygraph", resolution = 1.0)

p1 = DimPlot(aa, label = TRUE, label.size = 4, repel = TRUE)
p2 = DimPlot(aa, label = TRUE, label.size = 4, repel = TRUE, group.by = 'condition')

p1 + p2 

FeaturePlot(aa, features = c('nFeature_geneActivity'))

VlnPlot(aa, features = 'nFeature_geneActivity')

aa <- subset(x = aa, subset = geneactivitygraph_res.1 != 9)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
aa = ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1, 
              reduction.name = "umap_geneActivity", reduction.key = "UMAPGeneActivity_")

p1 = DimPlot(aa, group.by = 'condition', label = TRUE, 
        label.size = 2.5, repel = TRUE)

p2 = DimPlot(aa, group.by = 'Phase', label.size = 2.5, repel = TRUE)

p1 + p2

ggsave(filename = paste0(outDir, '/RAsamples_geneActivityMatrix_umap_phases.pdf'), height =6, width = 14)


## filter also RNA genes
ggs = rownames(aa)

DefaultAssay(aa) = 'RNA'
mm = match(ggs, rownames(aa))

saveRDS(aa, file = paste0(RdataDir,
                           'seuratObj_RAmultiome_snRNA.normalized.umap_scATAC.normalized',
                           'noMatureNeurons.smallClusters_wnn.umap.v2_geneActivityMatrix.rds'))



##########################################
# test MOFA 
# original code from https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/10x_scRNA_scATAC.html
##########################################
library(data.table)
library(ggplot2)
library(Seurat)
library(Signac)
# 
# # for GSEA analysis
# #library(msigdbr)
# 
# # For motif enrichment analysis
# #library(JASPAR2020)
# #library(TFBSTools)
# #library(BSgenome.Hsapiens.UCSC.hg38)
# 
# MOFA
library(MOFA2)

aa = readRDS(file = paste0(RdataDir,
              'seuratObj_RAmultiome_snRNA.normalized.umap_scATAC.normalized',
              'noMatureNeurons.smallClusters_wnn.umap.v2_geneActivityMatrix.rds'))


DefaultAssay(aa) = 'RNA'
aa <- FindVariableFeatures(aa,
                               selection.method = "vst",
                               nfeatures = 5000,
                               assay = "RNA",
                               verbose = FALSE
)

DefaultAssay(aa) = 'ATAC'

aa <- FindTopFeatures(aa, assay='ATAC', min.cutoff = 2000)
print(length(aa[['ATAC']]@var.features))
# 
# mofa <- create_mofa_from_Seurat(aa, assays = c("RNA","ATAC"), groups = 'condition',
#                                 extract_metadata = FALSE)
# mofa

## manually make mofa object from Seurat 
## original code from https://github.com/bioFAM/MOFA2/blob/master/R/create_mofa.R
.split_seurat_into_groups <- function(seurat, groups, assay = "RNA", layer = "data", features = NULL) {
  data <- SeuratObject::GetAssayData(object = seurat, assay = assay, layer = layer)
  if(is.null(data) | any(dim(data) == 0)){
    stop(paste("No data present in the layer",layer, "of the assay",assay ,"in the Seurat object."))
  }
  if (!is.null(features)) data <- data[features, , drop=FALSE]
  .split_data_into_groups(list(data), groups)[[1]]
}
.split_data_into_groups <- function(data, groups) {
  group_indices <- split(seq_along(groups), factor(groups, exclude = character(0))) # factor call avoids dropping NA
  lapply(data, function(x) {
    lapply(group_indices, function(idx) {
      x[, idx, drop = FALSE]
    })
  })
}

.create_samples_metadata <- function(object) {
  # TO-DO: CHECK SAMPLE AND GROUP COLUMN IN PROVIDED METADATA
  foo <- lapply(object@data[[1]], colnames)
  tmp <- data.frame(
    sample = unname(unlist(foo)),
    group = unlist(lapply(names(foo), function(x) rep(x, length(foo[[x]])) )),
    stringsAsFactors = FALSE
  )
  if (.hasSlot(object, "samples_metadata") && (length(object@samples_metadata) > 0)) {
    object@samples_metadata <- cbind(tmp, object@samples_metadata[match(tmp$sample, rownames(object@samples_metadata)),, drop = FALSE])
  } else {
    object@samples_metadata <- tmp
  }
  return(object)
}

.create_features_metadata <- function(object) {
  tmp <- data.frame(
    feature = unname(unlist(lapply(object@data, function(x) rownames(x[[1]])))),
    view = unlist(lapply(seq_len(object@dimensions$M), function(x) rep(views_names(object)[[x]], object@dimensions$D[[x]]) )),
    stringsAsFactors = FALSE
  )
  if (.hasSlot(object, "features_metadata") && (length(object@features_metadata) > 0)) {
    object@features_metadata <- cbind(tmp, object@features_metadata[match(tmp$feature, rownames(object@features_metadata)),])
  } else {
    object@features_metadata <- tmp
  }
  return(object)
}

.rename_duplicated_features <- function(object) {
  feature_names <- unname(unlist(lapply(object@data, function(x) rownames(x[[1]]))))
  duplicated_names <- unique(feature_names[duplicated(feature_names)])
  if (length(duplicated_names)>0) 
    warning("There are duplicated features names across different views. We will add the suffix *_view* only for those features 
            Example: if you have both TP53 in mRNA and mutation data it will be renamed to TP53_mRNA, TP53_mutation")
  # Rename data matrices
  for (m in names(object@data)) {
    for (g in names(object@data[[m]])) {
      tmp <- which(rownames(object@data[[m]][[g]]) %in% duplicated_names)
      if (length(tmp)>0) {
        rownames(object@data[[m]][[g]])[tmp] <- paste(rownames(object@data[[m]][[g]])[tmp], m, sep="_")
      }
    }
  }
  
  # Rename features metadata
  tmp <- object@features_metadata[["feature"]] %in% duplicated_names
  object@features_metadata[tmp,"feature"] <- paste(object@features_metadata[tmp,"feature"], object@features_metadata[tmp,"view"], sep="_")
  return(object)
}

.quality_control <- function(object, verbose = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Check views names
  if (verbose == TRUE) message("Checking views names...")
  stopifnot(!is.null(views_names(object)))
  stopifnot(!duplicated(views_names(object)))
  if (any(grepl("/", views_names(object)))) {
    stop("Some of the views names contain `/` symbol, which is not supported.
  This can be fixed e.g. with:
    views_names(object) <- gsub(\"/\", \"-\", views_names(object))")
  }
  
  # Check groups names
  if (verbose == TRUE) message("Checking groups names...")
  if (any(grepl("/", groups_names(object)))) {
    stop("Some of the groups names contain `/` symbol, which is not supported.
    This can be fixed e.g. with:
    groups_names(object) <- gsub(\"/\", \"-\", groups_names(object))")
  }
  stopifnot(!is.null(groups_names(object)))
  stopifnot(!duplicated(groups_names(object)))
  
  # Check samples names
  if (verbose == TRUE) message("Checking samples names...")
  stopifnot(!is.null(samples_names(object)))
  stopifnot(!duplicated(unlist(samples_names(object))))
  enc <- stringi::stri_enc_mark(unlist(samples_names(object)))
  if (any(enc!="ASCII")) {
    tmp <- unname(unlist(samples_names(object))[enc!="ASCII"])
    stop(sprintf("non-ascii characters detected in the following samples names, please rename them and run again create_mofa():\n- %s ", paste(tmp, collapse="\n- ")))
    print()
  }
  
  # Check features names
  if (verbose == TRUE) message("Checking features names...")
  stopifnot(!is.null(features_names(object)))
  stopifnot(!duplicated(unlist(features_names(object))))
  enc <- stringi::stri_enc_mark(unlist(features_names(object)))
  if (any(enc!="ASCII")) {
    tmp <- unname(unlist(features_names(object))[enc!="ASCII"])
    stop(sprintf("non-ascii characters detected in the following features names, please rename them and run again create_mofa():\n- %s ", paste(tmp, collapse="\n- ")))
    print()
  }
  
  # Check dimensionalities in the input data
  if (verbose == TRUE) message("Checking dimensions...")
  N <- object@dimensions$N
  D <- object@dimensions$D
  for (i in views_names(object)) {
    for (j in groups_names(object)) {
      stopifnot(ncol(object@data[[i]][[j]]) == N[[j]])
      stopifnot(nrow(object@data[[i]][[j]]) == D[[i]])
      stopifnot(length(colnames(object@data[[i]][[j]])) == N[[j]])
      stopifnot(length(rownames(object@data[[i]][[j]])) == D[[i]])
    }
  }
  
  # Check that there are no features with complete missing values (across all groups)
  if (object@status == "untrained" || object@data_options[["loaded"]]) {
    if (verbose == TRUE) message("Checking there are no features with complete missing values...")
    for (i in views_names(object)) {
      if (!(is(object@data[[i]][[1]], "dgCMatrix") || is(object@data[[i]][[1]], "dgTMatrix"))) {
        tmp <- as.data.frame(sapply(object@data[[i]], function(x) rowMeans(is.na(x)), simplify = TRUE))
        if (any(unlist(apply(tmp, 1, function(x) mean(x==1)))==1))
          warning("You have features which do not contain a single observation in any group, consider removing them...")
      }
    }
  }
  
  # check dimensionalities of sample_covariates 
  if (verbose == TRUE) message("Checking sample covariates...")
  if(.hasSlot(object, "covariates") && !is.null(object@covariates)){
    stopifnot(ncol(object@covariates) == sum(object@dimensions$N))
    stopifnot(nrow(object@covariates) == object@dimensions$C)
    stopifnot(all(unlist(samples_names(object)) == colnames(object@covariates)))
  }
  
  # Sanity checks that are exclusive for an untrained model  
  if (object@status == "untrained") {
    
    # Check features names
    if (verbose == TRUE) message("Checking features names...")
    tmp <- lapply(object@data, function(x) unique(lapply(x,rownames)))
    for (x in tmp) stopifnot(length(x)==1)
    for (x in tmp) if (any(duplicated(x[[1]]))) stop("There are duplicated features names within the same view. Please rename")
    all_names <- unname(unlist(tmp))
    duplicated_names <- unique(all_names[duplicated(all_names)])
    if (length(duplicated_names)>0) 
      warning("There are duplicated features names across different views. We will add the suffix *_view* only for those features 
            Example: if you have both TP53 in mRNA and mutation data it will be renamed to TP53_mRNA, TP53_mutation")
    for (i in names(object@data)) {
      for (j in names(object@data[[i]])) {
        tmp <- which(rownames(object@data[[i]][[j]]) %in% duplicated_names)
        if (length(tmp)>0) {
          rownames(object@data[[i]][[j]])[tmp] <- paste(rownames(object@data[[i]][[j]])[tmp], i, sep="_")
        }
      }
    }
    
    # Sanity checks that are exclusive for a trained model  
  } else if (object@status == "trained") {
    # Check expectations
    if (verbose == TRUE) message("Checking expectations...")
    stopifnot(all(c("W", "Z") %in% names(object@expectations)))
    # if(.hasSlot(object, "covariates") && !is.null(object@covariates)) stopifnot("Sigma" %in% names(object@expectations))
    stopifnot(all(sapply(object@expectations$W, is.matrix)))
    stopifnot(all(sapply(object@expectations$Z, is.matrix)))
    
    # Check for intercept factors
    if (object@data_options[["loaded"]]) { 
      if (verbose == TRUE) message("Checking for intercept factors...")
      if (!is.null(object@data)) {
        factors <- do.call("rbind",get_factors(object))
        r <- suppressWarnings( t(do.call('rbind', lapply(object@data, function(x) 
          abs(cor(colMeans(do.call("cbind",x),na.rm=TRUE),factors, use="pairwise.complete.obs"))
        ))) )
        intercept_factors <- which(rowSums(r>0.75)>0)
        if (length(intercept_factors)) {
          warning(sprintf("Factor(s) %s are strongly correlated with the average expression of features for at least one of your omics. Such factors appear when there are differences in the total 'levels' between your samples, *sometimes* because of poor normalisation in the preprocessing steps.\n",paste(intercept_factors,collapse=", ")))
        }
      }
    }
    
    # Check for correlated factors
    if (verbose == TRUE) message("Checking for highly correlated factors...")
    Z <- do.call("rbind",get_factors(object))
    op <- options(warn=-1) # suppress warnings
    
    noise <- matrix(rnorm(n=length(Z), mean=0, sd=1e-10), nrow(Z), ncol(Z))
    tmp <- cor(Z+noise); diag(tmp) <- NA
    options(op) # activate warnings again
    if (max(tmp,na.rm=TRUE)>0.5) {
      warning("The model contains highly correlated factors (see `plot_factor_cor(MOFAobject)`). \nWe recommend that you train the model with less factors and that you let it train for a longer time.\n")
    }
    
  }
  
  return(object)  
}


groups = 'condition'
assays = c("RNA","ATAC")
layer = "data"

groups <- aa@meta.data[, groups]
features <- lapply(assays, function(i) aa@assays[[i]]@var.features)
names(features) <- assays

# Extract data matrices
data_matrices <- lapply(assays, function(i) 
  .split_seurat_into_groups(aa, groups = groups, assay = i, layer = layer, 
                            features = features[[i]]))
names(data_matrices) <- assays

# Initialise MOFA object
object <- new("MOFA")
object@status <- "untrained"
object@data <- data_matrices

# Define dimensions
object@dimensions[["M"]] <- length(assays)
object@dimensions[["D"]] <- vapply(data_matrices, function(m) nrow(m[[1]]), 1L)
object@dimensions[["G"]] <- length(data_matrices[[1]])
object@dimensions[["N"]] <- vapply(data_matrices[[1]], function(g) ncol(g), 1L)
object@dimensions[["K"]] <- 0


# Set views & groups names
groups_names(object) <- as.character(names(data_matrices[[1]]))
views_names(object)  <- assays

object@samples_metadata <- aa@meta.data


# Create sample metadata
object <- .create_samples_metadata(object)

# Create features metadata
object <- .create_features_metadata(object)

# Rename duplicated features
object <- .rename_duplicated_features(object)

# Do quality control
object <- .quality_control(object)

mofa = object 
rm(object)

model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15

mofa <- prepare_mofa(mofa,
                     model_options = model_opts
)

mofa <- run_mofa(mofa, use_basilisk = TRUE)

saveRDS(mofa, file = paste0(RdataDir,
                                       'seuratObj_RAmultiome_snRNA.normalized.umap_scATAC.normalized',
                                       'noMatureNeurons.smallClusters_mofa.rds'))



