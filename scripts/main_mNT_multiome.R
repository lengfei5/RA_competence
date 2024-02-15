##########################################################################
##########################################################################
# Project:
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Feb  6 10:49:04 2024
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R16597_mNT_10xmultiome_20240206'

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


########################################################
########################################################
# Section I : process the sample information and prepare for 10x cellranger
# 
########################################################
########################################################
dataDir = '/scratch/jiwang/Hannah_multiome/'

## process the sample info csv file
meta = read.csv2(file = paste0(dataDir, '/sample_info_parsed.csv'))
jj = unique(meta$sample_id)

meta = meta[match(jj, meta$sample_id), ]
meta$time = sapply(meta$sample_description, function(x){unlist(strsplit(as.character(x), ' '))[1]})
meta$time[which(meta$time == 'd2.0')] = 'd2'
meta$time[which(meta$time == 'd2.5')] = 'd2_12h'
meta$time[which(meta$time == 'd3.5')] = 'd3_12h'
meta$time[which(meta$time == 'd3.0')] = 'd3'
meta$time[which(meta$time == 'd4.0')] = 'd4'
meta$time[which(meta$time == 'd5.0')] = 'd5'

meta$treatment = 'RA'
meta$treatment[grep('beforeRA', meta$sample_description)] = 'beforeRA'
meta$treatment[grep('noRA', meta$sample_description)] = 'noRA'
meta$condition = paste0(meta$time, '_', meta$treatment)

write.csv2(meta, file = paste0(dataDir, 'sampleInfos.csv'), row.names = FALSE)

## save library file for each RNA+ATAC pair
meta$modality = 'RNA'
meta$modality[grep('ATAC', meta$sample_description)] = 'ATAC'

#meta$time = gsub('.', '_', meta$time)

cc = unique(meta$condition)

for(n in 1:length(cc))
{
  # n = 1
  cat(n, ' -- ', cc[n], '\n')
  libs = data.frame(matrix(NA, nrow = 2, ncol = 3), stringsAsFactors = FALSE)
  colnames(libs) = c('fastqs', 'sample', 'library_type')
  libs$library_type[1] = 'Chromatin Accessibility'
  libs$library_type[2] = 'Gene Expression'
  libs$sample[1] = meta$sample_id[which(meta$condition == cc[n] & meta$modality == 'ATAC')]
  libs$sample[2] = meta$sample_id[which(meta$condition == cc[n] & meta$modality == 'RNA')]
  libs$fastqs[1] = paste0('/scratch/jiwang/Hannah_multiome/raw_merged/', libs$sample[1])
  libs$fastqs[2] = paste0('/scratch/jiwang/Hannah_multiome/raw_merged/', libs$sample[2])
  write.csv(libs, file = paste0('/scratch/jiwang/Hannah_multiome/multiome_library/', 
                                'multiome_', cc[n], '.csv'), row.names = FALSE, quote = FALSE)
  
}

saveRDS(meta, file = paste0(RdataDir, 'meta_data.rds'))

########################################################
########################################################
# Section I : merge all cellranger peaks as peak consensus
# and quantify the count tables 
# 
# 
########################################################
########################################################
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)

# get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

meta = readRDS(paste0(RdataDir, 'meta_data.rds'))
dataDir = '../mNT_scmultiome_R16597'

##########################################
## first merge called peaks from different samples
##########################################
design = meta[which(meta$modality == 'ATAC'), ]

for(n in 1:nrow(design))
{
  # n = 1
  cat('----------- : ', n, ':',  design$condition[n], '-------------\n')
  
  topdir = paste0(dataDir, '/multiome_', design$condition[n], '/outs')
  
  p = read.table(paste0(topdir, "/atac_peaks.bed"), 
                 col.names = c("chr", "start", "end"))
  p = makeGRangesFromDataFrame(p)
  cat('---', length(p), ' peaks \n')
  ## atac_peaks.bed have the same peaks as in filtered_feature_bc_matrix.h5
  #counts <- Read10X_h5(paste0(topdir, "/filtered_feature_bc_matrix.h5"))
  
  if(n == 1){
    peaks = p
  }else{
    peaks = union(peaks, p)
  }
}

length(peaks)
combined.peaks = peaks
rm(peaks)

## a quick filtering the combined peaks
peakwidths = width(combined.peaks)
combined.peaks = combined.peaks[peakwidths > 50]

# there is a problem with coordinates starting at 0 for some reason...
combined.peaks = restrict(combined.peaks, start = 1)
cat(length(combined.peaks), ' combined peaks \n')

##########################################
# creat seurat object with combined peaks
##########################################
srat_cr = list()

for(n in 1:nrow(design))
#for(n in 2:nrow(design))
{
  # n = 1
  cat('----------- : ', n, ':',  design$condition[n], '-------------\n')
  
  # load nf output and process
  topdir = paste0(dataDir, '/multiome_', design$condition[n], '/outs')
  
  tic()
  counts <- Read10X_h5(filename = paste0(topdir, "/filtered_feature_bc_matrix.h5"))
  #counts <- Read10X_h5(filename = paste0(topdir, "/raw_feature_bc_matrix.h5"))
  fragpath <- paste0(topdir, "/atac_fragments.tsv.gz")
  toc()
  
  bb <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA"
  )
  
  cells_peak = colnames(counts$Peaks)
  cells_rna = colnames(counts$Peaks)
  cat(length(cells_peak), ' cells from atac \n')
  cat(ncol(bb), ' cell from rna \n')
  sum(!is.na(match(cells_rna, cells_peak)))
  
  #frags_l = CreateFragmentObject(path = fragpath, cells = colnames(counts$Peaks))
  frags_l = CreateFragmentObject(path = fragpath, cells = cells_peak)
  
  # slow step  mins without parall computation
  tic()
  feat = FeatureMatrix(fragments = frags_l, 
                       features = combined.peaks, 
                       cells = cells_peak,
                       process_n = 80000)
  
  saveRDS(feat, file = paste0(RdataDir, 'snATAC_FeatureMatrix_', design$condition[n], '.rds'))
  toc()
  
  # create ATAC assay and add it to the object
  bb[["ATAC"]] <- CreateChromatinAssay(
    counts = feat,
    sep = c(":", "-"),
    fragments = frags_l,
    annotation = annotation
  )
  
  # further subset the peaks to just focus on those at the standard chromosomes. 
  # This step is not a must but recommended as it makes some following analysis 
  # that requires calculation of local background GC content much easier.
  standard_chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
  idx_standard_chroms <- which(as.character(seqnames(granges(bb[['ATAC']]))) %in% standard_chroms)
  
  bb[["ATAC"]] <- subset(bb[["ATAC"]],
                             features = rownames(bb[["ATAC"]])[idx_standard_chroms])
  
  seqlevels(bb[['ATAC']]@ranges) <- intersect(seqlevels(granges(bb[['ATAC']])),
                                                  unique(seqnames(granges(bb[['ATAC']]))))
  
  bb$condition = design$condition[n]
  bb$time = design$time[n]
  bb$sampleID = design$sample_id[n]
  
  metadata <- read.csv(
    file = paste0(topdir, '/per_barcode_metrics.csv'),
    header = TRUE,
    row.names = 1
  )
  mm = match(colnames(bb), metadata$gex_barcode)
  metadata = metadata[mm, c(1,2, 3, 21:30)]
  bb = AddMetaData(bb, metadata = metadata)
  
  bb$pct_reads_in_peaks <- bb$atac_peak_region_fragments / bb$atac_fragments * 100
  bb$pct_usable_fragments = bb$atac_fragments/bb$atac_raw_reads
  
  bb = PercentageFeatureSet(bb, pattern = "^mt-",
                            col.name = "percent.mt", assay = "RNA")
  
  DefaultAssay(bb) <- "ATAC"
  
  bb = NucleosomeSignal(bb)
  bb = TSSEnrichment(bb, fast = TRUE)
  
  srat_cr[[n]] = bb
  
}

saveRDS(srat_cr, file = (paste0(RdataDir, 'seuratObj_scATAC_beforeMerged.peaks.cellranger_v1.rds')))

srat_cr = readRDS(file = paste0(RdataDir, 'seuratObj_scATAC_beforeMerged.peaks.cellranger_v1.rds'))
srat_reduced = Reduce(merge, srat_cr)

saveRDS(srat_reduced, file = (paste0(RdataDir, 'seuratObj_scATAC_merged.peaks.cellranger_v1.rds')))

########################################################
########################################################
# Section II : QCs of scATAC and snRNA
########################################################
########################################################
srat_cr = readRDS(file = paste0(RdataDir, 'seuratObj_scATAC_merged.peaks.cellranger_v1.rds'))

meta = readRDS(paste0(RdataDir, 'meta_data.rds'))
design = meta[which(meta$modality == 'ATAC'), ]
#design$condition = gsub('_12h_', '.5_', design$condition)
#design$condition = gsub('^d', 'day', design$condition)

srat_cr$condition = gsub('_12h_', '.5_', srat_cr$condition)
srat_cr$condition = gsub('^d', 'day', srat_cr$condition)

table(srat_cr$condition)


srat_cr$condition = factor(srat_cr$condition, levels = levels_sels)
Idents(srat_cr) = srat_cr$condition

make.QC.plots = FALSE
if(make.QC.plots){
  
  outDir = paste0(resDir, '/QCs_multiome')
  if(!dir.exists(outDir)) dir.create(outDir)
  
  #TSSPlot(srat_cr, group.by = 'high.tss') + NoLegend()
  
  VlnPlot(srat_cr, features = "nCount_RNA", ncol = 1, y.max = 10000, group.by = 'condition', pt.size = 0., 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(1000, 3000, 5000))
  
  ggsave(filename = paste0(outDir, '/QCs_nCount_RNA_cellRanger.pdf'), height =8, width = 12)
  
  VlnPlot(srat_cr, features = "nFeature_RNA", ncol = 1, y.max = 7000, group.by = 'condition', pt.size = 0., 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(500, 1000, 2000))
  
  ggsave(filename = paste0(outDir, '/QCs_nFeature_RNA_cellRanger.pdf'), height = 8, width = 12)
  
  VlnPlot(srat_cr, features = "nCount_ATAC", ncol = 1, y.max = 100000, group.by = 'condition', pt.size = 0., 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(1000, 5000, 10000))
  ggsave(filename = paste0(outDir, '/QCs_nCount_ATAC_cellRanger.pdf'), height =8, width = 12)
  
  
  VlnPlot(srat_cr, features = "nFeature_ATAC", ncol = 1, y.max = 50000, group.by = 'condition', pt.size = 0., 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(1000, 5000, 10000))
  
  ggsave(filename = paste0(outDir, '/QCs_nFeature_ATAC_cellRanger.pdf'), height =8, width = 12)
  
  
  VlnPlot(srat_cr, features = c("pct_reads_in_peaks"), raster=FALSE,  pt.size = 0.1)
  ggsave(filename = paste0(outDir, '/QCs_pct_readsWithinPeaks_cellRangerPeaks.pdf'), 
         height =8, width = 12 )
  
  VlnPlot(srat_cr, features = "percent.mt", ncol = 1, y.max = 100, group.by = 'condition', 
          pt.size = 0.0, 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(10, 20, 30))
  
  ggsave(filename = paste0(outDir, '/QCs_percent.mt_RNA_cellRanger.pdf'), height =8, width = 12)
  
  srat_cr$atac_percent.mt = srat_cr$atac_mitochondrial_reads/srat_cr$nCount_ATAC*100
  VlnPlot(srat_cr, features = "atac_percent.mt", ncol = 1, group.by = 'condition', 
          pt.size = 0.0, 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(10, 20, 30))
  
  
  VlnPlot(object = srat_cr, features = c("TSS.enrichment"), pt.size = 0, y.max = 20, raster=FALSE) +
    geom_hline(yintercept = c(1, 2))
  
  ggsave(filename = paste0(outDir, '/QCs_TSS.enrichment_cellRangerPeaks.pdf'), height =8, width = 12)
  
  VlnPlot(object = srat_cr, features = c("nucleosome_signal"), pt.size = 0, raster=FALSE) + 
  geom_hline(yintercept = c(2))
  
  ggsave(filename = paste0(outDir, '/QCs_nucleosome_signal_cellRangerPeaks.pdf'), height =8, width = 12)
  
  
}

# quick filtering 
srat_cr <- subset(
  x = srat_cr,
  subset = nFeature_RNA < 7500 &
    nFeature_RNA > 500 &
    percent.mt < 30 &
    nFeature_ATAC < 30000 &
    nFeature_ATAC > 1000 &
    TSS.enrichment > 1 &
    nucleosome_signal < 2
)

table(srat_cr$condition)

##########################################
# process snRNÅ-seq data
##########################################
DefaultAssay(srat_cr) <- "RNA"

srat_cr <- NormalizeData(srat_cr) %>%
  FindVariableFeatures(nfeatures = 5000) %>%
  CellCycleScoring(s.features = firstup(cc.genes.updated.2019$s.genes),
                   g2m.features = firstup(cc.genes.updated.2019$g2m.genes)) %>%
  ScaleData() %>%
  RunPCA(npcs = 50)

ElbowPlot(srat_cr, ndims = 30)

srat_cr <- RunUMAP(srat_cr, dims = 1:20, n.neighbors = 20, min.dist = 0.05, 
                   reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', cols = cols_sel) + 
  NoLegend()

p1 = DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', cols = cols_sel) + 
  NoLegend()
p2 <- FeaturePlot(srat_cr,
                  c("Foxa2","Pax6"),
                  reduction = "umap_rna", ncol = 1) 
p1 | p2

ggsave(filename = paste0(resDir, '/UMAP_RNA_20pcs_20neighbors_0.05dist_v1.pdf'), height =8, width = 12)

#DimPlot(srat_cr, label = TRUE, group.by = 'celltypes', repel = TRUE, reduction = 'umap') + NoLegend()

saveRDS(srat_cr, file = paste0(RdataDir, 
                               'seuratObj_multiome_snRNA.normalized.umap.rds'))

##########################################
# process and normalize the ATAC-seq data 
##########################################
srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.normalized.umap.rds'))

# normalize ATAC and UMAP
DefaultAssay(srat_cr) <- "ATAC"
srat_cr <- FindTopFeatures(srat_cr, min.cutoff = 'q5')
srat_cr <- FindTopFeatures(srat_cr, min.cutoff = 'q10')

srat_cr <- RunTFIDF(srat_cr, method = 1)
srat_cr <- RunSVD(srat_cr)

saveRDS(srat_cr, file = paste0(RdataDir,
                               'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized.rds'))

p1 <- ElbowPlot(srat_cr, ndims = 30, reduction="lsi")
p2 <- DepthCor(srat_cr, n = 30)
p1 | p2

ggsave(filename = paste0(resDir, '/lsi_Elbowplot_depthcor.pdf'), height =8, width = 12)

srat_cr <- RunUMAP(object = srat_cr, reduction = 'lsi', 
                   dims = 2:30, 
                   n.neighbors = 50, 
                   min.dist = 0.1, 
                   reduction.name = "umap_atac",
                   reduction.key = "UMAPATAC_"
                   )

DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_atac', cols = cols_sel) + NoLegend()

p1 <- DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_atac', cols = cols_sel) + NoLegend()
p2 <- FeaturePlot(srat_cr,
                  c("Foxa2","Pax6"),
                  reduction = "umap_atac", ncol = 1) 
p1 | p2

ggsave(filename = paste0(resDir, '/UMAP_ATAC_30pcs_50neighbors_0.1dist_v1.pdf'), height =8, width = 12)

srat_cr = FindNeighbors(object = srat_cr, reduction = 'lsi', dims = 2:30, 
                        force.recalc = T, graph.name = "thegraph")

srat_cr = FindClusters(object = srat_cr, verbose = FALSE, algorithm = 3, 
                       graph.name = "thegraph", resolution = 0.7)

DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_atac') + NoLegend()

ggsave(filename = paste0(resDir, '/UMAP_ATAC_quickClustering.pdf'), height =8, width = 12)


srat_cr = FindClusters(object = srat_cr, verbose = FALSE, algorithm = 3, 
                       graph.name = "thegraph", resolution = 1.0)

DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_atac') + NoLegend()

ggsave(filename = paste0(resDir, '/UMAP_ATAC_quickClustering_resolution.1.0.pdf'), height =8, width = 12)

saveRDS(srat_cr, file = paste0(RdataDir, 
                               'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized.umap.rds'))


##########################################
# link peaks and coveragePlots
##########################################
library(BSgenome.Mmusculus.UCSC.mm10)
library(presto)

DefaultAssay(srat_cr) = 'ATAC'

tic()
srat_cr <- RegionStats(srat_cr, genome = BSgenome.Mmusculus.UCSC.mm10)
toc()

# link peaks to genes
tic()
srat_cr <- LinkPeaks(
  object = srat_cr,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = c("Foxa2", "Pax6")
)
toc()

idents.plot <- c('0', '3', '2', '13', '1', '12', '9', '11', '6', '5', '10', '4', '7', '8')

CoveragePlot(
  object = srat_cr,
  region = "Foxa2",
  features = "Foxa2",
  expression.assay = "RNA",
  idents = idents.plot,
  #group.by = "celltype",
  extend.upstream = 20000,
  extend.downstream = 30000
)

ggsave(filename = paste0(resDir, '/FoxA2_enhancers_v1.pdf'), height = 12, width = 12)


CoveragePlot(
  object = srat_cr,
  region = "Pax6",
  features = "Pax6",
  expression.assay = "RNA",
  idents = idents.plot,
  #group.by = "celltype",
  extend.upstream = 20000,
  extend.downstream = 20000
)

ggsave(filename = paste0(resDir, '/Pax6_enhancers_v1.pdf'), height = 12, width = 12)


########################################################
########################################################
# Section III : symmetry breaking by subsetting only RA treatment
# 
########################################################
########################################################
Subset_RAtreated = FALSE
if(Subset_RAtreated){
  srat_cr = readRDS(file = paste0(RdataDir, 
                                  'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized.umap.rds'))
  
  Idents(srat_cr) = factor(srat_cr$condition)
  levels_sels = c('day2_beforeRA', 'day2.5_RA', 'day3_RA', 'day3.5_RA', 'day4_RA', 'day5_RA')
  srat_cr = subset(srat_cr, idents = levels_sels)
  
  outDir = paste0(resDir, '/RA_treated')
  if(!dir.exists(outDir)) dir.create(outDir)
  
  
  ## process RNA
  DefaultAssay(srat_cr) <- "RNA"
  
  srat_cr <- NormalizeData(srat_cr) %>%
    FindVariableFeatures(nfeatures = 5000) %>%
    CellCycleScoring(s.features = firstup(cc.genes.updated.2019$s.genes),
                     g2m.features = firstup(cc.genes.updated.2019$g2m.genes)) %>%
    ScaleData() %>%
    RunPCA(npcs = 50)
  
  ElbowPlot(srat_cr, ndims = 30)
  
  srat_cr <- RunUMAP(srat_cr, dims = 1:20, n.neighbors = 50, min.dist = 0.1, 
                     reduction.name = "umap_rna", reduction.key = "UMAPRNA_")
  
  DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', cols = cols_sel) + 
    NoLegend()
  
  p1 = DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', cols = cols_sel) + 
    NoLegend()
  p2 <- FeaturePlot(srat_cr,
                    c("Foxa2","Pax6"),
                    reduction = "umap_rna", ncol = 1) 
  p1 | p2
  
  ggsave(filename = paste0(outDir, '/UMAP_RNA_20pcs_50neighbors_0.1dist.pdf'), height = 8, width = 16)
  
  
  # normalize ATAC and UMAP
  DefaultAssay(srat_cr) <- "ATAC"
  #srat_cr <- FindTopFeatures(srat_cr, min.cutoff = 'q5')
  srat_cr <- FindTopFeatures(srat_cr, min.cutoff = 'q10')
  
  tic()
  srat_cr <- RunTFIDF(srat_cr, method = 1)
  toc()
  tic()
  srat_cr <- RunSVD(srat_cr)
  toc()
  
  
  p1 <- ElbowPlot(srat_cr, ndims = 30, reduction="lsi")
  p2 <- DepthCor(srat_cr, n = 30)
  p1 | p2
  
  ggsave(filename = paste0(outDir, '/lsi_Elbowplot_depthcor.pdf'), height =8, width = 12)
  
  srat_cr <- RunUMAP(object = srat_cr, reduction = 'lsi', 
                     dims = 2:20, 
                     n.neighbors = 50, 
                     min.dist = 0.1, 
                     reduction.name = "umap_atac",
                     reduction.key = "UMAPATAC_"
  )
  
  DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_atac', cols = cols_sel) + NoLegend()
  
  p1 <- DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_atac', group.by = 'condition',
                cols = cols_sel) + NoLegend()
  p2 <- FeaturePlot(srat_cr,
                    c("Foxa2","Pax6"),
                    reduction = "umap_atac", ncol = 1) 
  p1 | p2
  
  ggsave(filename = paste0(outDir, '/UMAP_ATAC_20pcs_50neighbors_0.1dist_v1.pdf'), height =8, width = 16)
  
  
  srat_cr = FindNeighbors(object = srat_cr, reduction = 'lsi', dims = 2:20, 
                          force.recalc = T, graph.name = "thegraph")
  
  srat_cr = FindClusters(object = srat_cr, verbose = FALSE, algorithm = 3, 
                         graph.name = "thegraph", resolution = 0.5)
  
  p3 = DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_atac') + NoLegend()
  
  p1|p2|p3
  
  ggsave(filename = paste0(outDir, '/UMAP_ATAC_quickClustering_resolution.0.5_FoxA2.Pax6.pdf'), 
         height =8, width = 24)
  
  
  saveRDS(srat_cr, file = paste0(outDir,
                                 '/seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized.rds'))
  
  ##########################################
  # link peaks and coveragePlots
  ##########################################
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(presto)
  
  DefaultAssay(srat_cr) = 'ATAC'
  
  tic()
  srat_cr <- RegionStats(srat_cr, genome = BSgenome.Mmusculus.UCSC.mm10)
  toc()
  
  # link peaks to genes
  tic()
  srat_cr <- LinkPeaks(
    object = srat_cr,
    peak.assay = "ATAC",
    expression.assay = "RNA",
    genes.use = c("Foxa2", "Pax6")
  )
  toc()
  
  
  idents.plot <- c('1', '3', '2', '7', '0', '6', '4', '5', '8')
  
  Idents(srat_cr) = factor(srat_cr$seurat_clusters, levels = idents.plot)
  
  CoveragePlot(
    object = srat_cr,
    region = "Foxa2",
    features = "Foxa2",
    expression.assay = "RNA",
    idents = idents.plot,
    #group.by = "celltype",
    extend.upstream = 50000,
    extend.downstream = 20000
  )
  
  ggsave(filename = paste0(outDir, '/FoxA2_enhancers_v1.pdf'), height = 10, width = 18)
  
  
  CoveragePlot(
    object = srat_cr,
    region = "Pax6",
    features = "Pax6",
    expression.assay = "RNA",
    idents = idents.plot,
    #group.by = "celltype",
    extend.upstream = 20000,
    extend.downstream = 20000
  )
  
  ggsave(filename = paste0(outDir, '/Pax6_enhancers_v1.pdf'), height = 10, width = 18)
  
}


