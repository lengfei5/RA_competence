##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: main script to analyze the mNT multiome data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Feb  6 10:49:04 2024
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
dataDir = '../mNT_scmultiome_R16597/'

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

write.csv2(meta, file = paste0(resDir, '/sampleInfos.csv'), row.names = FALSE)

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
dataDir = '../mNT_scmultiome_R16984'

design = meta[which(meta$modality == 'ATAC'), ]

##########################################
## first merge called peaks from different samples
##########################################

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
reload_processed = TRUE

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
  if(reload_processed){
    
    feat = readRDS(file = paste0(RdataDir, 'snATAC_FeatureMatrix_', design$condition[n], '.rds'))
  
  }else{
    tic()
    feat = FeatureMatrix(fragments = frags_l, 
                         features = combined.peaks, 
                         cells = cells_peak,
                         process_n = 80000)
    
    saveRDS(feat, file = paste0(RdataDir, 'snATAC_FeatureMatrix_', design$condition[n], '.rds'))
    toc()
    
  }
 
  
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
  
  saveRDS(bb, file = paste0(RdataDir, 'seuratObj_snATAC_', design$condition[n], '.rds'))
  
 
  
}

Merge_allSamples = FALSE
if(Merge_allSamples){
  srat_cr = list()
    
  for(n in 1:nrow(design))
  {
    cat('----------- : ', n, ':',  design$condition[n], '-------------\n')
    bb = readRDS(file = paste0(RdataDir, 'seuratObj_snATAC_', design$condition[n], '.rds'))
    
    srat_cr[[n]] = bb
    
  }
  
}

saveRDS(srat_cr, file = (paste0(RdataDir, 'seuratObj_scATAC_beforeMerged.peaks.cellranger_v1.rds')))

srat_cr = readRDS(file = paste0(RdataDir, 'seuratObj_scATAC_beforeMerged.peaks.cellranger_v1.rds'))
srat_cr = Reduce(merge, srat_cr)
saveRDS(srat_cr, file = (paste0(RdataDir, 'seuratObj_scATAC_merged.peaks.cellranger_v1.rds')))

########################################################
########################################################
# Section II : QCs of scATAC and snRNA
########################################################
########################################################
srat_cr = readRDS(file = paste0(RdataDir, 'seuratObj_scATAC_merged.peaks.cellranger_v1.rds'))

meta = readRDS(paste0(RdataDir, 'meta_data.rds'))
design = meta[which(meta$modality == 'ATAC'), ]

srat_cr$condition = design$condition[match(srat_cr$sampleID, design$sample_id)]
table(srat_cr$condition)
#design$condition = gsub('_12h_', '.5_', design$condition)
#design$condition = gsub('^d', 'day', design$condition)

srat_cr$condition = gsub('_12h_', '.5_', srat_cr$condition)
srat_cr$condition = gsub('^d', 'day', srat_cr$condition)

table(srat_cr$condition)

srat_cr$condition = factor(srat_cr$condition, levels = levels_sels)
Idents(srat_cr) = srat_cr$condition


make.QC.plots = TRUE
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
    geom_hline(yintercept = c(1000, 1500, 2000))
  
  ggsave(filename = paste0(outDir, '/QCs_nFeature_RNA_cellRanger.pdf'), height = 8, width = 12)
  
  VlnPlot(srat_cr, features = "nCount_ATAC", ncol = 1, y.max = 60000, group.by = 'condition', pt.size = 0., 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(1000, 10000, 20000))
  
  ggsave(filename = paste0(outDir, '/QCs_nCount_ATAC_cellRanger_toCompareFernando.pdf'), height =8, width = 12)
  
  VlnPlot(srat_cr, features = "nCount_ATAC", ncol = 1, y.max = 100000, group.by = 'condition', pt.size = 0., 
          log = FALSE, raster=FALSE) +
    geom_hline(yintercept = c(1000, 5000, 10000))
  
  ggsave(filename = paste0(outDir, '/QCs_nCount_ATAC_cellRanger_to1stBatch.pdf'), height =8, width = 12)
  
  
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
# process snRNA-seq data
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


p1 <- DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_rna', cols = cols_sel) + NoLegend()
p2 <- DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_atac', cols = cols_sel) + NoLegend()
p1 | p2

ggsave(filename = paste0(resDir, '/UMAP_RNA_vs_ATAC.pdf'), height =8, width = 12)


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



##########################################
# Identify the doublets  
##########################################
library(DoubletFinder)
require(Seurat)

srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized.umap.rds'))

DefaultAssay(srat_cr) <- "RNA"


run_DF = FALSE
if(run_DF){
  outDir = paste0(resDir, '/mNTs_multiome_DF_out')
  if(!dir.exists(outDir)) dir.create(outDir)
  
  Idents(srat_cr) = srat_cr$condition
  cc = unique(srat_cr$condition)
  srat_cr$DF_out = NA
  
  for(n in 1:length(cc))
  {
    # n = 1
    cat(n, ' -- ', cc[n], '\n')
    subs <- subset(srat_cr, condition == cc[n])
    
    subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 2000)
    subs <- ScaleData(subs)
    
    subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE)
    
    subs <- FindNeighbors(subs, dims = 1:30)
    subs <- FindClusters(subs, resolution = 1)
    
    subs <- RunUMAP(subs, dims = 1:30)
    
    sweep.res.list_nsclc <- paramSweep_v3(subs)
    
    sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
    bcmvn_nsclc <- find.pK(sweep.stats_nsclc)
    
    pK <- bcmvn_nsclc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
      filter(BCmetric == max(BCmetric)) %>%
      select(pK) 
    
    pK <- as.numeric(as.character(pK[[1]]))
    annotations <- subs@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations) 
    
    
    nExp_poi <- round(0.076*nrow(subs@meta.data))  
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    subs <- doubletFinder_v3(subs, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj,  
                             reuse.pANN = FALSE, sct = FALSE)
    
    df_out = subs@meta.data
    subs$DF_out = df_out[, grep('DF.classification', colnames(df_out))]
    
    DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'DF_out',
            raster=FALSE)
    ggsave(filename = paste0(outDir, '/subs_doubletFinder_out_', cc[n], '.pdf'), 
           width = 12, height = 8)
    
    saveRDS(subs, file = paste0(outDir, '/subs_doubletFinder_out_', cc[n], '.rds'))
    srat_cr$DF_out[match(colnames(subs), colnames(srat_cr))] = subs$DF_out
    
  }
  
  ## save the DF outputs 
  for(n in 1:length(cc))
  {
    # n = 1
    cat(n, '--', as.character(cc[n]), '\n')
    subs = readRDS(file = paste0(outDir, '/subs_doubletFinder_out_', cc[n], '.rds'))
    srat_cr$DF_out[match(colnames(subs), colnames(srat_cr))] = subs$DF_out
    
  }
  
  saveRDS(srat_cr, file = paste0(RdataDir, 
                       'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_DFout.rds'))
  
  
  srat_cr = readRDS(file = paste0(RdataDir, 
                                  'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_DFout.rds'))
  
  DefaultAssay(srat_cr) <- "RNA"
  
  DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', group.by = 'DF_out')
  
  srat_cr = subset(srat_cr, cells = colnames(srat_cr)[which(srat_cr$DF_out == 'Singlet')])
  
  srat_cr <- NormalizeData(srat_cr) %>%
    FindVariableFeatures(nfeatures = 5000) %>%
    ScaleData() %>%
    RunPCA(npcs = 100)
  
  ElbowPlot(srat_cr, ndims = 50)
  
  srat_cr <- RunUMAP(srat_cr, dims = 1:30, n.neighbors = 50, min.dist = 0.1, 
                     reduction.name = "umap_rna", reduction.key = "UMAPRNA_")
  
  DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', cols = cols_sel) + 
    NoLegend()
  
  
  DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', group.by = 'Phase')
  
  FeaturePlot(srat_cr, 
              features = c('nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'nFeature_ATAC', 'percent.mt'), 
              reduction = 'umap_rna')
  
  
  saveRDS(srat_cr, file = paste0(RdataDir, 
                                 'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                                 'DFoutSinglet.rds'))
  
                                
}

##########################################
# further filtering of cells due to UMAP structures
##########################################
srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                                'DFoutSinglet.rds'))

outDir = paste0(resDir, '/mNTs_multiome_processing')
if(!dir.exists(outDir)) dir.create(outDir)

p1 = DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', cols = cols_sel) + 
  NoLegend()

p2 = DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', group.by = 'Phase')

p1 + p2

ggsave(filename = paste0(outDir, '/mNTs_multiome_snRNA_condition_phase.pdf'), 
       width = 18, height = 8)

FeaturePlot(srat_cr, 
            features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 
                         'nCount_ATAC', 'nFeature_ATAC'), 
            reduction = 'umap_rna')

ggsave(filename = paste0(outDir, '/mNTs_multiome_snRNA_features4QCs.pdf'), 
       width = 16, height = 18)


p1 = FeatureScatter(srat_cr, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", group.by = 'condition', 
               cols = cols_sel)
p2 = FeatureScatter(srat_cr, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = 'condition', 
                    cols = cols_sel)
p3 = FeatureScatter(srat_cr, feature1 = "nFeature_RNA", feature2 = "nCount_ATAC", group.by = 'condition', 
                    cols = cols_sel)
p4 = FeatureScatter(srat_cr, feature1 = "nFeature_RNA", feature2 = "nFeature_ATAC", group.by = 'condition', 
                    cols = cols_sel)

(p1 + p2)/(p3+p4)


VlnPlot(srat_cr, features = "nCount_RNA", ncol = 1, y.max = 20000, group.by = 'condition', pt.size = 0., 
        log = FALSE, raster=FALSE) +
  geom_hline(yintercept = c(2000, 3000, 5000, 15000, 12000))


VlnPlot(srat_cr, features = "nFeature_RNA", ncol = 1, y.max = 7000, group.by = 'condition', pt.size = 0., 
        log = FALSE, raster=FALSE) +
  geom_hline(yintercept = c(1000, 1500, 2000, 5000))

VlnPlot(srat_cr, features = "percent.mt", ncol = 1, y.max = 35, 
             group.by = 'condition', pt.size = 0.0, 
             log = FALSE, raster=FALSE)


# quick filtering 
srat_cr <- subset(
  x = srat_cr,
  subset = nCount_RNA > 2000 &
    nCount_RNA < 15000 &
    nFeature_RNA < 5000 &
    nFeature_RNA > 1500 &
    percent.mt < 30 
)


## remove Rp and mt
head(grep('^Rp[sl]|^mt-', rownames(srat_cr))) # double check if ribo and mito genes 
srat_cr = subset(srat_cr, features = rownames(srat_cr)[grep('^Rp[sl]|^mt-', rownames(srat_cr), 
                                                            invert = TRUE)])

srat_cr <- NormalizeData(srat_cr) %>%
  FindVariableFeatures(nfeatures = 5000) %>%
  ScaleData() %>%
  RunPCA(npcs = 100)
ElbowPlot(srat_cr, ndims = 50)

srat_cr <- RunUMAP(srat_cr, dims = 1:50, n.neighbors = 100, min.dist = 0.1, 
                   reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', cols = cols_sel) + 
  NoLegend()


saveRDS(srat_cr, file = paste0(RdataDir, 
                               'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                               'DFoutSinglet_cell.gene.filtered.rds'))


##########################################
# try to regress out 'nFeature_RNA',  "percent.mt"   
##########################################
srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                                'DFoutSinglet_cell.gene.filtered.rds'))

p1 = FeatureScatter(srat_cr, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", group.by = 'condition', 
                    cols = cols_sel)
p2 = FeatureScatter(srat_cr, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = 'condition', 
                    cols = cols_sel)

p1 + p2

VlnPlot(srat_cr, features = "nCount_RNA", ncol = 1, y.max = 20000, group.by = 'condition', pt.size = 0., 
        log = FALSE, raster=FALSE) +
  geom_hline(yintercept = c(2000, 3000, 5000, 15000, 12000))


VlnPlot(srat_cr, features = "nFeature_RNA", ncol = 1, y.max = 7000, group.by = 'condition', pt.size = 0., 
        log = FALSE, raster=FALSE) +
  geom_hline(yintercept = c(1000, 1500, 2000, 5000))

VlnPlot(srat_cr, features = "percent.mt", ncol = 1, y.max = 35, 
        group.by = 'condition', pt.size = 0.0, 
        log = FALSE, raster=FALSE)


srat_cr = NormalizeData(srat_cr, normalization.method = "LogNormalize", scale.factor = 10000)

srat_cr <- FindVariableFeatures(srat_cr, selection.method = "vst", nfeatures = 5000)

all.genes <- rownames(srat_cr)
srat_cr <- ScaleData(srat_cr, features = all.genes, vars.to.regress = c('nFeature_RNA',  "percent.mt"))


saveRDS(srat_cr, file = paste0(RdataDir, 
                               'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                               'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT.rds'))


srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                                'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT.rds'))


srat_cr <- FindVariableFeatures(srat_cr, nfeatures = 5000, selection.method = "vst") %>%
  RunPCA(npcs = 100)
ElbowPlot(srat_cr, ndims = 50)

srat_cr <- RunUMAP(srat_cr, dims = 1:30, n.neighbors = 50, min.dist = 0.1, 
                   reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', cols = cols_sel) + 
  NoLegend()

FeaturePlot(srat_cr, features = c('Foxa2', 'Pax6'))


##########################################
# try to identify some weird clusters and discard them 
##########################################
srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                                'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT.rds'))

outDir = paste0(resDir, '/mNTs_multiome_processing')
if(!dir.exists(outDir)) dir.create(outDir)


srat_cr <- FindVariableFeatures(srat_cr, nfeatures = 5000, selection.method = "vst") %>%
  RunPCA(npcs = 100)
ElbowPlot(srat_cr, ndims = 50)

srat_cr <- RunUMAP(srat_cr, dims = 1:30, n.neighbors = 50, min.dist = 0.1, 
                   reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', cols = cols_sel) + 
  NoLegend()

DimPlot(srat_cr, label = FALSE, repel = TRUE, reduction = 'umap_rna', group.by = 'Phase')


srat_cr <- FindVariableFeatures(srat_cr, nfeatures = 3000, selection.method = "vst") %>%
  RunPCA(npcs = 50)
srat_cr <- FindNeighbors(srat_cr, dims = 1:30)

srat_cr <- FindClusters(srat_cr, verbose = FALSE, algorithm = 3, resolution = 0.7)

p1 = DimPlot(srat_cr, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(srat_cr, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

ggsave(filename = paste0(outDir, '/mNTs_multiome_snRNA_seuratClusters.res0.7.pdf'), 
       width = 14, height = 6)



## subsetting only RA
aa = subset(srat_cr, cells = colnames(srat_cr)[grep('_RA|_beforeRA', srat_cr$condition)])
aa$condition = droplevels(aa$condition)

aa <- FindVariableFeatures(aa, nfeatures = 3000, selection.method = "vst") %>%
  RunPCA(npcs = 50)

ElbowPlot(aa, ndims = 30)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.1, 
                   reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

aa$clusters = aa$seurat_clusters

all.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.4)
all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(aa, features = top10$gene) + NoLegend()

ggsave(filename = paste0(outDir, '/mNTs_multiome_snRNA_seuratClusters.res0.7_heatmap_markerGenes.pdf'), 
       width = 14, height = 20)


all.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

## cluster 12: Tubb3
saveRDS(aa, file = paste0(RdataDir, 
                               'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                               'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT_RAsamples.rds'))


aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                           'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT_RAsamples.rds'))


cells_drops = colnames(aa)[which(!is.na(match(aa$clusters, c('12', '13', '14', '15'))))]


### subset noRA samples
rm(aa)

aa = subset(srat_cr, cells = colnames(srat_cr)[grep('_noRA|_beforeRA', srat_cr$condition)])
aa$condition = droplevels(aa$condition)

aa <- FindVariableFeatures(aa, nfeatures = 3000, selection.method = "vst") %>%
  RunPCA(npcs = 50)

ElbowPlot(aa, ndims = 30)

aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1, 
              reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

aa$clusters = aa$seurat_clusters

all.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.4)
all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(aa, features = top10$gene) + NoLegend()

ggsave(filename = paste0(outDir, '/mNTs_multiome_snRNA_seuratClusters.res0.7_noRA_heatmap_markerGenes.pdf'), 
       width = 14, height = 20)


all.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

## cluster 12: Tubb3
saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                          'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT_noRAsamples.rds'))

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                           'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT_noRAsamples.rds'))

cells_drops = c(cells_drops, colnames(aa)[which(!is.na(match(aa$clusters, c('15', '16'))))])

saveRDS(cells_drops, file = paste0(outDir, 'smallClusters_toDrop_RA_noRA.rds'))


## drop cells in small clusteres of RA or noRA samples
srat_cr = subset(srat_cr, cells = colnames(srat_cr)[which(is.na(match(colnames(srat_cr), cells_drops)))])

srat_cr <- FindVariableFeatures(srat_cr, nfeatures = 5000, selection.method = "vst") %>%
  RunPCA(npcs = 50, weight.by.var = FALSE)
ElbowPlot(srat_cr, ndims = 50)

srat_cr <- RunUMAP(srat_cr, dims = 1:30, n.neighbors = 50, min.dist = 0.1, 
                   reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', cols = cols_sel) + 
  NoLegend()

DimPlot(srat_cr, label = FALSE, repel = TRUE, reduction = 'umap_rna', group.by = 'Phase')


aa = subset(srat_cr, cells = colnames(srat_cr)[grep('_RA|_beforeRA', srat_cr$condition)])
aa$condition = droplevels(aa$condition)

aa <- FindVariableFeatures(aa, nfeatures = 3000, selection.method = "vst") 

aa = RunPCA(aa, weight.by.var = FALSE)

ElbowPlot(aa, ndims = 30)

aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 30, min.dist = 0.1, 
              reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

DimPlot(aa, label = FALSE, repel = TRUE, reduction = 'umap_rna', group.by = 'Phase')

saveRDS(srat_cr, file = paste0(RdataDir, 
                               'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                               'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT_',
                               'discardTinyClusters.RAnoRA.rds'))

##########################################
# test snRNA umap parameters
##########################################
#srat_cr = readRDS(file = paste0(RdataDir, 
#              'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_DFout.rds'))

srat_cr = readRDS(paste0(RdataDir, 
       'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
       'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT_',
       'discardTinyClusters.RAnoRA.rds'))

DefaultAssay(srat_cr) <- "RNA"

Explore.umap.parameters = FALSE
if(Explore.umap.parameters){
  source(paste0(functionDir, '/functions_scRNAseq.R'))
  explore.umap.params.combination(sub.obj = srat_cr, 
                                  resDir = resDir, 
                                  pdfname = 'mNTs_multiome_snRNA_umap_afterFiltering_weight.by.var.TRUE.pdf',
                                  use.parallelization = FALSE,
                                  group.by = 'condition',
                                  cols = cols_sel, 
                                  weight.by.var = TRUE,
                                  nfeatures.sampling = c(3000, 5000, 8000),
                                  nb.pcs.sampling = c(30, 50), 
                                  n.neighbors.sampling = c(30, 50, 100),
                                  min.dist.sampling = c(0.1, 0.3)
                                  
  )
  
  
}


ElbowPlot(srat_cr, ndims = 50)

srat_cr <- RunUMAP(srat_cr, dims = 1:20, n.neighbors = 20, min.dist = 0.05, 
                   reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap_rna', cols = cols_sel) + 
  NoLegend()



########################################################
########################################################
# Section II: analysis of RA samples
# 
########################################################
########################################################
# srat_cr = readRDS(paste0(RdataDir, 
#                          'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
#                          'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT_',
#                          'discardTinyClusters.RAnoRA.rds'))
# 
# DefaultAssay(srat_cr) <- "RNA"
# aa = subset(srat_cr, cells = colnames(srat_cr)[grep('_RA', srat_cr$condition)])

rm(srat_cr)

outDir = paste0(resDir, '/mNTs_multiome_RA_processing')
if(!dir.exists(outDir)) dir.create(outDir)


aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                           'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT_RAsamples.rds'))

aa$condition = droplevels(aa$condition)

aa <- FindVariableFeatures(aa, nfeatures = 3000, selection.method = "vst") %>%
  RunPCA(npcs = 50, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 30)

aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1, 
              reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

DimPlot(aa, label = FALSE, repel = TRUE, group.by = 'Phase', raster=FALSE)

#aa <- FindNeighbors(aa, dims = 1:30)
#aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
#p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
#p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
#p1 + p2
#aa$clusters = aa$seurat_clusters

load(file = '../results/Rdata/tfs_sps_geneExamples_4scRNAseq.Rdata')
knownGenes =  c('Pou5f1', 'Sox2', 'Zfp42', 'Utf1',
                'Otx2', 'Cyp26a1', 'Stra8',
                'Hoxa1', 'Hoxa3', 'Hoxb4',
                'Sox1', 'Pax6', 'Tubb3', 'Elavl4', 'Neurod4',
                'Foxa2', 'Shh', 'Arx', 'Vtn', "Spon1", 'Slit2', "Ntn1",
                'Nkx2-2', 'Olig2', 'Pax3', 'Pax7', 'Nkx6-1')

Discard_cellCycle.corrrelatedGenes = TRUE
if(Discard_cellCycle.corrrelatedGenes){
  library(scater)
  Idents(aa) = aa$condition
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  
  # Identifying the likely cell cycle genes between phases,
  # using an arbitrary threshold of 5%.
  scaledMatrix = GetAssayData(aa, slot = c("scale.data"))
  
  diff <- getVarianceExplained(scaledMatrix, data.frame(phase = aa$Phase))
  diff = data.frame(diff, gene = rownames(diff))
  diff = diff[order(-diff$phase), ]
  
  hist(diff$phase, breaks = 100); 
  abline(v = c(1:5), col = 'red')
  
  rm(scaledMatrix)
  
  genes_discard = diff$gene[which(diff$phase > 5)]
  cat(length(genes_discard), 'genes to discard \n')
  
  tfs_sels = intersect(genes_discard, gene_examples)
  print(tfs_sels)
  
  knownGenes_sels = intersect(genes_discard, knownGenes)
  print(knownGenes_sels)
  
  if(length(knownGenes_sels)>0) genes_discard = setdiff(genes_discard, knownGenes_sels)
  
  tfs_sels = intersect(genes_discard, gene_examples)
  print(tfs_sels)
  
  aa = subset(aa, features = setdiff(rownames(aa), genes_discard))
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000) # find subset-specific HVGs
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 50, min.dist = 0.1, 
                reduction.name = "umap_rna", reduction.key = "UMAPRNA_")
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  
  FeaturePlot(aa, features = c('Foxa2', 'Pax6'))
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  
  aa <- FindNeighbors(aa, dims = 1:20)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)
  
  ggsave(filename = paste0(outDir, '/mNTs_multiome_snRNA_RAsamples_rmCellCycleGenes_',
                           'seuratClusters.res0.7.pdf'), 
         width = 14, height = 8)
  
}

aa$clusters = aa$seurat_clusters

all.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(aa, features = top10$gene) + NoLegend()

ggsave(filename = paste0(outDir, '/mNTs_multiome_snRNA_seuratClusters.res0.7_RA_heatmap_markerGenes.pdf'), 
       width = 14, height = 20)


aa = subset(aa, cells = colnames(aa)[which(aa$clusters != '12' & aa$clusters != '14')])
aa$clusters = droplevels(aa$clusters)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000) # find subset-specific HVGs

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)

ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 50, min.dist = 0.1, 
              reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)


saveRDS(aa, file = paste0(RdataDir, 
                           'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized_',
                           'DFoutSinglet_cell.gene.filtered_regressed.nFeature.pctMT_RAsamples',
                           '_rmCCgenes_rmMatureNeurons.rds'))




########################################################
########################################################
# Section III : Integrate the scRNA-seq data and multiome
# 
########################################################
########################################################
Integrate_scRNAseq_Multiome = FALSE
if(Integrate_scRNAseq_Multiome){
  
  outDir = paste0(resDir, '/Integration_scRNAseq_multiome/')
  system(paste0('mkdir -p ', outDir))
  
  
  ## multiome
  srat_cr = readRDS(file = paste0(RdataDir, 
                'seuratObj_multiome_snRNA.normalized.umap_scATAC.normalized.rds'))
  
  DefaultAssay(srat_cr) = 'RNA'
  srat_cr[['ATAC']] = NULL
  
  srat_cr$chemistry = 'multiome'
  
  ## scRNA-seq
  aa =  readRDS(file = paste0('../results/Rdata/',  
                              'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                              'cellCycleScoring_annot.v1_', 'mNT_scRNAseq',
                              '_R13547_10x_mNT_20220813', '.rds'))
  
  aa$chemistry = 'scRNA'
  
  features.common = intersect(rownames(aa), rownames(srat_cr))
  
  aa = subset(aa, features = features.common)
  srat_cr = subset(srat_cr, features = features.common)
  
  multi_combine = merge(aa, y = srat_cr, add.cell.ids = c("scRNA", ""), project = "RA_competence")
  
  saveRDS(multi_combine, file = paste0(outDir, 
                                 'seuratObj_multiome_scRNAseq_combined.rds'))
  
  ## test the harmony for integration
  multi_combine = readRDS(file = paste0(outDir, 
                                        'seuratObj_multiome_scRNAseq_combined.rds'))
  
  table(multi_combine$chemistry)
  
  multi_combine <- NormalizeData(multi_combine)
  multi_combine <- FindVariableFeatures(multi_combine, nfeatures = 5000)
  
  VariableFeaturePlot(object = multi_combine)
  
  tic()
  multi_combine <- ScaleData(multi_combine, vars.to.regress = c("nCount_RNA"))
  toc()
  
  saveRDS(multi_combine, file = paste0(outDir, 
                                       'seuratObj_multiome_scRNAseq_combined_regress.nCountRNA_pca.rds'))
  
  
  ## select only the overlapping conditions
  Idents(multi_combine) = multi_combine$condition
  multi_combine = subset(multi_combine, idents = c("day0_beforeRA", "day1_beforeRA", "day3_RA.rep2", 
                                          "day6_noRA", "day6_RA", "day5_noRA"), invert = TRUE)
  
  multi_combine$chem = multi_combine$chemistry
  
  multi_combine <- FindVariableFeatures(multi_combine, nfeatures = 5000)
  multi_combine <- RunPCA(multi_combine, features = VariableFeatures(object = multi_combine))
  
  main_pc_chem <- DimPlot(object = multi_combine, reduction = 'pca', group.by = 'chem', 
                          cols = c("turquoise4", "gray"), raster=FALSE)
  
  main_pc <- FeaturePlot(object = multi_combine, reduction = 'pca', 
                         features = c('Foxa2', 'Pax6', 'Shh', 'Sox17'))
  #main_pc1 <- FeaturePlot(object = multi_combine, reduction = 'pca', 
  #                        features = c('SLC17A6', 'GAD1'))
  main_pc_chem
  
  #Run harmony to correct for chemistry
  library(harmony)
  tic()
  multi_combine <- RunHarmony(multi_combine, 
                              "chem",
                              #nclust = 10, 
                              max.iter.harmony = 10, 
                              epsilon.harmony = -Inf,
                              verbose = TRUE,
                              reference_values = 'scRNA',
                              plot_convergence = TRUE)
  toc()
  
  #Identify the highest contributing PCs
  ElbowPlot(multi_combine, ndims = 50)
  
  
  #multi_combine <- FindNeighbors(multi_combine, dims = 1:npcs, reduction = 'harmony')
  
  #Use clustree to find stability of res for clustering
  # Set different resolutions 
  
  #res.used <- seq(0.1,1,by=0.2)
  # Loop over and perform clustering of different resolutions 
  #for(i in res.used){
  #  multi_combine <- FindClusters(multi_combine, resolution = i)}
  # Make Plot
  #clustree(multi_combine, layout="sugiyama") + theme(legend.position = "bottom") + scale_color_brewer(palette = "Set1") + scale_edge_color_continuous(low = "grey80", high = "red")
  
  #Find clusters 
  #res = 0.7
  #multi_combine <- FindClusters(multi_combine, resolution = res)
  npcs = 30
  multi_combine <- RunUMAP(multi_combine, reduction = 'harmony', dims = 1:npcs, 
                           reduction.name = 'umap_harmony')
  
  DimPlot(multi_combine, reduction = "umap_harmony", label = TRUE, group.by = 'condition', cols = cols_sel)
  
  DimPlot(multi_combine, reduction = "umap_harmony", label = TRUE, group.by = 'chem')
  
  
  
  #Make UMAPs
  dimplot_main_clus <- DimPlot(multi_combine, reduction = "umap_harmony", label = T, group.by = 'seurat_clusters') +
    NoLegend()
  dimplot_main_sample <- DimPlot(multi_combine, reduction = "umap_harmony", label = F, group.by = 'condition')
  
  dimplot_main_chem <- DimPlot(multi_combine, reduction = "umap_harmony", label = F, group.by = 'sample', 
                               cols = c("turquoise4", "turquoise3", "gray", "gray", "gray", "gray", "gray", "gray"), 
                               shuffle = T)
  dimplot_main_clus
  dimplot_main_sample
  dimplot_main_chem
  
  
}

########################################################
########################################################
# Section IV : symmetry breaking by subsetting only RA treatment
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



