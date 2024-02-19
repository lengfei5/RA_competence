rm(list = ls())

version.analysis = '_R16597_mNT_10xmultiome_20240206'

resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/'

#source(paste0(functionDir,  'functions_scATAC.R'))
#source(paste0(functionDir, 'functions_scRNAseq.R'))
#source(paste0(functionDir, 'functions_Visium.R'))

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

options(future.globals.maxSize = 160 * 1024^3)
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


########################################################
########################################################
# Section I : process the sample information and prepare for 10x cellranger
# 
########################################################
########################################################
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
                       process_n = 100000)
  
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

