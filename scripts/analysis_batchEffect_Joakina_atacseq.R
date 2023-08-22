##########################################################################
##########################################################################
# Project: RA competence  
# Script purpose: check the batch effect of Joakina's atac-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Mar 28 15:12:46 2023
##########################################################################
##########################################################################
rm(list = ls())

RNA.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

# setup for data import and sequencing QCs
version.analysis = '_Joakina_atacseq_20230328'

resDir = paste0("../results/bulk_ATACseq_analysis_", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../Joaquina_atac/'

source('/groups/tanaka/People/current/jiwang/scripts/functions/Functions_atac.R')
source('/groups/tanaka/People/current/jiwang/scripts/functions/Functions_chipSeq.R')

require(ggplot2)
require(GenomicFeatures)


########################################################
########################################################
# Section I: import the sample design and make peak consensus 
# 
########################################################
########################################################
design = read.table(file = paste0(dataDir, 'filereport_read_run_PRJNA841970_tsv.txt'), 
                    sep = '\t', header = TRUE)

design = design[, c(1:4, 7,10,13:14)]
design = design[, c(1:5, 8)]
design$sample = sapply(design$sample_title, function(x){unlist(strsplit(as.character(x), ','))[1]})
design$replicate = sapply(design$sample_title, function(x){unlist(strsplit(as.character(x), ','))[2]})
design$replicate = gsub('Replicate', 'rep', design$replicate)

design = design[, c(4,7,8)]

colnames(design)[1:2] = c('SampleID', 'condition')

peakDir = paste0(dataDir,  'nf_out/peaks_macs2')
peak.files = list.files(path = peakDir,
                        pattern = '*_peaks.xls', full.names = TRUE)

index  = c()
for(n in 1:nrow(design))
{
  test = grep(design$SampleID[n], peak.files)
  if(length(test) != 1) {
    cat(length(test), 'peak files Found \n')
  }else{
    index = c(index, test)
  }
}

peak.files = peak.files[index]

# union of all peaks from all replicates
peak.merged = merge.peaks.macs2(peak.files, pcutoff = 10)

xx = peak.merged[grep('^chrUn_', seqnames(peak.merged), invert = TRUE)]

peak.merged = xx;
rm(xx)

## clean peaks
peaks = peak.merged
peaks = data.frame(peaks)
colnames(peaks)[c(1:3)] = c("chr", "start", "end")

dim(peaks)

peaks$peak.name = paste0(peaks$chr, ":", peaks$start, "_", peaks$end)

jj = match(unique(peaks$peak.name), peaks$peak.name)
peaks = peaks[jj, ];
dim(peaks)

peaks = peaks[which(peaks$chr != 'chrM'), ]

dim(peaks)


# preapre SAF input for featureCounts
require(Rsubread)

#counts = quantify.signals.within.peaks(peaks = peak.merged, bam.list=bam.list, rpkm.normalization = FALSE, isPairedEnd = FALSE)
SAF = data.frame(GeneID=peaks$peak.name, 
                 Chr=peaks$chr, 
                 Start=peaks$start, 
                 End=peaks$end, 
                 Strand=peaks$strand, stringsAsFactors = FALSE)

write.table(SAF, file = paste0(dataDir, 'merge_peak.saf'), sep = '\t', row.names = FALSE, 
            col.names = TRUE, quote = FALSE) 


saveRDS(design, file = paste0(RdataDir, 'design_sampleInfos.rds'))

########################################################
########################################################
# Section :
# # Run DESeq2 for QC (PCA) and the scaling factors for track normalization
# The consensus peaks were not filtered here
########################################################
########################################################
design = readRDS(file = paste0(RdataDir, 'design_sampleInfos.rds'))

RNA.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

xlist<-list.files(path=paste0(dataDir, 'featurecounts_peaks.Q30'),
                  pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge

all = cat.countTable(xlist, countsfrom = 'featureCounts')

colnames(design)[1] = 'SampleID'

counts = process.countTable(all=all, design = design[, c(1,2)])

save(design, counts, file = paste0(RdataDir, 'design_rawCounts.Rdata'))

##########################################
#  peak signal normalization
##########################################
load(file = paste0(RdataDir, 'design_rawCounts.Rdata'))

design$sampleID = design$SampleID

ss = apply(as.matrix(counts[, -1]), 1, mean)

par(mfrow=c(1,2))
hist(log10(as.matrix(counts[, -1])), breaks = 100, 
     xlab = 'log10(nb of reads within peaks)', 
     main = 'distribution')
plot(ecdf(log10(as.matrix(counts[, -1]) + 0.1)), 
     xlab = 'log10(nb of reads within peaks)', 
     main = 'cumulative distribution')

ss = apply(as.matrix(counts[, -1]), 2, sum)
design$usable.reads.withinPeaks = ss

ss = apply(as.matrix(counts[, -1]), 1, max)

cutoff = 20
hist(log10(ss), breaks = 200)
abline(v = log10(cutoff), col = 'red')
kk = which(ss>cutoff)
length(which(ss>cutoff))


require(ggplot2)
require(DESeq2)

rownames(counts) = counts$gene
dds <- DESeqDataSetFromMatrix(as.matrix(counts[kk, -1]), DataFrame(design), design = ~ condition)

#ss = rowSums(counts(dds))
#length(which(ss > quantile(ss, probs = 0.6)))
#dd0 = dds[ss > quantile(ss, probs = 0.6) , ]
dds = estimateSizeFactors(dds)

#sizefactors.UQ = sizeFactors(dd0)
#sizeFactors(dds) <- sizefactors.UQ
fpm = fpm(dds, robust = TRUE)

#dds <- estimateDispersions(dds, fitType = 'parametric')

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

save(design, vsd, file = paste0(RdataDir, 'design_vsd.Rdata'))

load(file = paste0(RdataDir, 'design_vsd.Rdata'))

pca=plotPCA(vsd, intgroup = colnames(design)[2], ntop = 3000, returnData = FALSE)
print(pca)


pca2save = as.data.frame(plotPCA(vsd, intgroup = colnames(design)[c(2)], returnData = TRUE, ntop = 3000))
#pca2save$name = paste0(design$condition, '_', design$SampleID, '_', design$batch)
pca2save$name = paste0(design$condition, '_', design$replicate)

pca2save$time = sapply(pca2save$group, function(x){unlist(strsplit(as.character(x), '_'))[1]}) 
pca2save$sag = as.numeric(sapply(pca2save$group, function(x){unlist(strsplit(as.character(x), '_'))[2]})) 
pca2save$celltype = sapply(pca2save$group, function(x){unlist(strsplit(as.character(x), '_'))[3]}) 
#pca2save$batch[grep('1361|1373', pca2save$name)] = 'new'
#pca2save$batch = as.factor(pca2save$batch)

ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color=time, shape = celltype, size = sag)) + 
  geom_point() + 
  geom_text(hjust = 0.2, nudge_y = 1.2, size=3)

plot(ggp)

ggsave(paste0(resDir, "/PCA_allatacseq_ntop3000_allSamples",  version.analysis, ".pdf"), 
       width = 16, height = 10)
