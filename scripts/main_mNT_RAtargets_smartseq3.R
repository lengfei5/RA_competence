##########################################################################
##########################################################################
# Project: RA competence
# Script purpose: analyze the RA receptor targets
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Sat Jun 19 16:34:59 2021
##########################################################################
##########################################################################
rm(list = ls())
RNA.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

# setup for data import and sequencing QCs
version.analysis = '_L118404_smartseq3_20221117'

resDir = paste0("../results/RA_targets", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../L118404_smartseq3/'

Collect.QCs.stat = TRUE
Counts.to.Use = "UMI"

########################################################
########################################################
# Section I : data processing and sequencing quality controls
# 
########################################################
########################################################

##########################################
# prepare design matrix, count table and QC table
##########################################
design = openxlsx::read.xlsx(paste0(dataDir, '20221112_sequencingInformation_hannahs.xlsx'), sheet = 1)
design = design[, c(2,4, 10, 1)]
design$unid = design$File 
design$unid = gsub('_R1.fastq.gz', '', design$unid)
design$unid = gsub('_R2.fastq.gz', '', design$unid)
design = design[match(unique(design$unid), design$unid), ]
design = design[, -c(4)]

xx = openxlsx::read.xlsx(paste0(dataDir, '20221112_sequencingInformation_hannahs.xlsx'), sheet = 3)
xx = xx[, c(1, 11)]
mm = match(design$SampleID, xx$SampleID)
design = data.frame(design, xx[mm, ], stringsAsFactors = FALSE)

## discard samples with issues
kk = which(!is.na(match(design$LibraryID, c('L118404', 'L118421', 'L118444'))))
kk = kk[which(design$FC_Code[kk] == 'HNN5VDSX3')]
design = design[-kk, ]

## clean the conditions
design$condition = as.character(design$Notes)
design$rep = sapply(design$condition, function(x) {unlist(strsplit(as.character(x), '-'))[2]}) 
design$time = sapply(design$condition, function(x) {unlist(strsplit(as.character(x), '_'))[1]}) 
design$treat = sapply(design$condition, function(x) {unlist(strsplit(as.character(x), '_'))[2]}) 
design$time[which(design$time == '48h')] = 'd2'
design$time[which(design$time == '72h')] = 'd3'
design$time[which(design$time == '96h')] = 'd4'
design$time[which(design$time == '48plus18h')] = 'd2.18h'
design$time[which(design$time == '72plus18h')] = 'd3.18h'
design$time[which(design$time == '96plus18h')] = 'd4.18h'

design$treat = gsub('plus', '.', design$treat)

design$condition = paste0(design$treat, '_', design$time)

design = design[, -c(5)]
design$unid = sapply(design$unid, function(x) {unlist(strsplit(as.character(x), '-'))[2]}) 

saveRDS(design, file = paste0(RdataDir, 'desing_sampleInfo.rds'))

if(Collect.QCs.stat){
  MultiQC.Dir = paste0(dataDir, 'multiqc_data/')
  
  stats = read.delim(paste0(MultiQC.Dir, 'multiqc_general_stats.txt'), sep = '\t', header = TRUE)
  alignment = read.delim(paste0(MultiQC.Dir, 'multiqc_star.txt'), sep = '\t')
  
  ##########################################
  # process the stat table
  ##########################################
  keep = data.frame(matrix(NA, nrow = nrow(design), ncol = 11))
  colnames(keep) = c('sampleID', 'pct.duplication', 'pct.GC', 'avg.seq.length', 'total.reads', 
                     'pct.assign', 'assigned.reads', 'alignment.rate', 'trimmed.reads', 'unique.aligned', 'multimapper')
  keep$sampleID = design$sampleID
  
  stats = data.frame(stats, stringsAsFactors = FALSE)
  alignment = data.frame(alignment, stringsAsFactors = FALSE)
  
  for(n in 1:nrow(keep))
  {
    jj = grep(keep$sampleID[n], stats$Sample)
    jj = jj[grep('R1', stats$Sample[jj])]
    jj1 = jj[grep('R1_trimmed_sorted_umiDedup', stats$Sample[jj], invert = TRUE)]
    jj2 = jj[grep('R1_trimmed_sorted_umiDedup', stats$Sample[jj])]
    keep[n, c(2, 3, 4, 5)] = stats[jj1, c(2, 3, 4, 6)]
    keep[n, c(6, 7)] = stats[jj2, c(9, 10)]
    jj3 = grep(keep$sampleID[n], alignment$Sample)
    keep$alignment.rate[n] = alignment$uniquely_mapped_percent[jj3]
    keep$trimmed.reads[n] = alignment$total_reads[jj3]
    keep$unique.aligned[n] = alignment$uniquely_mapped[jj3]
    keep$multimapper[n] = alignment$multimapped[jj3]
    
  }
  
  xx = data.frame(design, keep[, -1], stringsAsFactors = FALSE)
  #xx = xx[order(xx$fileName), ]
  
  write.csv(xx, file = paste0(resDir, '/sampleInfos_QCs.stats.csv'), row.names = FALSE)
  
  design = xx
  
  saveRDS(design, file = paste0(RdataDir, 'sampleInfo_QC.stats.rds'))
  
}

##################################################
## Import UMI count table
##################################################
source(RNA.functions)
source(RNA.QC.functions)

Import.onlyUMI = TRUE

design = readRDS(file = paste0(RdataDir, 'desing_sampleInfo.rds'))
colnames(design)[1] = 'sampleID'
colnames(design)[4] = 'SampleID'

dataDir = '../L118404_smartseq3/featurecounts_Q30'
xlist = list.files(path=dataDir,
                   pattern = "*featureCounts.txt$", full.names = TRUE) ## list of data set to merge

all = cat.countTable(xlist, countsfrom = 'featureCounts')

counts = process.countTable(all=all, design = design[, c(4, 6)], merge.technicalRep.sameID = FALSE)

save(design, counts, file=paste0(RdataDir, 'design_rawCounts', version.analysis, '.Rdata'))

##########################################
# gene names converted from ensID to gene symbol 
##########################################
load(file=paste0(RdataDir, 'design_rawCounts', version.analysis, '.Rdata'))

colnames(counts)[-1] = paste0(design$condition, '_', design$LibraryID, '_', design$SampleID)

## convert gene names
annot = read.delim(paste0('/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/', 
                          'ens_BioMart_GRCm38.p6.txt'), sep = '\t', header = TRUE)

# annot = annot[which(annot$Gene.type == 'protein_coding'), ]
annot = annot[!is.na(match(annot$Chromosome.scaffold.name, as.character(c(1:19, 'X', 'Y')))), ]

mm = match(counts$gene, annot$Gene.stable.ID)
cat(length(which(!is.na(mm))), ' genes in the count table \n')

counts = counts[!is.na(mm), ]

rownames(counts) = counts$gene

counts$gene = annot$Gene.name[match(counts$gene, annot$Gene.stable.ID)]

counts = counts[!is.na(counts$gene), ]
gg.uniq = unique(counts$gene)
counts = counts[match(gg.uniq, counts$gene), ]
rownames(counts) = counts$gene
counts = counts[, -1]

Add.more.sample.details = FALSE
if(Add.more.sample.details){
  library("openxlsx")
  xx = read.xlsx(paste0(dataDir, 'Detailed_sample_information_with_SeqIDs.xlsx'), sheet = 1)
  xx = data.frame(design, xx[match(design$SampleID, xx$Seq.ID), ], stringsAsFactors = FALSE)
  design = xx
  
  colnames(design)[grep('Triplicate', colnames(design))] = 'triplicate.nb'
  
  xx = read.csv(paste0(dataDir, '20210429_Venus_F02-Batch_Analysis_29042021122425.csv'))
  saveRDS(xx, file = paste0(RdataDir, 'facs_positive_negative_ratios.rds'))
  
}

save(design, counts, file = paste0(RdataDir, 'design_rawCounts_geneSymbols', version.analysis, '.Rdata'))

##########################################
# QCs of replicates and conditions
##########################################
load(file = paste0(RdataDir, 'design_rawCounts_geneSymbols', version.analysis, '.Rdata'))

QC.for.cpm = FALSE
if(QC.for.cpm){
  
  source(RNA.functions)
  source(RNA.QC.functions)
  
  raw = as.matrix(counts)
  
  #kk = which(design$SampleID != '161040' & design$condition != 'N2B27')
  #raw = raw[, -kk]
  
  ss = apply(as.matrix(raw), 1, sum)
  raw = raw[which(ss > 0), ]
  
  pdfname = paste0(resDir, "/Data_qulity_assessment", version.analysis, ".pdf")
  pdf(pdfname, width = 20, height = 16)
  
  Check.RNAseq.Quality(read.count=raw, design.matrix = design[ ,c(4, 6)], 
                       lowlyExpressed.readCount.threshold=20)
  
  dev.off()
  
}


########################################################
########################################################
# Section : normalization and pairwise comparison with DESeq2 
# 
########################################################
########################################################
require(DESeq2)
require(ggplot2)
library(ggrepel)
require(gridExtra)
library(dplyr)
library(patchwork)
require(pheatmap)
library(org.Mm.eg.db)
library(enrichplot)
library(clusterProfiler)
library(stringr)

load(file = paste0(RdataDir, 'design_rawCounts_geneSymbols', version.analysis, '.Rdata'))

# add batch information
xx = openxlsx::read.xlsx(paste0(dataDir, '20221112_sequencingInformation_hannahs.xlsx'), sheet = 3)
design$batch = xx$batch[match(design$sampleID, xx$SampleID)]

# not consider sample 163095
sels = which(design$SampleID != '163095')

design.matrix = design[sels, ]
rm(design)

raw = counts[, sels]

dds <- DESeqDataSetFromMatrix(raw, DataFrame(design.matrix), design = ~ condition)

ss = rowSums(counts(dds))

hist(log10(ss), breaks = 100, main = 'log10(sum of reads for each gene)')

cutoff.gene = 10
cat(length(which(ss > cutoff.gene)), 'genes selected \n')

dds <- dds[ss > cutoff.gene, ]

# normalization and dimensionality reduction
dds = estimateSizeFactors(dds)
fpm = fpm(dds, robust = TRUE)

#save(dds, design.matrix,  file = paste0(RdataDir, 'TM3_dds_normalized.Rdata'))
#save(fpm, design, file = paste0(tfDir, '/RNAseq_fpm_fitered.cutoff.', cutoff.gene, '.Rdata'))
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('treat', 'time'), returnData = FALSE)
print(pca)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('treat', 'time'), returnData = TRUE, ntop = 3000))
pca2save$time = as.factor(pca2save$time)

ggplot(data=pca2save, aes(PC1, PC2, label = name, color= treat, shape = time))  + 
  geom_point(size=4) + 
  scale_shape_manual(values=1:nlevels(pca2save$time)) +
  geom_text(hjust = 0.5, nudge_y = 0.5, size=2.5)

ggsave(paste0(resDir, '/PCA_smartseq3_allSamples_filteredOneSAG.pdf'),  width=12, height = 10)

library(dplyr)

saveFigure = FALSE
if(saveFigure){
  pca2save$treat[which(pca2save$treat == '2iLIF')] = 'noRA'
  pca2save %>% filter(time %in% c('d0', 'd2', 'd2.18h', 'd3', 'd3.18h', 'd4', 'd4.18h')) %>%
    filter(treat %in% c('noRA', 'RA')) %>%
    ggplot(aes(PC1, PC2, label = name, color= time, shape = treat))  + 
    geom_point(size=4, aes(fill=time)) + 
    scale_shape_manual(values=1:nlevels(pca2save$time)) +
    geom_text(hjust = 0.5, nudge_y = 0.5, size=2.0) + 
    #geom_text_repel(size = 1.0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, size = 14), 
          axis.text.y = element_text(angle = 0, size = 14), 
          axis.title =  element_text(size = 14),
          legend.text = element_text(size=12),
          legend.title = element_text(size = 14)
    ) 
  
  ggsave(paste0("../results/plots_MondaySeminar", 
                '/PCA_smartseq3_RAresponse.pdf'),  width=8, height = 6) 
  
  
}


fpm = fpm[, order(colnames(fpm))]

#write.csv2(fpm, file = paste0(resDir, '/fpm_allSamples.csv'))
           
save(dds, fpm, design.matrix, file = paste0(RdataDir, 'design_dds_fpm_allSamples_filteredSAG.Rdata'))

########################################################
########################################################
# Section : Pairwise comparisons 
# Here we first identify significantly different genes in RA positive and negative genes
# RA positive vs negative
# In addition, the different response of Foxa2 positive and negative cells are also containing information of 
# which genes are expressed in positive or negative cells
########################################################
########################################################
load(file = paste0(RdataDir, 'design_dds_fpm_allSamples_filteredSAG.Rdata'))
cpm = fpm(dds)
colnames(cpm) = paste0(colnames(cpm), '.fpm')

fdr.cutoff = 0.05
log2fc.cutoff = 1

Compare.diffRAstimulationTime.Batch1 = FALSE
if(Compare.diffRAstimulationTime.Batch1){
  
  kk = which(design.matrix$batch == '1' & design.matrix$condition != '2iLIF_d0')
  tableDir = paste0(resDir, '/Compare.diffRAstimulationTime.Batch1')
  if(!dir.exists(tableDir)) dir.create(tableDir)
  
  dds1 = dds[, kk]
  design = design.matrix[kk, ]
  dds1$condition <- droplevels(dds1$condition)
  cpm1 = fpm(dds1)
  
  dds1 <- estimateDispersions(dds1, fitType = 'parametric')
  plotDispEsts(dds1, ymin = 10^-3, main = 'diff RA simulation time -- batch 1')
  
  dds1 <- nbinomWaldTest(dds1)
  resultsNames(dds1)  
  
  res = c()
  for(t in c(2, 3, 4))
  {
    # t = 2
    cat('-- day', t, '\n')
    ii = grep(paste0('_d', t), design$condition)
    cc = unique(design$condition[ii])
    cf = cc[grep(paste0('_d', t, '$'), cc)]
    cc = setdiff(cc, cf)
    c1 = cc[grep('^RA_', cc)]
    c2 = cc[grep('^noRA_', cc)]
    cc = c(c1, c2, setdiff(cc, c(c1, c2, cf)))
    cat(cf, ' and ', cc, '\n')
    
    for(n in 1:length(cc))
    {
      # n = 1
      res.ii = results(dds1, contrast = c('condition', cc[n], cf), test = 'Wald')
      res.ii <- lfcShrink(dds1, contrast = c('condition', cc[n], cf), type = 'normal')
      
      colnames(res.ii) = paste0(colnames(res.ii), "_", cc[n], "vs.", cf)
      
      if(is.null(res)){
        res = data.frame(res.ii, stringsAsFactors = FALSE)
      }else{
        res = data.frame(res, res.ii)
      }
      
      # add cpm signals
      xx = cpm1[, c(which(design$condition == cc[n]), which(design$condition == cf))]
      #xx = xx[, order(colnames(xx))]
      
      xx = data.frame(xx, res.ii[, c(2, 5, 6)])
      xx = xx[order(xx[,grep('pvalue_', colnames(xx))]), ]
      
      select = which(xx[,grep('^padj_', colnames(xx))] < fdr.cutoff & 
                       abs(xx[, grep('^log2FoldChange_', colnames(xx))]) > log2fc.cutoff)
      
      write.csv2(xx[select, ], file = paste0(tableDir, '/DESeq2_DEgenes_pairwiseComparison_', cc[n], '.vs.', cf, 
                 '.csv'), row.names = TRUE)
      
      if(n == 1){
        res.ii = results(dds1, contrast = c('condition', c1, c2), test = 'Wald')
        res.ii <- lfcShrink(dds1, contrast = c('condition', c1, c2), type = 'normal')
        
        colnames(res.ii) = paste0(colnames(res.ii), "_", c1, "vs.", c2)
        
        if(is.null(res)){
          res = data.frame(res.ii, stringsAsFactors = FALSE)
        }else{
          res = data.frame(res, res.ii)
        }
        
        # add cpm signals
        xx = cpm1[, c(which(design$condition == c1), which(design$condition == c2))]
        #xx = xx[, order(colnames(xx))]
        
        xx = data.frame(xx, res.ii[, c(2, 5, 6)])
        xx = xx[order(xx[,grep('pvalue_', colnames(xx))]), ]
        
        select = which(xx[,grep('^padj_', colnames(xx))] < fdr.cutoff & 
                         abs(xx[, grep('^log2FoldChange_', colnames(xx))]) > log2fc.cutoff)
        
        write.csv2(xx[select, ], file = paste0(tableDir, '/DESeq2_DEgenes_pairwiseComparison_', c1, '.vs.', c2, 
                                               '.csv'), row.names = TRUE)
      }
      
      if(n == 3){
        # RA_d2.18h vs SAG_d2.18h
        res.ii = results(dds1, contrast = c('condition', c1, cc[n]), test = 'Wald')
        res.ii <- lfcShrink(dds1, contrast = c('condition', c1, cc[n]), type = 'normal')
        
        colnames(res.ii) = paste0(colnames(res.ii), "_", c1, "vs.", cc[n])
        
        if(is.null(res)){
          res = data.frame(res.ii, stringsAsFactors = FALSE)
        }else{
          res = data.frame(res, res.ii)
        }
        
        # add cpm signals
        xx = cpm1[, c(which(design$condition == c1), which(design$condition == cc[n]))]
        #xx = xx[, order(colnames(xx))]
        
        xx = data.frame(xx, res.ii[, c(2, 5, 6)])
        xx = xx[order(xx[,grep('pvalue_', colnames(xx))]), ]
        
        select = which(xx[,grep('^padj_', colnames(xx))] < fdr.cutoff & 
                         abs(xx[, grep('^log2FoldChange_', colnames(xx))]) > log2fc.cutoff)
        
        write.csv2(xx[select, ], file = paste0(tableDir, '/DESeq2_DEgenes_pairwiseComparison_', c1, '.vs.', cc[n], 
                                               '.csv'), row.names = TRUE)
        
        # SAG_d2.18h vs. noRA_d2.18
        res.ii = results(dds1, contrast = c('condition', cc[n], c2), test = 'Wald')
        res.ii <- lfcShrink(dds1, contrast = c('condition', cc[n], c2), type = 'normal')
        
        colnames(res.ii) = paste0(colnames(res.ii), "_", cc[n], "vs.", c2)
        
        if(is.null(res)){
          res = data.frame(res.ii, stringsAsFactors = FALSE)
        }else{
          res = data.frame(res, res.ii)
        }
        
        # add cpm signals
        xx = cpm1[, c(which(design$condition == cc[n]), which(design$condition == c2))]
        #xx = xx[, order(colnames(xx))]
        
        xx = data.frame(xx, res.ii[, c(2, 5, 6)])
        xx = xx[order(xx[,grep('pvalue_', colnames(xx))]), ]
        
        select = which(xx[,grep('^padj_', colnames(xx))] < fdr.cutoff & 
                         abs(xx[, grep('^log2FoldChange_', colnames(xx))]) > log2fc.cutoff)
        
        write.csv2(xx[select, ], file = paste0(tableDir, '/DESeq2_DEgenes_pairwiseComparison_', cc[n], '.vs.', c2, 
                                               '.csv'), row.names = TRUE)
         
      }
      
    }
    
  }
  
}

saveRDS(res, file = paste0(RdataDir, 'res_12pairwise_comparisons_batch1.rds'))
  
##########################################
## perturbations in batch2
## 
##########################################
load(file = paste0(RdataDir, 'design_dds_fpm_allSamples_filteredSAG.Rdata'))
cpm = fpm(dds)
colnames(cpm) = paste0(colnames(cpm), '.fpm')

fdr.cutoff = 0.05
log2fc.cutoff = 1

Compare.diff.perturbationCompare.Batch2 = FALSE
if(Compare.diff.perturbationCompare.Batch2){
  
  kk = which(design.matrix$batch == '2')
  tableDir = paste0(resDir, '/Compare.diffRAperturbation.Batch2')
  if(!dir.exists(tableDir)) dir.create(tableDir)
  
  dds1 = dds[, kk]
  design = design.matrix[kk, ]
  dds1$condition <- droplevels(dds1$condition)
  cpm1 = fpm(dds1)
  
  dds1 <- estimateDispersions(dds1, fitType = 'parametric')
  plotDispEsts(dds1, ymin = 10^-3, main = 'diff RA perturbation -- batch 2')
  
  dds1 <- nbinomWaldTest(dds1)
  resultsNames(dds1)  
  
  res = c()
  
  cc = unique(design$condition)
  cf1 = "noRA_d2.18h" 
  cf2 = "RA_d2.18h"
  
  cc = setdiff(cc, c(cf1, cf2))
  cat('ref1 --', cf1, '\n')
  cat('ref2 --', cf2, '\n')
  cat('conditions to compare --', cc, '\n')
  
  rm(c1); rm(c2)
  
  # comapre the cf2 vs cf1 in the batch2
  c1 = cf2; c2 = cf1;
  cat('pairwise compare -- ', c1, 'vs.', c2, '\n')
  
  res.ii = results(dds1, contrast = c('condition', c1, c2), test = 'Wald')
  res.ii <- lfcShrink(dds1, contrast = c('condition', c1, c2), type = 'normal')
  
  colnames(res.ii) = paste0(colnames(res.ii), "_", c1, "vs.", c2)
  
  if(is.null(res)){
    res = data.frame(res.ii, stringsAsFactors = FALSE)
  }else{
    res = data.frame(res, res.ii)
  }
  
  # add cpm signals
  xx = cpm1[, c(which(design$condition == c1), which(design$condition == c2))]
  #xx = xx[, order(colnames(xx))]
  
  xx = data.frame(xx, res.ii[, c(2, 5, 6)])
  xx = xx[order(xx[,grep('pvalue_', colnames(xx))]), ]
  
  select = which(xx[,grep('^padj_', colnames(xx))] < fdr.cutoff & 
                   abs(xx[, grep('^log2FoldChange_', colnames(xx))]) > log2fc.cutoff)
  
  write.csv2(xx[select, ], file = paste0(tableDir, '/DESeq2_DEgenes_pairwiseComparison_', c1, '.vs.', c2, 
                                         '.csv'), row.names = TRUE)
  
  for(n in 1:length(cc))
  {
    # n = 1
    c1 = cc[n]
    c2 = cf1
    cat('pairwise compare -- ', c1, 'vs.', c2, '\n')
    
    res.ii = results(dds1, contrast = c('condition', c1, c2), test = 'Wald')
    res.ii <- lfcShrink(dds1, contrast = c('condition', c1, c2), type = 'normal')
    
    colnames(res.ii) = paste0(colnames(res.ii), "_", c1, "vs.", c2)
    
    if(is.null(res)){
      res = data.frame(res.ii, stringsAsFactors = FALSE)
    }else{
      res = data.frame(res, res.ii)
    }
    
    # add cpm signals
    xx = cpm1[, c(which(design$condition == c1), which(design$condition == c2))]
    #xx = xx[, order(colnames(xx))]
    
    xx = data.frame(xx, res.ii[, c(2, 5, 6)])
    xx = xx[order(xx[,grep('pvalue_', colnames(xx))]), ]
    
    select = which(xx[,grep('^padj_', colnames(xx))] < fdr.cutoff & 
                     abs(xx[, grep('^log2FoldChange_', colnames(xx))]) > log2fc.cutoff)
    
    write.csv2(xx[select, ], file = paste0(tableDir, '/DESeq2_DEgenes_pairwiseComparison_', c1, '.vs.', c2, 
                                           '.csv'), row.names = TRUE)
    
    c1 = cc[n]
    c2 = cf2
    cat('pairwise compare -- ', c1, 'vs.', c2, '\n')
    
    res.ii = results(dds1, contrast = c('condition', c1, c2), test = 'Wald')
    res.ii <- lfcShrink(dds1, contrast = c('condition', c1, c2), type = 'normal')
    
    colnames(res.ii) = paste0(colnames(res.ii), "_", c1, "vs.", c2)
    
    if(is.null(res)){
      res = data.frame(res.ii, stringsAsFactors = FALSE)
    }else{
      res = data.frame(res, res.ii)
    }
    
    # add cpm signals
    xx = cpm1[, c(which(design$condition == c1), which(design$condition == c2))]
    #xx = xx[, order(colnames(xx))]
    
    xx = data.frame(xx, res.ii[, c(2, 5, 6)])
    xx = xx[order(xx[,grep('pvalue_', colnames(xx))]), ]
    
    select = which(xx[,grep('^padj_', colnames(xx))] < fdr.cutoff & 
                     abs(xx[, grep('^log2FoldChange_', colnames(xx))]) > log2fc.cutoff)
    
    write.csv2(xx[select, ], file = paste0(tableDir, '/DESeq2_DEgenes_pairwiseComparison_', c1, '.vs.', c2, 
                                           '.csv'), row.names = TRUE)
    
    
  }
  
}

saveRDS(res, file = paste0(RdataDir, 'res_15pairwise_comparisons_batch2.rds'))


##########################################
# Intersect the DE genes with the RA targets identified from RAR chip-seq and chia-pet  
##########################################
RARtargetDir = '../results/RAR_targets/Rdata/'
peaks = readRDS(file = paste0(RARtargetDir, 
                              'RAR_chipseq.peak.assignment_promoters.genebody.downstream.chiapet.closestTSS.rds'))

ggs = unlist(peaks[, c(13:17)])
ggs = ggs[!is.na(ggs)]
ggs = unique(ggs)

#write.table(ggs, file = paste0(resDir, '/RAtargets_chipseq.peak.assignment_chiapet.closestTSS.txt'), 
#            sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)


xx = c()
for(n in 1:length(ggs))
{
  xx = c(xx, unlist(strsplit(as.character(ggs[n]), ';')))
}
ggs = unique(xx)

saveRDS(ggs, file = paste0(RdataDir, 'RAR_targets_chip_chiapet.rds'))

tableDir = '../results/RA_targets_L118404_smartseq3_20221117/Compare.diffRAperturbation.Batch2/'

file.list = list.files(path = paste0(tableDir, 'all'), pattern = '*.csv',full.names = TRUE)
for(n in 1:length(file.list))
{
  test = read.csv2(file = file.list[n], row.names = 1)
  ttest = test[which(!is.na(match(rownames(test), ggs))), ]
  cat(n, ' :  before ', nrow(test), ' -- after ', nrow(ttest), '\n')
  fname = basename(file.list[n])
  write.csv2(ttest, file = paste0(tableDir, 'RARtarget_intersection/', fname), 
             row.names = TRUE)
  
}


##########################################
# plot DE genes and RAR targets
##########################################
make.plot.summary = FALSE
if(make.plot.summary){
  library(tidyr)
  library(dplyr)
  
  cal_sample_means = function(cpm, conds = c("mUA", "mLA", "mHand") )
  {
    sample.means = c()
    for(n in 1:length(conds)) 
    {
      kk = grep(conds[n], colnames(cpm))
      if(length(kk)>1) {
        sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
      }else{
        sample.means = cbind(sample.means, cpm[, kk])
      }
      
    }
    colnames(sample.means) = conds
    return(sample.means)
  }
  
  load(file = paste0(RdataDir, 'design_dds_fpm_allSamples_filteredSAG.Rdata'))
  cpm = fpm(dds)
  colnames(cpm) = paste0(colnames(cpm), '.fpm')
  
  fdr.cutoff = 0.05
  log2fc.cutoff = 1
  
  res1 = readRDS(file = paste0(RdataDir, 'res_12pairwise_comparisons_batch1.rds'))
  res2 = readRDS(file = paste0(RdataDir, 'res_15pairwise_comparisons_batch2.rds'))
  res = data.frame(res1, res2)
  
  res = res[, grep('baseMean_|lfcSE|stat_', colnames(res), invert = TRUE)]
  
  res$pval.min = -log10(apply(res[, grep('pvalue_', colnames(res))], 1, min))
  res$log2fc.max = apply(res[, grep('log2FoldChange_', colnames(res))], 1, function(x) max(abs(x)))
  
  
  res = data.frame(cpm, res, stringsAsFactors = FALSE)  
  res = res[order(-res$log2fc.max), ]
  
  # write.csv2(res, file = paste0(resDir, '/fpm_stat.allpairwseCompares.csv'), row.names = TRUE)
  
  ##########################################
  # WNT, FGF, BMP genes in Foxa2 positive and negative cells
  ##########################################
  ggs = rownames(res)[which(res$pval.min > 3)]
  
  ggs = c(ggs, 'Shh')
  ggs = c(ggs, c('Lef1', 'Wnt3', 'Wnt3a', 'Wnt4', 'Wnt5b', 'Wnt6', 'Wnt7a', 'Wnt7b', 'Wnt8a', 'Dkk1', 'Dkk2', 
               'Dkk3', 'Tcf15', 'Tcf19', 
               "Wnt1", 'Sost', 'Sfrp5', 'Lypd6', 'Peg12', 'Notum', 'Draxin'))
  ggs= c(ggs, c('Spry4', 'Spry2', 'Etv4', 'Etv5', 'Fgf10', 'Fgf17','Fgf8', 'Fgf5', 'Fgf2','Fgf21', 
               'Fgf11','Fgf1', 'Fgf4', 'Fgfbp3',
               'Dusp1', 'Dusp10', 'Dusp27', 'Dusp4', 'Dusp5'))
  ggs = c(ggs, c('Id1', 'Id3', 'Smad6', 'Nog', 'Fst', 'Bambi', 'Bmp7', 'Bmp4', 'Bmp1', 'Bmp6', 
               'Bmpr2', 'Bmpr1b', 'Bmp2', 'Bmp3', 'Bmp5', 
               'Runx1', 'Runx2'))
  ggs = unique(ggs)
  
  Intersect.RAR.targets = FALSE
  if(Intersect.RAR.targets){
    targets = readRDS(file = paste0(RdataDir, 'RAR_targets_chip_chiapet.rds'))
    ggs = unique(intersect(ggs, targets))
    
  }
  
  cpm = cpm[!is.na(match(rownames(cpm), ggs)), ]
  rm(design)
  
  design.matrix$cc = paste0(design.matrix$condition, '_bc', design.matrix$batch)
  
  colnames(cpm) = design.matrix$cc
  
  cpm = cal_sample_means(cpm, conds = unique(design.matrix$cc))
  cpm = log2(cpm + 2^-6.5)
  
  # hist(rpkm.RA, breaks = 100, xlab = 'log2(RPKM)')
  # abline(v = c(0, 1), col = 'red', lwd = 2.0)
  # 
  # hist(rpkm.noRA, breaks = 100, xlab = 'log2(RPKM)')
  # abline(v = c(0, 1), col = 'red', lwd = 2.0)
  # 
  # hist(sorted, breaks = 100, xlab = 'log2(RPKM)')
  # abline(v = c(0, 1), col = 'red', lwd = 2.0)
  # 
  jj1 = grep('2iLIF|noRA_', colnames(cpm))
  jj1 = jj1[grep('_bc1', colnames(cpm)[jj1])]
  jj1 = jj1[order(colnames(cpm)[jj1])]
  
  kk1 = grep('_bc1', colnames(cpm))
  kk1 = kk1[grep('^RA_|^SAG', colnames(cpm)[kk1])]
  kk1 = kk1[order(colnames(cpm)[kk1])]
  
  jj2 = grep('_bc2', colnames(cpm))
  jj2 = jj2[order(colnames(cpm)[jj2])]
 
  
  pdfname = paste0(resDir, '/RAtargets_smartseq3_intersectingRAR_chipTargets_v1.5.pdf')
  pdf(pdfname,  width = 18, height = 8)
  par(cex = 1.0, las = 1, mgp = c(3,2,0), mar = c(6,6,2,0.2), tcl = -0.3)
  
  attach(mtcars)
  par(mfrow=c(1,2)) 
  
  for(n in 1:length(ggs))
  #for(n in 1:10)
  {
    # n = 1
    g = ggs[n]
    kk = which(rownames(cpm) == g)
    if(length(kk)){
      cat(n, '--', g, '\n')
      
      
      # batch 1 
      plot(c(0, 1), type = 'n', xlim = c(0, 120), ylim = range(c(cpm[kk, ], 0, 1)), main = g, 
           ylab = 'log2(CPM)', xlab = 'time')
      points(c(0, 48, 66, 72, 72+18, 96, 96+18), cpm[kk, jj1], col = 'darkblue', type = 'l', pch = 16, lwd = 2.0)
      points(c(0, 48, 66, 72, 72+18, 96, 96+18), cpm[kk, jj1], col = 'darkblue', type = 'p', pch = 1)
      
      #points(c(66, 66, 72+18, 96+18), cpm[kk, kk1], col = 'darkred', type = 'l', lwd = 1.0)
      points(c(48, 66), cpm[kk, c(jj1[2], kk1[1])], col = 'darkgreen', type = 'l', lwd = 2.5)
      points(c(48, 66), cpm[kk, c(jj1[2], kk1[1])], col = 'darkgreen', type = 'p', pch = 16, cex = 2)
      
      points(c(48, 66), cpm[kk, c(jj1[2], kk1[4])], col = 'magenta', type = 'l', lwd = 2.5)
      points(c(48, 66), cpm[kk, c(jj1[2], kk1[4])], col = 'magenta', type = 'p', pch = 16, cex = 2)
      
      points(c(72, 72+18), cpm[kk, c(jj1[4], kk1[2])], col = 'darkgreen', type = 'l', lwd = 2.5)
      points(c(72, 72+18), cpm[kk, c(jj1[4], kk1[2])], col = 'darkgreen', type = 'p', pch = 16, cex = 2)
      
      points(c(96, 96+18), cpm[kk, c(jj1[6], kk1[3])], col = 'darkgreen', type = 'l', lwd = 2.5)
      points(c(96, 96+18), cpm[kk, c(jj1[6], kk1[3])], col = 'darkgreen', type = 'p', pch = 16, cex = 2)
      # 
      # points(48, sorted[kk, 1], col = 'darkblue', type = 'p', cex = 2.0, pch = 21, bg = 'magenta')
      # points(48, sorted[kk, 2], col = 'darkblue', type = 'p', cex = 2.0, pch = 21, bg = 'darkgreen')
      # points(48, sorted[kk, 3], col = 'deepskyblue', type = 'p', cex = 2.0, pch = 18, bg = 'darkgreen')
      
      abline(h = c(0, 1), col = 'darkgray', lwd = 2.0)
      #abline(v = c(6, 32, 54), col = 'cornflowerblue', lwd = 2.0, lty=3)
      legend('topleft', legend = c('noRA', 'RA', 'SAG'), bty = 'n', 
             col = c('darkblue', 'darkgreen', 'magenta'),  lwd =2.0, 
             pch = c(1, 16, 16), lty = c(1, 1, 1))
      
      # batch 2
      names =  colnames(cpm)[jj2]
      names = sapply(names, function(x){unlist(strsplit(as.character(x), '_'))[1]})
      barplot(cpm[kk, jj2], names.arg = names, 
              horiz = FALSE, las = 2, col = c('darkblue', 'darkgreen', c(1:7)), ylab = 'log2(CPM)')
      
            
    }else{
      cat(g, 'Not Found \n')
      
    }
  }
  
  dev.off()
  
}


##########################################
# highlight signaling pathways
##########################################
Explore.perturbation.response.for.signaling.pathway.genes = FALSE
if(Explore.perturbation.response.for.signaling.pathway.genes){
  
  res = readRDS(file = paste0(RdataDir, '/TM3_res_pairwiseComparisons_perturbation.vs.RA_positive.negative.pooled.rds'))
  ggs = readRDS(file = paste0(RdataDir, '/TM3_examplesGenes_withGOterm.rds'))
  
  Compare.positive.vs.pooled.samples = FALSE
  
  cc = unique(design.matrix$condition)
  cc = cc[which(cc != "N2B27" & cc != 'RA')]
  
  res = res[!is.na(match(rownames(res), ggs$gene)), grep('pvalue_|log2FoldChange_', colnames(res))]
  
  library(ggrepel)
  library(ggplot2)
  
  examples = c('Bmp4', 'Bmp1', 'Bmp7', 'Bmp6', 'Nog', 'Dkk1', 'Fgfbp3', 'Lef1', 'Tcf15', 'Tcf19', 'Wnt3', 'Wnt3a', 'Wnt4', 
               'Wnt5b', 'Wnt6', 'Wnt7a', 'Wnt7b', 'Wnt8a', 'Sfrp5', 'Lypd6', 'Spry4', 'Etv5', 'Etv4', 'Fgf10', 'Fgf17', 
               'Fgf8', 'Fgf5', 'Fgf2', 'Fgf21', 'Fgf11', 'Fgf1', 'Fgf4')
  
  
  
  pdfname = paste0(resDir, '/TM3_comparing_response_signalingPathways_positve.negative_v2.pdf')
  pdf(pdfname,  width = 16, height = 12)
  #par(cex = 1.0, las = 1, mgp = c(3,2,0), mar = c(6,6,2,0.2), tcl = -0.3)
  
  for(n in 1:length(cc))
  {
    # n = 5
    cat(n, ' -- ', cc[n], '\n')
    jj = grep(cc[n], colnames(res))
    yy = res[, jj]
    colnames(yy) = gsub(paste0(cc[n]), '', colnames(yy))
    colnames(yy) = gsub('_vs_RA.', '', colnames(yy))
    colnames(yy) = gsub('.pos|.neg|.pooled', '', colnames(yy))
    colnames(yy) = gsub('log2FoldChange', 'lfc', colnames(yy))
    
    yy[, grep('pvalue_', colnames(yy))] = -log10(yy[, grep('pvalue_', colnames(yy))])
    yy = data.frame(yy, gene = rownames(yy))
    
    p1 = ggplot(data = yy,  aes(y = pvalue_pos, x = lfc_pos,  label = gene)) + 
      geom_point(size = 1) + 
      labs(title = paste0(cc[n], " - positive  "), x = '', y = '-log10(pval)') + 
      theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12)) + 
      geom_hline(yintercept = 2, colour = "red") +
      #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
      geom_text_repel(data=subset(yy, pvalue_pos > 2), size = 4)
    
    p2 = ggplot(data = yy,  aes(y = pvalue_neg, x = lfc_neg,  label = gene)) + 
      geom_point(size = 1) + 
      labs(title = paste0(cc[n], " - neg  "), x = '', y = '-log10(pval)') + 
      theme(axis.text.x = element_text(size = 12), 
            axis.text.y = element_text(size = 12)) + 
      geom_hline(yintercept = 2, colour = "red") +
      #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
      geom_text_repel(data=subset(yy, pvalue_neg > 2), size = 4)
    
    if(Compare.positive.vs.pooled.samples){
      p3 = ggplot(data = yy,  aes(y = pvalue_pooled, x = lfc_pooled,  label = gene)) + 
        geom_point(size = 1) + 
        labs(title = paste0(cc[n], " - pooled  "), x = '', y = '-log10(pval)') + 
        theme(axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12)) + 
        geom_hline(yintercept = 2, colour = "red") +
        #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
        geom_text_repel(data=subset(yy, pvalue_pooled > 2), size = 4)
    }
    
    examples.sel = unique(c(examples, rownames(yy)[which(yy$pvalue_pos > 2 | yy$pvalue_neg > 2)]))
    
    if(Compare.positive.vs.pooled.samples){
      p4 = ggplot(data = yy, aes(x = pvalue_pooled, y = pvalue_neg , label = gene)) +
        geom_point(size = 1) + 
        xlim(0, 6) + ylim(0, 6) + geom_hline(yintercept = 1.3, colour = "blue") +
        geom_vline(xintercept = 1.3, colour = "blue") + 
        labs(title = paste0(cc[n], " - neg vs pooled : -log10(pval) "), x = 'pooled', y = 'negative') +
        geom_abline(slope = 1, intercept = 0, colour = 'red') + 
        geom_label_repel(data=  as.tibble(yy) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel),
                         size = 3)
      #examples.sel = unique(c(examples, rownames(yy)[which(yy$pvalue_pos > 2 | yy$pvalue_neg > 2)]))
      p5 = ggplot(data = yy, aes(x = pvalue_pooled, y = pvalue_pos, label = gene)) +
        geom_point(size = 1) + 
        xlim(0, 6) + ylim(0, 6) + geom_hline(yintercept = 1.3, colour = "blue") +
        geom_vline(xintercept = 1.3, colour = "blue") + 
        labs(title = paste0(cc[n], " - pos vs pooled : -log10(pval) "), x = 'pooled', y = 'positive') + 
        geom_abline(slope = 1, intercept = 0, colour = 'red') + 
        geom_label_repel(data=  as.tibble(yy) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel),
                         size = 3)
    }
    
    p6 = ggplot(data = yy, aes(x = pvalue_neg, y = pvalue_pos, label = gene)) +
      geom_point(size = 1) + 
      xlim(0, 6) + ylim(0, 6) + geom_hline(yintercept = 2, colour = "blue") +
      geom_vline(xintercept = 2, colour = "blue") + 
      labs(title = paste0(cc[n], " - pos vs neg : -log10(pval) "), x = 'negative', y = 'positive') + 
      geom_abline(slope = 1, intercept = 0, colour = 'red') + 
      geom_label_repel(data=  as.tibble(yy) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel),
                       size = 3)
    
    plot(p1)
    plot(p2)
    if(Compare.positive.vs.pooled.samples){
      plot(p3)
      plot(p4)
      plot(p5)
    }
    plot(p6)
    
    #grid.arrange(p1, p2, p3, p4, p5, nrow = 2, ncol = 3)
  }
  
  dev.off()
  
}


##########################################
# rename the bigwig file of smarstq3 data
##########################################
Rename_bigwig_files = FALSE
if(Rename_bigwig_files){
  load(file = paste0(RdataDir, 'design_rawCounts_geneSymbols', version.analysis, '.Rdata'))
  
  design$sampleName = paste0(design$time, '_', design$treat, '_', design$SampleID)
  
  bwDir = '../L118404_smartseq3/bigwigs_deeptools'
  cwd = getwd()
  
  setwd(bwDir)
  
  bwFiles = list.files('.', pattern = '*.bw')
  
  for(n in 1:nrow(design))
  {
    # n = 1
    kk = grep(design$SampleID[n], bwFiles)
    if(length(kk) != 1) {
      cat('no file for ', design$SampleID[n], '\n')
    }else{
      system(paste0('mv ', bwFiles[kk], ' ', design$sampleName[n], '.bw'))
    }
    
  }
  
  setwd(cwd)
  
  
}

Rename_bigwig_ElenaFull.length = FALSE
if(Rename_bigwig_ElenaFull.length){
  
  cwd = getwd()
  
  design = readRDS(file = '../RNAseq_Elenea_old/sampleInfo_design.rds')
  design$SampleID = sapply(design$Linking_id, function(x) {unlist(strsplit(as.character(x), '-'))[2]})
  
  design$sampleName = paste0(design$condition, '_', design$SampleID)
  design = design[grep('AF|GFPp', design$Linking_id, invert = TRUE), ]
  
  bwDir = '../RNAseq_Elenea_old/bigwigs'
  
  setwd(bwDir)
  
  bwFiles = list.files('.', pattern = '*markDups.bam.all.bw')
  
  for(n in 1:nrow(design))
  {
    # n = 1
    kk = grep(design$SampleID[n], bwFiles)
    if(length(kk) != 1) {
      cat('no file for ', design$SampleID[n], '\n')
    }else{
      cat(paste0('mv ', bwFiles[kk], ' ', design$sampleName[n], '.bw\n'))
      system(paste0('mv ', bwFiles[kk], ' ', design$sampleName[n], '.bw'))
    }
    
  }
  
  setwd(cwd)
  
}

