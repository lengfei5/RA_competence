##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: analyze the hNT quantseq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Jun 14 13:06:24 2023
##########################################################################
##########################################################################
rm(list = ls())

RNA.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

# setup for data import and sequencing QCs
version.analysis = 'R15282_hNT_quantseq_20230614'

resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R15282_hNT_quantseq/'

Collect.QCs.stat = TRUE
Counts.to.Use = "UMI"


sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)

tfs = toupper(unique(c(tfs,sps)))

source('functions_utility.R')

########################################################
########################################################
# Section I : data processing and sequencing quality controls
# 
########################################################
########################################################

##########################################
# prepare design matrix, count table and QC table
##########################################
design = read.csv(paste0(dataDir, 'sampleInfos_hNT_quantseq.csv'))
design = design[, c(8, 3, 4, 5, 6, 7)]

colnames(design)[c(1,2)] = c('sampleID', 'timepoint')
design$timepoint = gsub('[+]','.', design$timepoint)

##################################################
## Import UMI count table
##################################################
source(RNA.functions)
source(RNA.QC.functions)

Import.onlyUMI = TRUE

#design = readRDS(file = paste0(RdataDir, 'desing_sampleInfo.rds'))
colnames(design)[1] = 'SampleID'
#colnames(design)[4] = 'SampleID'

xlist = list.files(path=paste0(dataDir, 'htseq_counts_BAMs_umi'),
                   pattern = "*.txt$", full.names = TRUE) ## list of data set to merge

all = cat.countTable(xlist, countsfrom = 'htseq')

counts = process.countTable(all=all, design = design[, c(1, 2:5)], merge.technicalRep.sameID = FALSE)
rownames(counts) = counts$gene
counts = as.matrix(counts[, -1])

save(design, counts, file=paste0(RdataDir, 'design_rawCounts', version.analysis, '.Rdata'))


##########################################
# QCs of replicates and conditions
##########################################
load(file=paste0(RdataDir, 'design_rawCounts', version.analysis, '.Rdata'))
design$condition = paste0(design[, c(2)], "_", design[,3], "_", design[,4], '_', design[,5])

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
  
  Check.RNAseq.Quality(read.count=raw, design.matrix = design[ , c(1, 7)], 
                       lowlyExpressed.readCount.threshold=20)
  
  dev.off()
  
}

########################################################
########################################################
# Section II : normalization and pairwise comparison with DESeq2 
# 
########################################################
########################################################
require(DESeq2)
require(ggplot2)
library(ggrepel)
require(gridExtra)
library(dplyr)
library(ggrepel)
library(patchwork)
require(pheatmap)
library(org.Mm.eg.db)
library(enrichplot)
library(clusterProfiler)
library(stringr)


load(file = paste0(RdataDir, 'design_rawCounts', version.analysis, '.Rdata'))
design$condition = paste0(design[, c(2)], "_", design[,3], "_", design[,4], '_', design[,5])

dds <- DESeqDataSetFromMatrix(counts, DataFrame(design), design = ~ condition)

ss = rowSums(counts(dds))

hist(log10(ss), breaks = 100, main = 'log10(sum of reads for each gene)')

cutoff.gene = 20
cat(length(which(ss > cutoff.gene)), 'genes selected \n')

dds <- dds[ss > cutoff.gene, ]

# normalization and dimensionality reduction
dds = estimateSizeFactors(dds)
fpm = fpm(dds, robust = TRUE)

#save(dds, design.matrix,  file = paste0(RdataDir, 'TM3_dds_normalized.Rdata'))
#save(fpm, design, file = paste0(tfDir, '/RNAseq_fpm_fitered.cutoff.', cutoff.gene, '.Rdata'))
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('timepoint', 'Wash', 'Culture', 'RA'), returnData = FALSE)
print(pca)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('timepoint', 'Wash', 'Culture', 'RA'), 
                                 returnData = TRUE, ntop = 3000))
pca2save$time = as.factor(pca2save$timepoint)
pca2save$cc = paste0(pca2save$Wash, '_', pca2save$Culture, '_', pca2save$RA)

ggplot(data=pca2save, aes(PC1, PC2, label = name, color= cc, shape = time))  + 
  geom_point(size=3) + 
  scale_shape_manual(values=1:nlevels(pca2save$time)) +
  #geom_text(hjust = 0.5, nudge_y = 0.5, size=3.0) + 
  geom_text_repel(data= pca2save, size = 3) 

ggsave(paste0(resDir, '/PCA_quantseq_allSamples.pdf'),  width=20, height = 14)

save(dds, design, file = paste0(RdataDir, 'design_dds_', version.analysis, '.Rdata'))


########################################################
########################################################
# Section III: pairwise comparison 
#  
########################################################
########################################################
load(file = paste0(RdataDir, 'design_dds_', version.analysis, '.Rdata'))
design$treatment = paste0(design$Wash, '_', design$Culture)

##########################################
# the control sample: mTesr_LDNSB 
##########################################
treatment = 'mTesr_LDNSB'
outDir = paste0(resDir, '/out_', treatment, '/')
if(!dir.exists(outDir)) dir.create(outDir)

jj = which(design$treatment == treatment)

dds1 = dds[, jj]

dds1$condition <- droplevels(dds1$condition)
dds1$condition = factor(gsub(' ','.', dds1$condition))
dds1 <- estimateDispersions(dds1, fitType = 'parametric')
plotDispEsts(dds1, ymin = 10^-3)


cpm = fpm(dds1, robust = TRUE)

vsd <- varianceStabilizingTransformation(dds1, blind = FALSE)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('timepoint', 'Wash', 'Culture', 'RA'), 
                                 returnData = TRUE, ntop = 1000))
pca2save$time = as.factor(pca2save$timepoint)
pca2save$treatment = paste0(pca2save$Wash, '_', pca2save$Culture, '_', pca2save$RA)

ggplot(data=pca2save, aes(PC1, PC2, label = name, color= treatment, shape = time))  + 
  geom_point(size=3) + 
  scale_shape_manual(values=1:nlevels(pca2save$time)) +
  #geom_text(hjust = 0.5, nudge_y = 0.5, size=3.0) + 
  geom_text_repel(data= pca2save, size = 3) 

ggsave(paste0(outDir, '/PCA_quantseq_samples_mTesr_LDNSB.pdf'),  width=20, height = 14)

#dev.off()

dds1 <- nbinomWaldTest(dds1)
resultsNames(dds1)  

compares = c()
cc = unique(design$timepoint[jj])

for(n in seq_len(length(cc)))
{
  if(n>1){
    compares = c(compares, list(c('condition', paste0(cc[n], '_mTesr_LDNSB_no'), 
                                  paste0(cc[n-1], '_mTesr_LDNSB_no'))))
  }
}
compares = c(compares, 
             list(c('condition', 'd2.18h_mTesr_LDNSB_d2.250nM', 'd2.18h_mTesr_LDNSB_no')), 
             list(c('condition', 'd4.18h_mTesr_LDNSB_d4.250nM', 'd4_mTesr_LDNSB_no')))




ggs = c()
for(n in seq_len(length(compares)))
{
  # n = 6
  res = results(dds1, contrast=compares[[n]], alpha = 0.05)
  res<- lfcShrink(dds1, contrast=compares[[n]], type="normal")
  res = as.data.frame(res)
  
  # save the significant genes
  gg.signif = rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)] 
  cat(length(gg.signif), ' significant genes for comparison: ', compares[[n]], '  \n')
  ggs = c(ggs, gg.signif)
  
  yy = res
  colnames(yy) = gsub('log2FoldChange', 'lfc', colnames(yy))
  yy$pvalue = -log10(yy$padj)
  yy = data.frame(yy, gene = rownames(yy))
  yy = yy[which(!is.na(yy$lfc & !is.na(yy$pvalue))), ]
  
  ggplot(data = yy,  aes(y = pvalue, x = lfc,  label = gene)) + 
    geom_point(size = 1) + 
    labs(title = paste0(compares[[n]][c(2:3)], collapse = '_vs._'), x = 'log2FC', y = '-log10(pval)') + 
    theme_classic()+
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12)) + 
    #geom_hline(yintercept = 2, colour = "red") +
    #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
    geom_text_repel(data=subset(yy, pvalue > 3|gene == 'Nog'), size = 4) 
  
  ggsave(filename = paste0(outDir, 'volcanoPlot_', paste0(compares[[n]][c(2:3)], collapse = '_vs._'),
                           '.pdf'), width=12, height = 8)
  
  write.csv(yy, file = paste0(outDir, 'DE_table_', paste0(compares[[n]][c(2:3)], collapse = '_vs._'), '.csv'))
  
}

ggs = unique(ggs)
saveRDS(ggs, file = paste0(outDir, 'DE_genes_fdr0.05_logfc.1.rds'))

ggs = readRDS(paste0(outDir,'DE_genes_fdr0.05_logfc.1.rds'))

cat(length(ggs), ' DE genes found in ', treatment, '\n')

library(grid)

paletteLength <- 50

cpm = fpm(dds1, robust = TRUE)
cpm = cpm[match(ggs, rownames(cpm)), ]

cc = unique(dds1$condition)

mat = matrix(NA, nrow = nrow(cpm), ncol = length(cc))
colnames(mat) = cc
rownames(mat) = rownames(cpm)

for(n in seq_len(ncol(mat)))
{
  mat[,n] = apply(cpm[, which(dds1$condition == colnames(mat)[n])], 1, median)
}

#mat = mat[, c(2,4,6,1,3,5)]
mat = mat[!is.na(match(rownames(mat), tfs)), ]

mat1 <- t(scale(t(mat)))

myColor <- viridis::viridis(paletteLength)
#myColor1 <- colorRampPalette(ArchRPalettes$coolwarm)(paletteLength)
#myColor2 <- colorRampPalette(c('lightgray', 'red'))(paletteLength)

myBreaks <- c(seq(min(mat1), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat1)/paletteLength, max(mat1), length.out=floor(paletteLength/2)))

# set annotation
anno_col <- data.frame(cbind(time = c('d0.8h', 'd1', 'd2', 'd2.18h', 'd2.18h', 'd4', 'd4.18h', 'd4.18h'), 
                             treatment = c('noRA', 'noRA', 'noRA', 'noRA', 'RA', 'noRA', 'noRA', 'RA')))
rownames(anno_col) <- colnames(mat1)

#col<- colorRampPalette(c("blue4", "white", "darkred"))(paletteLength)
#col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdGy")))(paletteLength)
#col = colorRampPalette(rev(brewer.pal(n = 7, name ="PuOr")))(16)
#col = palette(hcl.colors(8, "Viridis"))
#col = colorRampPalette(c("navy", "white", "red3"))(paletteLength)

library(pheatmap)
pheatmap::pheatmap(
  mat1,
  border_color = 'gray60',
  color = myColor,
  breaks = myBreaks,
  annotation_col = anno_col,
  angle_col = 45,
  cluster_cols = FALSE,
  fontsize_row = 4,
  treeheight_row = 30,
  filename = paste0(outDir, 'heatmap_', treatment, '_with.geneNames.pdf'), 
  width = 8, height = 30
)

pheatmap::pheatmap(
  mat1,
  border_color = 'gray60',
  color = myColor,
  breaks = myBreaks,
  annotation_col = anno_col,
  angle_col = 45,
  cluster_cols = FALSE,
  fontsize_row = 4,
  show_rownames = FALSE,
  treeheight_row = 30,
  filename = paste0(outDir, 'heatmap_sortedCells_conditions_without.geneNames.pdf'), 
  width = 5, height = 10
)


##########################################
# compare d2 samples
##########################################
treatment = 'compare_d2_ABCDE'
outDir = paste0(resDir, '/out_', treatment, '/')
if(!dir.exists(outDir)) dir.create(outDir)

jj = which(design$timepoint == 'd2')
dds1 = dds[, jj]

dds1$condition <- droplevels(dds1$condition)
dds1$condition = factor(gsub(' ','.', dds1$condition))
dds1 <- estimateDispersions(dds1, fitType = 'parametric')
plotDispEsts(dds1, ymin = 10^-3)

vsd <- varianceStabilizingTransformation(dds1, blind = FALSE)
pca2save = as.data.frame(plotPCA(vsd, intgroup = c('timepoint', 'Wash', 'Culture', 'RA'), 
                                 returnData = TRUE, ntop = 1000))

pca2save$time = as.factor(pca2save$timepoint)
pca2save$treatment = paste0(pca2save$Wash, '_', pca2save$Culture, '_', pca2save$RA)

ggplot(data=pca2save, aes(PC1, PC2, label = name, color= treatment, shape = time))  + 
  geom_point(size=3) + 
  scale_shape_manual(values=1:nlevels(pca2save$time)) +
  #geom_text(hjust = 0.5, nudge_y = 0.5, size=3.0) + 
  geom_text_repel(data= pca2save, size = 3) 

ggsave(paste0(outDir, '/PCA_quantseq_samples_', treatment, '.pdf'),  width=16, height = 10)

#dev.off()
dds1 <- nbinomWaldTest(dds1)
resultsNames(dds1)  

levels(dds1$condition)

compares = c()
compares = c(compares, 
             list(c('condition', 'd2_mTesr_SB_no', 'd2_mTesr_Adv_no')),  # C vs E
             list(c('condition', 'd2_mTesr_LDN_no', 'd2_mTesr_Adv_no')), # D vs E
             list(c('condition', 'd2_mTesr_SB_no', 'd2_mTesr_LDNSB_no')), # C vs B
             list(c('condition', 'd2_mTesr_LDN_no', 'd2_mTesr_LDNSB_no')), # D vs B
             list(c('condition', 'd2_mTesr_SB_no', 'd2_mTesr_LDNSB_no')), # C vs D
             list(c('condition', 'd2_mTesr_SB_no', 'd2_mTesr_LDN_no')), # C vs B
             list(c('condition', 'd2_mTesr_Adv_no', 'd2_mTesr_LDNSB_no')), # E vs D
             list(c('condition', 'd2_mTesr_Adv_no', 'd2_mTesr_LDN_no')), # E vs B
             list(c('condition', 'd2_mTesr_LDN_no', 'd2_mTesr_LDNSB_no'))
             )

ggs = c()
for(n in seq_len(length(compares)))
{
  # n = 6
  res = results(dds1, contrast=compares[[n]], alpha = 0.05)
  res<- lfcShrink(dds1, contrast=compares[[n]], type="normal")
  res = as.data.frame(res)
  
  # save the significant genes
  gg.signif = rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)] 
  cat(length(gg.signif), ' significant genes for comparison: ', compares[[n]], '  \n')
  ggs = c(ggs, gg.signif)
  
  yy = res
  colnames(yy) = gsub('log2FoldChange', 'lfc', colnames(yy))
  yy$pvalue = -log10(yy$padj)
  yy = data.frame(yy, gene = rownames(yy))
  yy = yy[which(!is.na(yy$lfc & !is.na(yy$pvalue))), ]
  
  ggplot(data = yy,  aes(y = pvalue, x = lfc,  label = gene)) + 
    geom_point(size = 1) + 
    labs(title = paste0(compares[[n]][c(2:3)], collapse = '_vs._'), x = 'log2FC', y = '-log10(pval)') + 
    theme_classic()+
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12)) + 
    #geom_hline(yintercept = 2, colour = "red") +
    #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
    geom_text_repel(data=subset(yy, pvalue > 3|gene == 'Nog'), size = 4) 
  
  ggsave(filename = paste0(outDir, 'volcanoPlot_', paste0(compares[[n]][c(2:3)], collapse = '_vs._'),
                           '.pdf'), width=12, height = 8)
  
  write.csv(yy, file = paste0(outDir, 'DE_table_', paste0(compares[[n]][c(2:3)], collapse = '_vs._'), '.csv'))
  
}

ggs = unique(ggs)
saveRDS(ggs, file = paste0(outDir, 'DE_genes_fdr0.05_logfc.1.rds'))

ggs = readRDS(file = paste0(outDir, 'DE_genes_fdr0.05_logfc.1.rds'))
cat(length(ggs), ' DE genes found in ', treatment, '\n')

library(grid)

paletteLength <- 50

cpm = fpm(dds1, robust = TRUE)
cpm = cpm[match(ggs, rownames(cpm)), ]

cc = unique(dds1$condition)

mat = matrix(NA, nrow = nrow(cpm), ncol = length(cc))
colnames(mat) = cc
rownames(mat) = rownames(cpm)

for(n in seq_len(ncol(mat)))
{
  mat[,n] = apply(cpm[, which(dds1$condition == colnames(mat)[n])], 1, median)
}

mat = mat[, c(2, 4, 3, 5, 1)]

# mat = mat[!is.na(match(rownames(mat), tfs)), ]
mat1 <- t(scale(t(mat)))

myColor <- viridis::viridis(paletteLength)
#myColor1 <- colorRampPalette(ArchRPalettes$coolwarm)(paletteLength)
#myColor2 <- colorRampPalette(c('lightgray', 'red'))(paletteLength)

myBreaks <- c(seq(min(mat1), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat1)/paletteLength, max(mat1), length.out=floor(paletteLength/2)))

# set annotation
anno_col <- data.frame(cbind(time = c('mTesr_LDNSB', 'mTesr_LDN', 'mTesr_Adv', 
                                      'mTesr_SB', 'Adv_LDNSB'), 
                             treatment = rep('noRA', time = ncol(mat1))))
rownames(anno_col) <- colnames(mat1)

#col<- colorRampPalette(c("blue4", "white", "darkred"))(paletteLength)
#col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdGy")))(paletteLength)
#col = colorRampPalette(rev(brewer.pal(n = 7, name ="PuOr")))(16)
#col = palette(hcl.colors(8, "Viridis"))
#col = colorRampPalette(c("navy", "white", "red3"))(paletteLength)

library(pheatmap)
pheatmap::pheatmap(
  mat1,
  border_color = 'gray60',
  color = myColor,
  breaks = myBreaks,
  annotation_col = anno_col,
  angle_col = 45,
  cluster_cols = FALSE,
  fontsize_row = 4,
  treeheight_row = 30,
  filename = paste0(outDir, 'heatmap_', treatment, '_with.geneNames.pdf'), 
  width = 8, height = 30
)

pheatmap::pheatmap(
  mat1,
  border_color = 'gray60',
  color = myColor,
  breaks = myBreaks,
  annotation_col = anno_col,
  angle_col = 45,
  cluster_cols = FALSE,
  fontsize_row = 4,
  show_rownames = FALSE,
  treeheight_row = 30,
  filename = paste0(outDir, 'heatmap_sortedCells_conditions_without.geneNames.pdf'), 
  width = 5, height = 10
)


##########################################
# Compare mTesr-LDNSB with mTesr_Adv and Adv_LDNSB
##########################################
treatment = 'mTesr_LDNSB.vs.mTesr_Adv.vs.Adv_LDNSB'
outDir = paste0(resDir, '/out_', treatment, '/')
if(!dir.exists(outDir)) dir.create(outDir)

jj = which(!is.na(match(design$treatment, c('mTesr_LDNSB', 'mTesr_Adv', 'Adv_LDNSB'))))
jj = jj[which(design$timepoint[jj] != 'd4' & design$timepoint[jj] != 'd4.18h')]

dds1 = dds[, jj]

dds1$condition <- droplevels(dds1$condition)
dds1$condition = factor(gsub(' ','.', dds1$condition))
dds1 <- estimateDispersions(dds1, fitType = 'parametric')
plotDispEsts(dds1, ymin = 10^-3)

vsd <- varianceStabilizingTransformation(dds1, blind = FALSE)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('timepoint', 'Wash', 'Culture', 'RA'), 
                                 returnData = TRUE, ntop = 1000))
pca2save$time = as.factor(pca2save$timepoint)
pca2save$treatment = paste0(pca2save$Wash, '_', pca2save$Culture)

ggplot(data=pca2save, aes(PC1, PC2, label = name, color= treatment, shape = time))  + 
  geom_point(size=3) + 
  scale_shape_manual(values=1:nlevels(pca2save$time)) +
  #geom_text(hjust = 0.5, nudge_y = 0.5, size=3.0) + 
  geom_text_repel(data= pca2save, size = 3) 

ggsave(paste0(outDir, '/PCA_quantseq_samples_mTesr_LDNSB.pdf'),  width=20, height = 14)


dds1 <- nbinomWaldTest(dds1)
resultsNames(dds1)  

levels(dds1$condition)

cc = unique(dds1$timepoint)
dds1$treatment = paste0(dds1$Wash, "_", dds1$Culture)

compares = c()
for(n in seq_len(length(cc)))
{
  compares = c(compares, 
               list(c('condition', paste0(cc[n], '_mTesr_Adv_no'), paste0(cc[n], '_mTesr_LDNSB_no'))),
               list(c('condition', paste0(cc[n], '_Adv_LDNSB_no'), paste0(cc[n], '_mTesr_LDNSB_no')))
               )
}

compares = c(compares, 
             list(c('condition', 'd2.18h_mTesr_Adv_d2.250nM', 'd2.18h_mTesr_LDNSB_d2.250nM')), 
             list(c('condition', 'd2.18h_Adv_LDNSB_d2.250nM', 'd2.18h_mTesr_LDNSB_d2.250nM')))



ggs = c()
for(n in seq_len(length(compares)))
{
  # n = 6
  res = results(dds1, contrast=compares[[n]], alpha = 0.05)
  res<- lfcShrink(dds1, contrast=compares[[n]], type="normal")
  res = as.data.frame(res)
  
  # save the significant genes
  gg.signif = rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)] 
  cat(length(gg.signif), ' significant genes for comparison: ', compares[[n]], '  \n')
  ggs = c(ggs, gg.signif)
  
  yy = res
  colnames(yy) = gsub('log2FoldChange', 'lfc', colnames(yy))
  yy$pvalue = -log10(yy$padj)
  yy = data.frame(yy, gene = rownames(yy))
  yy = yy[which(!is.na(yy$lfc & !is.na(yy$pvalue))), ]
  
  ggplot(data = yy,  aes(y = pvalue, x = lfc,  label = gene)) + 
    geom_point(size = 1) + 
    labs(title = paste0(compares[[n]][c(2:3)], collapse = '_vs._'), x = 'log2FC', y = '-log10(pval)') + 
    theme_classic()+
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12)) + 
    #geom_hline(yintercept = 2, colour = "red") +
    #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
    geom_text_repel(data=subset(yy, pvalue > 3|gene == 'Nog'), size = 4) 
  
  ggsave(filename = paste0(outDir, 'volcanoPlot_', paste0(compares[[n]][c(2:3)], collapse = '_vs._'),
                           '.pdf'), width=12, height = 8)
  
  write.csv(yy, file = paste0(outDir, 'DE_table_', paste0(compares[[n]][c(2:3)], collapse = '_vs._'), '.csv'))
  
}

ggs = unique(ggs)
saveRDS(ggs, file = paste0(outDir, 'DE_genes_fdr0.05_logfc.1.rds'))

ggs = readRDS(file = paste0(outDir, 'DE_genes_fdr0.05_logfc.1.rds'))
cat(length(ggs), ' DE genes found in ', treatment, '\n')

library(grid)

paletteLength <- 50

cpm = fpm(dds1, robust = TRUE)
cpm = cpm[match(ggs, rownames(cpm)), ]

cc = unique(dds1$condition)

mat = matrix(NA, nrow = nrow(cpm), ncol = length(cc))
colnames(mat) = cc
rownames(mat) = rownames(cpm)

for(n in seq_len(ncol(mat)))
{
  mat[,n] = apply(cpm[, which(dds1$condition == colnames(mat)[n])], 1, median)
}

mat = mat[, c(2,5,8,11,14, 3,6,9,12,15, 1,4,7,10,13)]

mat = mat[!is.na(match(rownames(mat), tfs)), ]

mat1 <- t(scale(t(mat)))

myColor <- viridis::viridis(paletteLength)
#myColor1 <- colorRampPalette(ArchRPalettes$coolwarm)(paletteLength)
#myColor2 <- colorRampPalette(c('lightgray', 'red'))(paletteLength)

myBreaks <- c(seq(min(mat1), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat1)/paletteLength, max(mat1), length.out=floor(paletteLength/2)))

# set annotation
anno_col <- data.frame(cbind(time = rep(c('d0.8h', 'd1', 'd2', 'd2.18h', 'd2.18h'), time = 3),
                             treatment =rep(c('noRA', 'noRA', 'noRA', 'noRA', 'RA'), time = 3)))
rownames(anno_col) <- colnames(mat1)

#col<- colorRampPalette(c("blue4", "white", "darkred"))(paletteLength)
#col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdGy")))(paletteLength)
#col = colorRampPalette(rev(brewer.pal(n = 7, name ="PuOr")))(16)
#col = palette(hcl.colors(8, "Viridis"))
#col = colorRampPalette(c("navy", "white", "red3"))(paletteLength)

library(pheatmap)
pheatmap::pheatmap(
  mat1,
  border_color = 'gray60',
  color = myColor,
  breaks = myBreaks,
  annotation_col = anno_col,
  angle_col = 45,
  cluster_cols = FALSE,
  fontsize_row = 4,
  treeheight_row = 30,
  gaps_col =  c(5, 10),
  filename = paste0(outDir, 'heatmap_tfs.sps_', treatment, '_with.geneNames.pdf'), 
  width = 8, height = 40
)

pheatmap::pheatmap(
  mat1,
  border_color = 'gray60',
  color = myColor,
  breaks = myBreaks,
  annotation_col = anno_col,
  angle_col = 45,
  cluster_cols = FALSE,
  fontsize_row = 4,
  show_rownames = FALSE,
  treeheight_row = 30,
  gaps_col =  c(5, 10),
  filename = paste0(outDir, 'heatmap_tfs.sps_', treatment, '_without.geneNames.pdf'), 
  width = 5, height = 10
)


##########################################
# plot individual genes  
##########################################
# tfs and sps annotations 
require(DESeq2)
load(file = paste0(RdataDir, 'design_dds_', version.analysis, '.Rdata'))
design$treatment = paste0(design$Wash, '_', design$Culture)

#dds <- DESeq(dds, test="LRT", reduced=~1)
#res <- results(dds)
#o1 = order(res$padj)

fpm = fpm(dds)
fpm = log2(fpm + 2^-4)

mm = match(rownames(fpm), tfs)
fpm = fpm[which(!is.na(mm)), ]

write.csv(fpm, file = paste0(resDir, '/tfs_sps_allCondtions_log2cpm.csv'), quote = FALSE,
          row.names = TRUE)

source("functions_utility.R")
cc = unique(design$condition)
treat = design$treatment[match(cc, design$condition)]
cc = gsub(' ', '', cc)
fpm = cal_sample_means(fpm, conds = cc)

jj0 = which(treat == 'ESCs_E8')
jj1 = which(treat == 'mTesr_LDNSB')
jj2 = which(treat == 'mTesr_Adv')
jj3 = which(treat == 'Adv_LDNSB')
jj4 = which(treat == 'mTesr_LDN')
jj5 = which(treat == 'mTesr_SB')

pdfname = paste0(resDir, '/Comparison_TFs_SPs_allconditions.pdf')
pdf(pdfname,  width = 12, height = 8)
#par(cex = 1.0, las = 1, mgp = c(3,2,0), mar = c(6,6,2,0.2), tcl = -0.3)

for(n in seq_len(nrow(fpm)))
#for(n in c(1:10))
{
  # n = 538
  cat(n, ' -- ', rownames(fpm)[n], '\n')
  plot(0, 0, type = 'n', xlim = c(0, 120), ylim = range(c(fpm[n, ])), 
       main = rownames(fpm)[n], 
       ylab = 'log2(CPM)', xlab = 'time')
  
  ## mTesr-LDNSB
  points(c(0, 8, 24, 48, 48+18, 96, 96+18), 
         fpm[n, c(jj0, jj1[-c(5,7)])], col = 'darkblue', type = 'l', pch = 16, 
         lwd = 2.0)
  points(c(0, 8, 24, 48, 48+18, 96, 96+18), 
         fpm[n, c(jj0, jj1[-c(5,7)])], col = 'darkblue', type = 'p', 
         pch = 1, cex = 1.5)
  points(c(48, 48+18), fpm[n, jj1[c(3,5)]], col = 'darkgreen', 
         type = 'l', pch = 16, 
         lwd = 2.0)
  points(c(48, 48+18), fpm[n, jj1[c(3,5)]], col = 'darkgreen', 
         type = 'p', pch = 16, 
         cex = 1.5)
  
  points(c(96, 96+18), fpm[n, jj1[c(6,7)]], col = 'darkgreen', 
         type = 'l', pch = 16, 
         lwd = 2.0)
  points(c(96, 96+18), fpm[n, jj1[c(6,7)]], col = 'darkgreen', 
         type = 'p', pch = 16, 
         cex = 1.5)
  
  # mTesr-Adv
  points(c(0, 8, 24, 48, 48+18), 
         fpm[n, c(jj0, jj2[-c(5)])], col = 'darkorange', type = 'l', 
         pch = 1, lwd = 1.5)
  points(c(0, 8, 24, 48, 48+18), 
         fpm[n, c(jj0, jj2[-c(5)])], col = 'darkorange', type = 'p', 
         pch = 1, cex = 1.5)
  points(c(48, 48+18), fpm[n, jj2[c(3,5)]], col = 'darkgreen', 
         type = 'l', pch = 16, 
         lwd = 2.0)
  points(c(48, 48+18), fpm[n, jj2[c(3,5)]], col = 'darkgreen', 
         type = 'p', pch = 16, 
         cex = 1.5)
  
  # Adv_LDNSB
  points(c(0, 8, 24, 48, 48+18), 
         fpm[n, c(jj0, jj3[-c(5)])], col = 'black', type = 'l', 
         pch = 1, lwd = 1.5)
  points(c(0, 8, 24, 48, 48+18), 
         fpm[n, c(jj0, jj3[-c(5)])], col = 'black', type = 'p', 
         pch = 1, cex = 1.5)
  points(c(48, 48+18), fpm[n, jj3[c(3,5)]], col = 'darkgreen', 
         type = 'l', pch = 16, 
         lwd = 2.0)
  points(c(48, 48+18), fpm[n, jj3[c(3,5)]], col = 'darkgreen', 
         type = 'p', pch = 16, 
         cex = 1.5)
  ## mTesr-LDN 
  points(c(48), fpm[n, c(jj4)], col = 'magenta', 
         type = 'p', pch = 2, 
         cex = 1.5)
  
  points(c(48), fpm[n, jj5], col = 'darkred', 
         type = 'p', pch = 4, 
         cex = 1.5)
  
  legend('topright', 
         legend = c('mTesr-LDNSB', 'mTesr-Adv', 
                               'Adv-LDNSB', 'mTesr-LDN', 
                               'mTesr-SB', 'RA'), 
         bty = 'n', 
         col = c('darkblue', 'darkorange','black', 'magenta', 'red', 
                 'darkgreen'),
         lwd =1.5, 
         pch = c(1, 1, 1,2,4,16), lty = rep(1, 6))
}

dev.off()


########################################################
########################################################
# Section IV: RA target genes and compare them with mouse RA target genes
# 
########################################################
########################################################
require(ggplot2)
library(ggrepel)
require(gridExtra)
library(dplyr)
library(ggrepel)

# import the mouse RA target by comparing the RA vs noRA
RARtargetDir = '../results/RA_targets_L118404_smartseq3_20221117/'
ggs = read.csv2(paste0(RARtargetDir, 'fpm_stat.allpairwseCompares.csv'),
                row.names = c(1))
ggs = ggs[, grep('RA_d2.18hvs.noRA_d2.18h', colnames(ggs))]
ggs = ggs[, c(1:3)] # select the batch 1
colnames(ggs) = gsub('RA_d2.18hvs.noRA_d2.18h', 'mm', colnames(ggs))

ggs = data.frame(ggs, gene = rownames(ggs), stringsAsFactors = FALSE)

#ggs = rownames(ggs)[which(abs(ggs$log2FoldChange_RA_d2.18hvs.noRA_d2.18h)>2 & 
#                            ggs$padj_RA_d2.18hvs.noRA_d2.18h<0.01)]

# import the human RA targets by comparing the RA vs noRA at day2.18h
res = read.csv(file = paste0("../results/R15282_hNT_quantseq_20230614/out_mTesr_LDNSB/", 
                'DE_table_d4.18h_mTesr_LDNSB_d4.250nM_vs._d4_mTesr_LDNSB_no.csv'), row.names = c(1))
res = res[, c(2,5,6,7)]
colnames(res) = paste0(colnames(res), '_hs')

names = toupper(ggs$gene)
mm = match(names, res$gene_hs)

hmc = data.frame(ggs[!is.na(mm), ], res[mm[!is.na(mm)], ], stringsAsFactors = FALSE)
hmc$pvalue_mm = -log10(hmc$pvalue_mm)

hmc = hmc[order(-hmc$pvalue_mm), ]

hmc = hmc[order(-hmc$pvalue_hs), ]

#hmc$pval.combine = apply(as.matrix(hmc[, grep('pvalue_', colnames(hmc))]), 1, 
#                         function(x){abs(x[1]-x[2])/(min(x) + 0.1)})
#hmc = hmc[order(hmc$pval.combine), ]

write.csv(hmc, file = paste0(resDir, '/RAtargets_combined_mouse_human.csv'), row.names = TRUE, quote = FALSE)

sels = which((hmc$padj_mm<0.05 & abs(hmc$log2FoldChange_mm) > 1) | (hmc$padj_hs <0.05 & abs(hmc$lfc_hs) >1))
hmc = hmc[sels, ]
hmc = hmc[order(hmc$pval.combine), ]

#hmc$pvalue_hs = -log10(hmc$pvalue_hs)

ggplot(data = hmc, aes(x = pvalue_mm, y = pvalue_hs, label = gene)) +
  geom_point(size = 1) + 
  labs(title = paste0('-log10 pvalues'), x = 'mNT', y = 'hNT') + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) + 
  xlim(0, 270) + ylim(0, 270) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  theme_bw() +
  #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
  geom_text_repel(data=subset(hmc, pvalue_mm > 10 & pvalue_hs > 10), size = 2, col = 'blue')

ggsave(paste0(resDir, '/RAtargets_comparison_mouse_human_pvalues.pdf'),  width=12, height = 12)

ggplot(data = hmc, aes(x = log2FoldChange_mm, y = lfc_hs, label = gene)) +
  geom_point(size = 1) + 
  labs(title = paste0('log2 fold-change'), x = 'mNT', y = 'hNT') + 
  theme(axis.text.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) + 
  xlim(-7, 10) + ylim(-7, 10) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  theme_bw() +
  #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
  geom_text_repel(data=subset(hmc, log2FoldChange_mm > 1 & lfc_hs > 1), size = 2, col = 'blue') +
  geom_text_repel(data=subset(hmc, log2FoldChange_mm < -1 & lfc_hs < -1), size = 2, col = 'orange') 

ggsave(paste0(resDir, '/RAtargets_comparison_mouse_human_log2FC.pdf'),  width=12, height = 12)


##########################################
# signaling pathway analysis 
##########################################



########################################################
########################################################
# Section V: compare the trajectory of N2B27 (Adv_LDNSB) vs 
# mTesr_LDNSB
########################################################
########################################################
treatment = 'compare_mTesr_LDNSB_vs._Adv_LDNSB'
outDir = paste0(resDir, '/out_', treatment, '/')
if(!dir.exists(outDir)) dir.create(outDir)

jj = which(design$RA == 'no' & (design$treatment == 'Adv_LDNSB'| design$treatment == 'mTesr_LDNSB'))
dds1 = dds[, jj]

dds1$condition <- droplevels(dds1$condition)
dds1$condition = factor(gsub(' ','.', dds1$condition))
dds1 <- estimateDispersions(dds1, fitType = 'parametric')
plotDispEsts(dds1, ymin = 10^-3)

vsd <- varianceStabilizingTransformation(dds1, blind = FALSE)
pca2save = as.data.frame(plotPCA(vsd, intgroup = c('timepoint', 'Wash', 'Culture'), 
                                 returnData = TRUE, ntop = 1000))

pca2save$time = as.factor(pca2save$timepoint)
pca2save$treatment = paste0(pca2save$Wash, '_', pca2save$Culture)

ggplot(data=pca2save, aes(PC1, PC2, label = name, color= time, shape = treatment))  + 
  geom_point(size=3) + 
  scale_shape_manual(values=1:nlevels(pca2save$time)) +
  #geom_text(hjust = 0.5, nudge_y = 0.5, size=3.0) + 
  geom_text_repel(data= pca2save, size = 2) 

ggsave(paste0(outDir, '/PCA_quantseq_samples_', treatment, '.pdf'),  width=16, height = 10)




