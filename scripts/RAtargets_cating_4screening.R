##########################################################################
##########################################################################
# Project: RA competence projects
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Jan 29 10:28:17 2024
##########################################################################
##########################################################################
rm(list = ls())

RNA.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')


# setup for data import and sequencing QCs
version.analysis = '_concatenation_4CrisprScreening/'

resDir = paste0("../results/RA_targets", version.analysis)
RdataDir = paste0(resDir, 'Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)


## import the tfs and sps annotation 
sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
#sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

##########################################
# RA direct target from ChIP-seq and ChIA-pet
##########################################
RARtargetDir = '../results/RAR_targets/Rdata/'
peaks = readRDS(file = paste0(RARtargetDir, 
                              'RAR_chipseq.peak.assignment_promoters.genebody.downstream.chiapet.closestTSS.rds'))

ggs = unlist(peaks[, c(13:17)])
ggs = ggs[!is.na(ggs)]
ggs = unique(ggs)

xx = c()
for(n in 1:length(ggs))
{
  xx = c(xx, unlist(strsplit(as.character(ggs[n]), ';')))
}
ggs = unique(xx)

saveRDS(ggs, file = paste0(RdataDir, 'RAR_targets_chip_chiapet.rds'))

write.csv2(ggs, file = paste0(resDir, '/RAtargets_chipseq.peak.assignment_chiapet.closestTSS.csv'), 
            row.names = FALSE, quote = FALSE)


########################################################
########################################################
# Section I :  
# 
########################################################
########################################################

##########################################
# Elena Gromberg bulk RNA-seq 
##########################################
tableDir = '../results/RAR_targets'

res = read.csv(file = paste0(tableDir, '/DESeq2_DEgenes_fdr.0.05_normalizedData.csv'), row.names = 1,
               header = TRUE)

res$padj = 10^(-res$padj)

## upregulated
xx = res[which(res$log2FoldChange > 1), ]
xx$tfs = NA
xx$sps = NA
xx$tfs[which(!is.na(match(xx$gene, tfs)))] = 'TF'
mm = match(xx$gene, sps$gene)
xx$sps[which(!is.na(mm))] = sps$pathway[mm[which(!is.na(mm))]]

write.csv2(xx, file = paste0(resDir, 'RAtargets_geneLists_Gromberg.bulkRNAseq_upregulated.csv'), 
          row.names = TRUE, quote = FALSE)

mm = match(xx$gene, ggs)
xx = xx[which(!is.na(mm)), ]

write.csv2(xx, file = paste0(resDir, 'RAtargets_geneLists_Gromberg.bulkRNAseq_upregulated_intersectedPeaks.csv'), 
           row.names = TRUE, quote = FALSE)

## downreuglated 
## upregulated
xx = res[which(res$log2FoldChange < (-1)), ]
xx$tfs = NA
xx$sps = NA
xx$tfs[which(!is.na(match(xx$gene, tfs)))] = 'TF'
mm = match(xx$gene, sps$gene)
xx$sps[which(!is.na(mm))] = sps$pathway[mm[which(!is.na(mm))]]

write.csv2(xx, file = paste0(resDir, 'RAtargets_geneLists_Gromberg.bulkRNAseq_downregulated.csv'), 
           row.names = TRUE, quote = FALSE)

mm = match(xx$gene, ggs)
xx = xx[which(!is.na(mm)), ]

write.csv2(xx, file = paste0(resDir, 'RAtargets_geneLists_Gromberg.bulkRNAseq_downregulated_intersectedPeaks.csv'), 
           row.names = TRUE, quote = FALSE)


##########################################
# Hannah's smart-seq3 data 
##########################################
tableDir = '../results/RA_targets_L118404_smartseq3_20221117/Compare.diffRAstimulationTime.Batch1/all'

res = read.csv2(file = paste0(tableDir, '/DESeq2_DEgenes_pairwiseComparison_RA_d2.18h.vs.noRA_d2.18h.csv'), 
               row.names = 1, header = TRUE)

#res$padj = 10^(-res$padj)
res$gene = rownames(res)

## upregulated
xx = res[which(res$log2FoldChange_RA_d2.18hvs.noRA_d2.18h > 1), ]
xx$tfs = NA
xx$sps = NA
xx$tfs[which(!is.na(match(xx$gene, tfs)))] = 'TF'
mm = match(xx$gene, sps$gene)
xx$sps[which(!is.na(mm))] = sps$pathway[mm[which(!is.na(mm))]]

write.csv2(xx, file = paste0(resDir, 'RAtargets_geneLists_Hannah.Laura.bulkRNAseq_upregulated.csv'), 
           row.names = TRUE, quote = FALSE)

mm = match(xx$gene, ggs)
xx = xx[which(!is.na(mm)), ]

write.csv2(xx, file = paste0(resDir, 
                             'RAtargets_geneLists_Hannah.Laura..bulkRNAseq_upregulated_intersectedPeaks.csv'), 
           row.names = TRUE, quote = FALSE)

## downreuglated 
## upregulated
xx = res[which(res$log2FoldChange < (-1)), ]
xx$tfs = NA
xx$sps = NA
xx$tfs[which(!is.na(match(xx$gene, tfs)))] = 'TF'
mm = match(xx$gene, sps$gene)
xx$sps[which(!is.na(mm))] = sps$pathway[mm[which(!is.na(mm))]]

write.csv2(xx, file = paste0(resDir, 'RAtargets_geneLists_Hannah.Laura.bulkRNAseq_downregulated.csv'), 
           row.names = TRUE, quote = FALSE)

mm = match(xx$gene, ggs)
xx = xx[which(!is.na(mm)), ]

write.csv2(xx, file = paste0(resDir, 
                             'RAtargets_geneLists_Hannah.Laura..bulkRNAseq_downregulated_intersectedPeaks.csv'), 
           row.names = TRUE, quote = FALSE)


##########################################
# scRNA-seq data for RA responsive targets and RA withdraw targets  
##########################################
require(Seurat)
aa = readRDS(file = paste0("../results/Rdata/", 
                           "seuratObject_RA.vs.noRA.bifurcation_doublet.rm_mt.ribo.filtered_regressout.nCounts",
                           "_cellCycleScoring_annot.v1_reduction.DM_princurves_",
                           "mNT_scRNAseq_R13547_10x_mNT_20220813.rds"))

Idents(aa) = factor(aa$condition)

res = FindMarkers(aa, 
                  ident.1 = "day2.5_RA",
                  ident.2 = 'day2.5_noRA',
                  only.pos = FALSE, 
                  test.use = "wilcox",
                  min.pct = 0.1, 
                  logfc.threshold = 0.25)

res$gene = rownames(res)

# upregulated
xx = res[which(res$avg_log2FC>0), ]
xx$tfs = NA
xx$sps = NA
xx$tfs[which(!is.na(match(xx$gene, tfs)))] = 'TF'
mm = match(xx$gene, sps$gene)
xx$sps[which(!is.na(mm))] = sps$pathway[mm[which(!is.na(mm))]]

write.csv2(xx, file = paste0(resDir, 'RAtargets_geneLists_scRNAseq_upregulated.csv'), 
           row.names = TRUE, quote = FALSE)

mm = match(xx$gene, ggs)
xx = xx[which(!is.na(mm)), ]

write.csv2(xx, file = paste0(resDir, 
                             'RAtargets_geneLists_scRNAseq_upregulated_intersectedPeaks.csv'), 
           row.names = TRUE, quote = FALSE)

# downregulated
xx = res[which(res$avg_log2FC < 0), ]
xx$tfs = NA
xx$sps = NA
xx$tfs[which(!is.na(match(xx$gene, tfs)))] = 'TF'
mm = match(xx$gene, sps$gene)
xx$sps[which(!is.na(mm))] = sps$pathway[mm[which(!is.na(mm))]]

write.csv2(xx, file = paste0(resDir, 'RAtargets_geneLists_scRNAseq_downregulated.csv'), 
           row.names = TRUE, quote = FALSE)

mm = match(xx$gene, ggs)
xx = xx[which(!is.na(mm)), ]

write.csv2(xx, file = paste0(resDir, 
                             'RAtargets_geneLists_scRNAseq_downregulated_intersectedPeaks.csv'), 
           row.names = TRUE, quote = FALSE)

### RA withdraw targets
res = FindMarkers(aa, 
                  ident.1 = "day3_RA.rep1",
                  ident.2 = 'day2.5_RA',
                  only.pos = FALSE, 
                  test.use = "wilcox",
                  min.pct = 0.1, 
                  logfc.threshold = 0.25)

res$gene = rownames(res)

# upregulated
xx = res[which(res$avg_log2FC>0), ]
xx$tfs = NA
xx$sps = NA
xx$tfs[which(!is.na(match(xx$gene, tfs)))] = 'TF'
mm = match(xx$gene, sps$gene)
xx$sps[which(!is.na(mm))] = sps$pathway[mm[which(!is.na(mm))]]

write.csv2(xx, file = paste0(resDir, 'RA_withdraw_geneLists_scRNAseq_upregulated.csv'), 
           row.names = TRUE, quote = FALSE)

