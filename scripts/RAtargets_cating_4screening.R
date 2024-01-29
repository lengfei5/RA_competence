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
sps = unique(sps$gene)

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

write.table(ggs, file = paste0(resDir, '/RAtargets_chipseq.peak.assignment_chiapet.closestTSS.txt'), 
            sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)





