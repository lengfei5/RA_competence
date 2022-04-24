##########################################################################
##########################################################################
# Project: RA competence project
# Script purpose: analyze FACS 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Sat Apr 23 23:09:17 2022
##########################################################################
##########################################################################
rm(list = ls())
resDir = '../results/FACS_analysis'
RdataDir = paste0(resDir, '/Rdata')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/groups/tanaka/People/current/jiwang/projects/RA_competence/FACS_timecourse/'

########################################################
########################################################
# Section : metadata and data processing 
# 
########################################################
########################################################
metadata = read.delim(file = paste0(dataDir, 'metadata.txt'), sep = '\t', header = FALSE)

colnames(metadata) = 'fileName'
metadata$fileName = gsub('.fcs', '', metadata$fileName)
metadata$condition = sapply(metadata$fileName, function(x) unlist(strsplit(as.character(x), '[#]'))[1])
metadata$sample = sapply(metadata$fileName, function(x) unlist(strsplit(as.character(x), '[#]'))[2])
metadata$sampleID = sapply(metadata$sample, function(x) unlist(strsplit(as.character(x), '_'))[4])

metadata$time = metadata$condition
metadata = metadata[, c(4, 2, 5, 1, 3)]

kk = grep('_noRA', metadata$condition)
metadata$condition[kk] = 'noRA'

kk = grep('_RA', metadata$condition)
metadata$condition[kk] = 'RA'

kk = grep('d2', metadata$condition)
metadata$condition[kk] = 'beforeRA'

metadata$sampleID = gsub('Sample', '', metadata$sampleID)

metadata$time = gsub('noRA_', '', metadata$time)
metadata$time = gsub('RA_', '', metadata$time)

metadata$time = gsub('_noRA', '', metadata$time)
metadata$time = gsub('_RA', '', metadata$time)

metadata = metadata[order(as.integer(metadata$sampleID)), ]

saveRDS(metadata, file = paste0(RdataDir, '/metadata.rds'))

##########################################
# data processing 
##########################################
aa = read.csv(file = paste0(dataDir, 'AIgood_AllSingleLive_AllParam_1_ScaleVal.csv'), header = TRUE)

jj = grep('Comp.*.H', colnames(aa))
aa = aa[, c(1, jj, 2:5, 31:32)]

saveRDS(aa, file = paste0(RdataDir, '/facs_data.rds'))

########################################################
########################################################
# Section : QCs and processing
# 
########################################################
########################################################



