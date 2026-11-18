##########################################################################
##########################################################################
# Project: RA competence for Hannah's project
# Script purpose: analysis for the first revision
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Nov 17 13:31:48 2026
##########################################################################
##########################################################################

rm(list = ls())

version.analysis = '_R13547_10x_mNT_20240522'
species = 'mNT_scRNAseq'

resDir = paste0("../results/figures_tables", version.analysis)
figureDir = '/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/plots_manuscript_1/'
RdataDir = '../results/Rdata/'
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

library(pryr) # monitor the memory usage
require(ggplot2)

require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(viridis)

#library(DropletUtils)
library(future)
options(future.globals.maxSize = 160000 * 1024^2)
mem_used()

levels = c("day0_beforeRA", "day1_beforeRA", 
           "day2_beforeRA",
           "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", 'day3.5_RA',
           "day4_RA", "day5_RA", "day6_RA",
           "day2.5_noRA", "day3_noRA", 'day3.5_noRA', "day4_noRA", "day5_noRA", "day6_noRA")

cols = readRDS(file = '../results/Rdata/color_scheme_4scRNAseq.rds')
load(file = '../results/Rdata/tfs_sps_geneExamples_4scRNAseq.Rdata')

feat_cols = c("#F0F0F0", "#EFFAB6", "#69C6BE", "#007BB7", "#121D60")

cols_cluster = c( "#7F7F7F", "#FFC000", "#70AD47", "#337f01", "#CD00CF", "#265401","#800080")



##########################################
#  
##########################################


