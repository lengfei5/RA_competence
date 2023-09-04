##########################################################################
##########################################################################
# Project: run the GENIE3 for FoxA2 regulators
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Sep  4 13:54:49 2023
##########################################################################
##########################################################################
rm(list =ls())

library(tictoc)

E = readRDS(file = paste0('analysis_examples/Genie3/ExprMatrix_4FoxA2.rds'))

tic()
source('myGENIE3.R')
wtm = GENIE3_withSpecificTargets(expr.matrix = E, priorTargets = c('Foxa2'), ncore = 1)
#saveRDS(wtm, file = paste0(RdataDir, '/first_test_Genie3_v3.rds'))

toc()

wtm = data.frame(wtm, regulator = rownames(wtm), stringsAsFactors = FALSE)
wtm = wtm[order(-wtm$Foxa2), ]

saveRDS(wtm, file = paste0('analysis_examples/Genie3/ranked_predictedRegulators_FoxA2.rds'))

