##########################################################################
##########################################################################
# Project: RA competence project
# Script purpose: detect doublet in scRNA-seq
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Aug 30 23:41:39 2022
##########################################################################
##########################################################################
library(DoubletFinder)
require(Seurat)

## import the pre-cleaned seurat object
aa = readRDS(file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_', species, version.analysis, '.rds'))

aa$condition[which(aa$condition == 'day0')] = 'day0_beforeRA'
levels = c("day0_beforeRA", "day1_beforeRA", 
           "day2_beforeRA",
           "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", 'day3.5_RA',
           "day4_RA", "day5_RA", "day6_RA",
           "day2.5_noRA", "day3_noRA", 'day3.5_noRA', "day4_noRA", "day5_noRA", "day6_noRA")

aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

cc = unique(aa$condition)

for(n in 1:length(cc))
{
  # n = 1
  subs <- subset(aa, condition == cc[n])
  
  subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 2000)
  subs <- ScaleData(subs)
  
  subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = TRUE)
  
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
  ggsave(filename = paste0(resDir, '/subs_doubletFinder_out_', cc[n], '.pdf'), 
         width = 12, height = 8)
  
  saveRDS(subs, file = paste0(RdataDir, 'subs_doubletFinder_out_', cc[n], '.rds'))
  
}

##########################################
# save the doubletFinder in the main table  
##########################################
cc = unique(aa$condition)
aa$DF_out = NA

for(n in 1:length(cc))
{
  # n = 1
  cat(n, '--', cc[n], '\n')
  subs = readRDS(file = paste0(RdataDir, 'subs_doubletFinder_out_', cc[n], '.rds'))
  aa$DF_out[match(colnames(subs), colnames(aa))] = subs$DF_out
  
}

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_doubletFinderOut_', 
                          species, version.analysis, '.rds'))

##########################################
# Visulize the doublet 
##########################################
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)

ElbowPlot(aa, ndims = 30)

Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'DF_out', raster=FALSE)
ggsave(filename = paste0(resDir, '/umap_doubletFinder_results.pdf'), width = 12, height = 8)

VlnPlot(aa, features = 'nCount_RNA', group.by = 'DF_out', pt.size = 0.1) +
  geom_hline(yintercept = c(25000, 12500), col = 'red')

ggsave(filename = paste0(resDir, '/nCounts_RNA_double.vs.singlet_doubletFinder_results.pdf'), 
       width = 12, height = 8)

as_tibble(data.frame(condition = aa$condition, group= aa$DF_out)) %>%
  group_by(condition, group) %>% tally() 

pcts = c()
for(n in 1:length(cc))
{
  # n =1
  pcts = c(pcts, length(which(aa$DF_out== 'Doublet' & aa$condition == cc[n]))/length(which(aa$condition == cc[n])))
  
}

data.frame(condition = cc, pct = pcts) %>%
  ggplot(aes(x = condition, y = pct, fill = condition)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none")  + 
  ggtitle('pct of doublets by DF ') + 
  theme(axis.text.x = element_text(angle = 90)) 

ggsave(filename = paste0(resDir, '/Percentages_doublet.vs.total_doubletFinder_results.pdf'), 
       width = 12, height = 8)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_doubletFinderOut.v2_', 
                          species, version.analysis, '.rds'))

