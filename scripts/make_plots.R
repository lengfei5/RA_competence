##########################################################################
##########################################################################
# Project: RA competence
# Script purpose: make plots
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jan 19 10:59:56 2023
##########################################################################
##########################################################################


########################################################
########################################################
# Section : make plots for Hannah's Monday seminar
# 
########################################################
########################################################
aa =  readRDS(file = paste0(RdataDir, 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', species, version.analysis, '.rds'))

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

#aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 50, min.dist = 0.1)
#aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.1)
# aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 50, min.dist = 0.2)
# aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 50, min.dist = 0.3)
#aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 100, min.dist = 0.2)
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)


DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE) +
  theme(axis.text.x = element_text(angle = 0, size = 14), 
        axis.text.y = element_text(angle = 0, size = 14), 
        axis.title =  element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)
        #legend.position=c(0.2, 0.8),
        #plot.margin = margin()
        #legend.key.size = unit(1, 'cm')
        #legend.key.width= unit(1, 'cm')
  )
ggsave(paste0("../results/plots_MondaySeminar", 
              '/UMAP_condition_toUse.pdf'),  width=8, height = 6) 


DimPlot(aa, cols = cols, group.by = 'condition', label = FALSE, repel = TRUE, raster = FALSE) +
  theme(axis.text.x = element_text(angle = 0, size = 14), 
        axis.text.y = element_text(angle = 0, size = 14), 
        axis.title =  element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14)
        #legend.position=c(0.2, 0.8),
        #plot.margin = margin()
        #legend.key.size = unit(1, 'cm')
        #legend.key.width= unit(1, 'cm')
  )
ggsave(paste0("../results/plots_MondaySeminar", 
              '/UMAP_condition_toUse_noLabel.pdf'),  width=8, height = 6) 


saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                          'cellCycleScoring_annot.v1_savedUMAP.v1_', species, version.analysis, '.rds'))


##########################################
# features and TF activity overlaying UMAP 
##########################################
levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1",  'day3_RA.rep2',
                "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
data_version = "_d2_d2.5_d3_d3.5_d4_d5"

names(cols) = levels
cols_sel = cols[match(levels_sels, names(cols))]

Idents(aa) = factor(aa$condition, levels = levels)
subs = subset(aa, idents = levels_sels)

DimPlot(subs, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)

library(viridis) 
(FeaturePlot(subs, features = c('Foxa2', 'Pax6'), cols = c("lightgrey", "blue"))) 
#&
#  scale_colour_viridis_c(option = "D")

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_Foxa2_Pax6.pdf'),  width=12, height = 6) 

# rerun the umap 
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(subs, ndims = 50)
Idents(subs) = subs$condition

subs <- RunUMAP(subs, dims = 1:30, n.neighbors = 100, min.dist = 0.2)

DimPlot(subs, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE, pt.size = 1.5)

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_updatedUMAP.pdf'),  width=8, height = 6) 

DimPlot(subs, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE, pt.size = 1.5) +
  NoLegend()

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_updatedUMAP.pdf'),  width=7, height = 6) 


DimPlot(subs, cols = cols, group.by = 'condition', label = FALSE, repel = TRUE, raster = FALSE, pt.size = 1.5) +
  NoLegend()

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_updatedUMAP.pdf'),  width=7, height = 6) 

FeaturePlot(subs, features = c('Foxa2', 'Pax6'))

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_Foxa2_Pax6_updatedUMAP.pdf'),  width=12, height = 6) 

FeaturePlot(object = subs,  features = c('Foxa2', 'Pax6'), blend = TRUE,
            cols =  c("lightgrey","#00ff00",  "magenta"))

ggsave(paste0("../results/plots_MondaySeminar", 
              '/umap_d2.to.d6_Foxa2_Pax6_updatedUMAP_blend.pdf'),  width=24, height = 6) 

saveRDS(subs, file = paste0(RdataDir, 
              'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
              'cellCycleScoring_annot.v1_savedUMAP.subs.v2_', species, version.analysis, '.rds'))

ggs = c('Cdh1', 'Crabp2', 'Zeb2', 'Rbpj', 'Gli1', 'Shh', 'Tead1', 'Tead2', 'Tead3', "Tead4",
        "Gli2", "Esr1", "Kdm5b", 'Lef1', 'Otx2','Pou5f1', 'Sox1', 'Nanog', 'Nkx6-1', 
        'Nkx2-2', "Olig2", "Pax3", "Pax7", "Arx", "Sox10", "Pbx3", "Prdm1",  "Sox2", "Tcf7l1",
        "Tcf7l2", "Zeb1", "Xbp1", "Tfdp1","Pax6", "Foxa2", "Smad1", "Smad4", "Smad2", "Smad5", "Smad3", "Smad7", 
        "Smad9", "Bmp4", 'Bmpr2', 'Bmpr1b', 'Bmpr1a', "Acvr2a", "Acvr2b", "Smo", "Notch1", 
        "Lrp5", "Lrp6", "Fzd4", 'Fzd8', 'Wnt6', "Wnt5a", "Wnt5b", "Wnt4", "Wnt1", "Wnt3a", "Ptch1",
        "Ptch2", "Lrp2", "Gpc1", "Gas1", "Cdon", "Hhip", "Fgf8", "Cdh2")
mm = match(ggs, rownames(subs))
ggs[which(is.na(mm))]

ggs = ggs[which(!is.na(mm))]

for(g in ggs)
{
  # g = "Foxa2"
  cat(g, "-- \n")
  FeaturePlot(object = subs,  features = g)
  
  ggsave(paste0("../results/plots_MondaySeminar", 
                '/umap_d2.to.d6_Foxa2_Pax6_updatedUMAP_geneExpression.example_', g, '_withLabel.pdf'),
         
         width=8, height = 6) 
  
  FeaturePlot(object = subs,  features = g) + NoLegend()
  
  ggsave(paste0("../results/plots_MondaySeminar", 
                '/umap_d2.to.d6_Foxa2_Pax6_updatedUMAP_geneExpression.example_', g, '.pdf'),
         width=7, height = 6) 
  
}


aa = readRDS(file = paste0(RdataDir, 
                   'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                   'cellCycleScoring_annot.v1_savedUMAP.subs.v2_savedTFactivities.decoupleR.rds'))

DefaultAssay(object = aa) <- "tfswmean"

ggs = c('Cdh1', 'Crabp2', 'Zeb2', 'Rbpj', 'Gli1', 'Shh', 'Tead1', 'Tead2', 'Tead3', "Tead4",
        "Gli2", "Esr1", "Kdm5b", 'Lef1', 'Otx2','Pou5f1', 'Sox1', 'Nanog', 'Nkx6-1', 
        'Nkx2-2', "Olig2", "Pax3", "Pax7", "Arx", "Sox10", "Pbx3", "Prdm1",  "Sox2", "Tcf7l1",
        "Tcf7l2", "Zeb1", "Xbp1", "Tfdp1","Pax6", "Foxa2", "Smad1", "Smad4", "Smad2", "Smad5", "Smad3", "Smad7", 
        "Smad9", "Bmp4", 'Bmpr2', 'Bmpr1b', 'Bmpr1a', "Acvr2a", "Acvr2b", "Smo", "Notch1", 
        "Lrp5", "Lrp6", "Fzd4", 'Fzd8', 'Wnt6', "Wnt5a", "Wnt5b", "Wnt4", "Wnt1", "Wnt3a", "Ptch1",
        "Ptch2", "Lrp2", "Gpc1", "Gas1", "Cdon", "Hhip", "Fgf8", "Cdh2")
mm = match(ggs, rownames(aa))
ggs[which(is.na(mm))]

ggs = ggs[which(!is.na(mm))]

for(gene in ggs)
{
  # gene = "Foxa2"
  cat("--", gene, "-- \n")
  #FeaturePlot(object = aa,  features = g)
  (FeaturePlot(aa, features = gene) & 
      scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) + 
    ggtitle(paste0(gene, ' activity'))
  ggsave(paste0("../results/plots_MondaySeminar/TF_activity", 
                '/umap_d2.to.d6_TFactivity_geneExamples_', gene, '_withLabel.pdf'),
         width=8, height = 6) 
  
  (FeaturePlot(aa, features = gene) & 
      scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) + 
    ggtitle(paste0(gene, ' activity')) + NoLegend()
  
  ggsave(paste0("../results/plots_MondaySeminar/TF_activity", 
                '/umap_d2.to.d6_TFactivity_geneExamples_', gene, '.pdf'),
         width=7, height = 6)
  
}

