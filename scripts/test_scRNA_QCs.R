##########################################################################
##########################################################################
# Project: RA competence with Hannah  
# Script purpose: further QC of scRNA-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Sep 29 10:12:18 2022
##########################################################################
##########################################################################
##########################################
# further QCs
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_lognormamlized_pca_umap_v1.rds'))


p1 = DimPlot(aa,  reduction = 'pca', dims = c(1, 2), cols = cols, 
             label = TRUE, repel = TRUE, group.by = 'condition') + ggtitle("PCA")

p2 = DimPlot(aa,  reduction = 'pca', dims = c(1, 3), cols = cols, 
             label = TRUE, repel = TRUE, group.by = 'condition') + ggtitle("PCA")
p3 = DimPlot(aa,  reduction = 'pca', dims = c(1, 4), cols = cols, 
             label = TRUE, repel = TRUE, group.by = 'condition') + ggtitle("PCA")

p4 = DimPlot(aa,  reduction = 'pca', dims = c(1, 5), cols = cols, 
             label = TRUE, repel = TRUE, group.by = 'condition') + ggtitle("PCA")


(p1 + p2) / (p3 + p4)
ggsave(filename = paste0(resDir, '/PCA_colorCoded_v1.pdf'), width = 14, height = 10)

DimPlot(aa,  reduction = 'pca', dims = c(4, 5), cols = cols, 
        label = TRUE, repel = TRUE, group.by = 'condition') + ggtitle("PCA")
ggsave(filename = paste0(resDir, '/PCA_pc4.vs.pc5.pdf'), width = 10, height = 8)

DimPlot(aa, label = TRUE,  cols = cols, repel = TRUE, group.by = 'condition') + ggtitle("Unsupervised clustering")

ggsave(filename = paste0(resDir, '/umap_colorCoded_v1.pdf'), width = 10, height = 8)

p1 = FeaturePlot(aa, features = 'nCount_RNA')
p2 = FeaturePlot(aa, features = 'nFeature_RNA')
p1 + p2

FeaturePlot(aa, features = 'percent.mt')
FeaturePlot(aa, dims = c(1, 3), reduction = 'pca', features = 'percent.mt')

# features = rownames(aa)[grep('Tubb3|Nkx2-2|Foxa2|Bmp4|Pax6|Shh|Sox2
#                              |Pou5f1|Sox1|Nkx6-1|Pax3|Pax7|Lef1|Zfp703|Arx', rownames(aa))]
ggs = sapply(rownames(aa), function(x) {x = unlist(strsplit(x, '-')); x = x[grep('ENSMUSG', x, invert = TRUE)];
paste0(x, collapse = '-')})

features = c('Pou5f1', 'Sox2', 'Lef1', 'Otx2', 'Zfp703', 'Pax6', 'Foxa2', 'Shh', 'Nkx6-1', 'Nkx2-2', 'Olig2', 
             'Sox1', 'Tubb3', 'Bmp4', 'Bmp7', 'Nog', 'Pax3', 'Pax7', 'Arx')
features = rownames(aa)[!is.na(match(ggs, features))]

FeaturePlot(aa, features = features)
FeaturePlot(aa, features =  "Hba-x-ENSMUSG00000055609")


ggsave(filename = paste0(resDir, '/markerGenes_inUMAP.pdf'), width = 20, height = 20)


FeaturePlot(aa, features = c( "Foxa2-ENSMUSG00000037025", "Bmp4-ENSMUSG00000021835",  "Bmp7-ENSMUSG00000008999", 
                              "Nog-ENSMUSG00000048616"))

Idents(aa) = aa$condition
FeatureScatter(subset(aa, idents = 'day4_RA'), cols = cols, feature1 = "Pax6-ENSMUSG00000027168", 
               feature2 = "Foxa2-ENSMUSG00000037025")

plot1 <- VariableFeaturePlot(aa)
top10 <- head(VariableFeatures(aa), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

