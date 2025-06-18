##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: predict the signaling pathways controlling the cell proportions 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Feb 17 14:25:06 2025
##########################################################################
##########################################################################
# only RA samples incl. dya2_beforeRA
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_',
                           'regressout.nCounts_',
                           'cellCycleScoring_annot.v1_savedUMAP.subs.v2_', 
                           species, version.analysis, '.rds'))

Idents(aa) = aa$condition
DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)

#ElbowPlot(aa, ndims = 50)
aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

ggsave(filename = paste0(outDir, '/UMAP_RA_symmetryBreaking.pdf'), 
       width = 14, height = 6)

levels_sels = c("day3.5_RA", "day4_RA")
data_version = "_d3.5_d4"

outDir = paste0(resDir, '/RA_symetryBreaking/signaling_pathway_cellProportions', data_version)
system(paste0('mkdir -p ', outDir))

ElbowPlot(aa, ndims = 50)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

ggsave(filename = paste0(outDir, '/UMAP_RA_symmetryBreaking.pdf'), 
       width = 14, height = 6)


##########################################
# test the scFates outcome to select the cells as senders and recevers 
##########################################
library(data.table)
library(plyr)
library(ggplot2)
library(scales)

sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))


dataDir = paste0("../results/scRNAseq_R13547_10x_mNT_20220813/RA_symetryBreaking/TF_modules/",
                 'd2.5_d5_TFs_SPs_regressed.CellCycle_v1/')

outDir = paste0("../results/scRNAseq_R13547_10x_mNT_20220813/RA_symetryBreaking/", 
                "signaling_pathway_cellProportions/")

aa = readRDS(file = paste0(dataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered',
                           '_regressout.nCounts_cellCycleScoring_annot.v2_newUMAP_clusters_time_',
                           'd2.5_d5_regressed.CellCycle_v1.rds'))

#fit = read.csv(file = paste0(dataDir, 'annData_layer_fitted.csv'), header = TRUE, row.names = c(1))
# fit = read.csv(file = paste0(dataDir, 'annData_magic_impuated.csv'), 
#                header = TRUE, row.names = c(1))
# 
# USE_MAGIC_allGenes = TRUE
# if(USE_MAGIC_allGenes){
#   fit.all = fread(file = paste0(dataDir, 'annData_magic_impuated_allGenes.csv'), header = TRUE)
#   fit.all = data.frame(fit.all)
#   rownames(fit.all) = fit.all$V1
#   fit.all = fit.all[, -1]
#   fit.sel = fit.all[, match(colnames(fit), colnames(fit.all))]
#   #fit.sel = data.frame(fit.sel)
#   fit = fit.sel
#   
#   rm(fit.all)
#   rm(fit.sel)
#   #rownames(impuated) = rownames(fit)
#   #colnames(impuated) = colnames(fit)
#   #fit = impuated
# }

pst = read.csv(file = paste0(dataDir, 'annData_pseudotime_segments_milestones.csv'), header = TRUE,
               row.names = c(1))

levels_sels = unique(aa$condition)

names(cols) = levels
cols_sel = cols[match(levels_sels, names(cols))]

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

ggsave(filename = paste0(outDir, '/UMAP_RA_scFates_used_conditions_clusters.pdf'), 
       width = 14, height = 6)


mm = match(colnames(aa), rownames(pst))

aa = AddMetaData(aa, metadata = pst[mm, ])

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE, cols = cols_sel)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seg', raster=FALSE)
p1 + p2

assign = read.csv(file = paste0(dataDir, 'annData_cellAssignment_to_nonIntersectingWindows_8.csv'),
                  header = TRUE, row.names = c(1))
#assign = t(assign)
assignment = rownames(assign)[apply(assign, 2, which.max)] 
names(assignment) = gsub('[.]','-', colnames(assign))

aa$windows = assignment[match(colnames(aa), names(assignment))]

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'milestones', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'windows', raster=FALSE)
p1 + p2

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'windows', raster=FALSE)
ggsave(filename = paste0(outDir, '/UMAP_RA_scFates_used_windows.pdf'), 
       width = 8, height = 6)

FeaturePlot(aa, features = c('Pax6', 'Foxa2'))
ggsave(filename = paste0(outDir, '/UMAP_RA_scFates_used_windows_Pax6_FoxA2.pdf'), 
       width = 14, height = 6)


##########################################
# ligand expression comparisons between cluster 4 and 6  
##########################################
Idents(aa) = as.factor(aa$windows)

markers = FindMarkers(aa,
                         ident.1 = 6, ident.2 = 4,
                         only.pos = FALSE,   
                         min.pct = 0.1,
                         logfc.threshold = 0.1)


#DotPlot(aa, idents = c('6', '4', '3'), features = features) + RotatedAxis()
#DotPlot(pbmc3k.final, features = features, split.by = "groups") + RotatedAxis()


features = sps$gene[which(sps$pathway == "WNT")]
features = features[which(!is.na(match(features, rownames(aa))))]

kk = match(features, rownames(markers))
features = features[which(!is.na(kk))]

xx = markers[which(!is.na(match(rownames(markers), features))), ]
xx = xx[grep('Wnt|Dkk|Tcf', rownames(xx)), ]

xx = data.frame(gene = rownames(xx), log2fc = xx$avg_log2FC)

xx$sign = 'positive'
xx$sign[which(xx$log2fc<0)] = 'negative'

xx = xx[c(grep('Wnt', xx$gene), 
          grep('Dkk', xx$gene),
          grep('Tcf', xx$gene)), ]

ggplot(xx, aes(log2fc, gene)) +
  geom_bar(stat = "identity", aes(fill = sign)) + 
  scale_fill_manual(values = c("darkgray", "green3")) + 
  scale_x_continuous(limits = c(-1, 1)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 10, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  labs( x = 'log2 fold-change', y = '' )

ggsave(filename = paste0(outDir, '/barplot_WNT.pdf'), 
       width = 8, height = 10)


features = sps$gene[which(sps$pathway == "BMP")]
features = features[which(!is.na(match(features, rownames(aa))))]

kk = match(features, rownames(markers))
features = features[which(!is.na(kk))]

xx = markers[which(!is.na(match(rownames(markers), features))), ]
xx = xx[grep('Bmp|Smad|Nogg|Run|Id|Fst', rownames(xx)), ]

xx = data.frame(gene = rownames(xx), log2fc = xx$avg_log2FC)

xx$sign = 'positive'
xx$sign[which(xx$log2fc<0)] = 'negative'

# xx = xx[c(grep('Wnt', xx$gene), 
#           grep('Dkk', xx$gene),
#           grep('Tcf', xx$gene)), ]

ggplot(xx, aes(log2fc, gene)) +
  geom_bar(stat = "identity", aes(fill = sign)) + 
  scale_fill_manual(values = c("darkgray", "green3")) + 
  scale_x_continuous(limits = c(-1.5, 1.5)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 10, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  labs( x = 'log2 fold-change', y = '' )

ggsave(filename = paste0(outDir, '/barplot_BMP.pdf'), 
       width = 8, height = 10)


features = sps$gene[which(sps$pathway == "FGF")]
features = features[which(!is.na(match(features, rownames(aa))))]

kk = match(features, rownames(markers))
features = features[which(!is.na(kk))]

xx = markers[which(!is.na(match(rownames(markers), features))), ]
xx = xx[grep('Dusp|Fgf|Mgpk|Spry|En2', rownames(xx)), ]

xx = data.frame(gene = rownames(xx), log2fc = xx$avg_log2FC)

xx$sign = 'positive'
xx$sign[which(xx$log2fc<0)] = 'negative'

# xx = xx[c(grep('Wnt', xx$gene), 
#           grep('Dkk', xx$gene),
#           grep('Tcf', xx$gene)), ]

ggplot(xx, aes(log2fc, gene)) +
  geom_bar(stat = "identity", aes(fill = sign)) + 
  scale_fill_manual(values = c("darkgray", "green3")) + 
  scale_x_continuous(limits = c(-1.5, 1.5)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, size = 10, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14)) +
  labs( x = 'log2 fold-change', y = '' )

ggsave(filename = paste0(outDir, '/barplot_FGF.pdf'), 
       width = 8, height = 10)



##########################################
# LIANA test 
##########################################
require(liana)
require(Seurat)
require(scater)
require(scran)

subref = subset(aa, cells = colnames(aa)[which(aa$windows == '6'|aa$windows == '4')]); 

additionalLabel = '_fixedCelltypes'; 
subref$celltypes = subref$windows
celltypes = unique(subref$celltypes)

system(paste0('mkdir -p ', paste0(outDir, '/LR_analysis_LIANA')))
# source('functions_scRNAseq.R') 

sce <- as.SingleCellExperiment(subref)
colLabels(sce) = as.factor(sce$celltypes)
rownames(sce) = toupper(rownames(sce))

ave.counts <- calculateAverage(sce, assay.type = "counts")

#hist(log10(ave.counts), breaks=100, main="", col="grey80",
#     xlab=expression(Log[10]~"average count"))

num.cells <- nexprs(sce, byrow=TRUE)
#smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
#              xlab=expression(Log[10]~"average count"))

# detected in >= 5 cells, ave.counts >=5 but not too high
genes.to.keep <- num.cells > 5 & ave.counts >= 10^-4  & ave.counts <10^2  
summary(genes.to.keep)

sce <- sce[genes.to.keep, ]

## run the liana wrap function by specifying resource and methods
# Resource currently included in OmniPathR (and hence `liana`) include:
show_resources()
# Resource currently included in OmniPathR (and hence `liana`) include:
show_methods()

liana_test <- liana_wrap(sce,  
                         method = c("natmi", "connectome", "logfc", "sca", "cytotalk"  
                                    #'cellphonedb'
                         ),
                         resource = c("Consensus", 'CellPhoneDB', "OmniPath", "LRdb",
                                      "CellChatDB",  "CellTalkDB"), 
                         assay.type = "logcounts", 
                         idents_col = 'celltypes')

# Liana returns a list of results, each element of which corresponds to a method
liana_test %>% glimpse

# We can aggregate these results into a tibble with consensus ranks
saveRDS(liana_test, file = paste0(outDir, '/LR_analysis_LIANA/res_lianaTest_Consensus', 
                                  additionalLabel, '.rds'))

liana_test = readRDS(file = paste0(outDir, '/LR_analysis_LIANA/res_lianaTest_Consensus', 
                                   additionalLabel, '.rds'))

liana_test <- liana_test %>%
  liana_aggregate(resource = 'Consensus')

if(is.na(receiver_cells)){ # loop over all cell type candidates
  celltypes = as.character(celltypes)
  receiver_cells = as.character(celltypes)
}

ntop = 300
#' manually_specifying_sender_receiver = FALSE
#' if(manually_specifying_sender_receiver){
#'   ntop = 100
#'   sender_cells = c('3', '4')
#'   receiver_cells = sender_cells
#'   
#'   liana_test %>%
#'     liana_dotplot(source_groups = sender_cells,
#'                   target_groups = receiver_cells,
#'                   ntop = ntop)
#'   ggsave(filename = paste0(outDir, '/LR_analysis_LIANA_v2/liana_LR_prediction_cluster3.vs.cluster4', 
#'                            additionalLabel, 
#'                            #'_receiverCells.', receiver_cells[m], 
#'                            '_ntop.', ntop, '.pdf'), 
#'          width = 15, height = 0.25*ntop, limitsize = FALSE)
#'   
#'   ntop = 500
#'   sender_cells = c('0', '6', '2', '8')
#'   receiver_cells = sender_cells
#'   
#'   liana_test %>%
#'     liana_dotplot(source_groups = sender_cells,
#'                   target_groups = receiver_cells,
#'                   ntop = ntop)
#'   ggsave(filename = paste0(outDir, '/LR_analysis_LIANA_v2/liana_LR_prediction_cluster_0_6_2_8', 
#'                            additionalLabel, 
#'                            #'_receiverCells.', receiver_cells[m], 
#'                            '_ntop.', ntop, '.pdf'), 
#'          width = 25, height = 0.25*ntop, limitsize = FALSE)
#'   
#'   
#' }

for(m in 1:length(receiver_cells))
{
  # m = 1
  cat(m, '-- receiver cells : ', receiver_cells[m], '\n')
  #liana_test %>%
  #  liana_dotplot(source_groups = celltypes[n],
  #                target_groups = celltypes,
  #                ntop = ntop)
  liana_test %>%
    liana_dotplot(source_groups = celltypes,
                  target_groups = receiver_cells[m],
                  ntop = ntop)
  #liana_test_save =  liana_test %>% filter()
  #  liana_dotplot(source_groups = celltypes,
  #                target_groups = receiver_cells[m],
  #                ntop = ntop)
  
  ggsave(filename = paste0(outDir, '/LR_analysis_LIANA/liana_LR_prediction_recieveCell', 
                           additionalLabel, 
                           '_receiverCells.', receiver_cells[m], 
                           '_ntop.', ntop, '.pdf'), 
         width = 30, height = 0.25*ntop, limitsize = FALSE)
  
}


liana_test = readRDS(file = paste0(outDir, '/LR_analysis_LIANA/res_lianaTest_Consensus', 
                                   additionalLabel, '.rds'))

liana_test <- liana_test %>%
  liana_aggregate(resource = 'Consensus')

celltypes = unique(subref$celltypes)

receivers = celltypes

df_test = liana_test %>% filter(target %in% receivers & source %in% celltypes) %>% as.data.frame() 

write.table(df_test, file = paste0(outDir, '/res_lianaTest_Consensus', additionalLabel, '.txt'), 
            sep = '\t', quote = FALSE)

saveRDS(liana_test, file = paste0(outDir, '/res_lianaTest_Consensus', 
                                  additionalLabel, '_saved.rds'))


res = df_test
res = res[, c(1:5, which(colnames(res) == 'natmi.edge_specificity'), 
              which(colnames(res) == 'sca.LRscore'))]

colnames(res)[1:4] = c('sender', 'receiver', 'ligand', 'receptor')

#require(cellcall)
library(SeuratData)
library(Connectome)
library(cowplot)

#res = res[order(-res$sca.LRscore), ]

colnames(res)[1:2] = c('source', 'target')
res$weight_norm = res$sca.LRscore
res$pair = paste0(res$ligand, ' - ', res$receptor)
res$vector = paste0(res$source, ' - ', res$target)
res$edge = paste0(res$source, ' - ', res$ligand, ' - ', res$receptor, ' - ', res$target)
res$source.ligand = paste0(res$source, ' - ', res$ligand)
res$receptor.target = paste0(res$receptor, ' - ', res$target)

write.table(res, 
            file = paste0(outDir, '/LR_interactions_allPairs_LIANA.txt'), 
            quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')

saveRDS(res, file = paste0(outDir, '/res_lianaTest_for_circosplot.rds'))

##########################################
### plot circosplot
##########################################
#res = readRDS(file = paste0(outDir, '/res_lianaTest_for_circosplot.rds'))
res = readRDS(file = paste0(outDir, '/res_lianaTest_for_circosplot.rds'))

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
#res = readRDS(file = paste0(outDir, '/res_lianaTest_for_circosplot.rds'))
source(paste0(functionDir, '/functions_cccInference.R'))

celltypes = unique(c(res$source, res$target))
additionalLabel = '_fixedCelltypes'
receivers = celltypes

sender_cells = receivers
receiver_cells = receivers

print(as.character(sender_cells))
print(as.character(receiver_cells))


head(grep('BMP', res$ligand))
#res = res[!is.na(match(res$source, sender_cells)) & !is.na(match(res$target, receiver_cells)), ]

res = res[order(res$aggregate_rank), ]

head(grep('BMP', res$ligand))

#cell_color = randomcoloR::distinctColorPalette(length(cells.of.interest))

cells.of.interest = unique(c(res$source, res$target))
print(cells.of.interest)

cell_color = c("magenta2", "green3")
names(cell_color) <- c('4', '6')
cell_color = cell_color[match(cells.of.interest, names(cell_color))]

# Mac_BL_early #2AD0B7
# Neu_BL_early #98B304
# Epidermis_BL_early #16005e
# CT_BL_early_1 #005e45
# CT_BL_early_3 #2f7c67
# Epidermis_BL.CSD_early #3200F5

pdfname = paste0(outDir, '/LR_interactions_LIANA_tops.pdf')
pdf(pdfname, width=12, height = 8)
for(ntop in c(100, 200, 300))
{
  # ntop = 100
  cat('top LR -- ', ntop, '\n')
  test = res[c(1:ntop), ]
  
  # jj = which(test$ligand == 'RGMB'|test$receptor == 'RGMB'|
  #              test$ligand == "FGFR3")
  # if(length(jj) >0){
  #   test = test[-jj, ]
  # }
  
  #test = test[-which(test$ligand == 'SPON1'), ] 
  
  my_CircosPlot(test, 
                weight.attribute = 'weight_norm',
                cols.use = cell_color,
                sources.include = cells.of.interest,
                targets.include = cells.of.interest,
                lab.cex = 0.5,
                title = paste('LR scores top :', ntop))
  
}

dev.off()

res = res[which(res$source != res$target), ]

pdfname = paste0(outDir, '/LR_interactions_LIANA_tops_noAutoregulation.pdf')
pdf(pdfname, width=12, height = 8)
for(ntop in c(100, 200))
{
  # ntop = 100
  cat('top LR -- ', ntop, '\n')
  if(ntop > nrow(res)) ntop = nrow(res)
  test = res[c(1:ntop), ]
  
  # jj = which(test$ligand == 'RGMB'|test$receptor == 'RGMB'|
  #              test$ligand == "FGFR3")
  # if(length(jj) >0){
  #   test = test[-jj, ]
  # }
  
  #test = test[-which(test$ligand == 'SPON1'), ] 
  
  my_CircosPlot(test, 
                weight.attribute = 'weight_norm',
                cols.use = cell_color,
                sources.include = cells.of.interest,
                targets.include = cells.of.interest,
                lab.cex = 0.5,
                title = paste('LR scores top :', ntop))
  
}

dev.off()

FeaturePlot(aa, features = c('Bmp7', 'Bmp4', 'Bmp1'))

ggsave(filename = paste0(outDir, '/FeaturePlot_Bmps.pdf'), 
       width = 12, height = 8)

FeaturePlot(aa, features = c('Wnt4', 'Wnt1', 'Wnt5b', 'Fzd4', 'Fzd2', 'Fzd7'))

ggsave(filename = paste0(outDir, '/FeaturePlot_Wnt.pdf'), 
       width = 12, height = 12)
