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

version.analysis = '_R13547_10x_mNT_20261118'
species = 'mNT_scRNAseq'

resDir = paste0("../results/revision_figures", version.analysis)

figureDir = '/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/plots_manuscript_1/'
RdataDir = '../results/Rdata/'
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')
source(paste0(functionDir, '/functions_dataIntegration.R'))


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



########################################################
########################################################
# Section I : compare the scRNA-seq embryo data with the mNTO
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam_clusters_4save.rds'))

levels_sels = c("day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
cols_sel = cols[match(levels_sels, names(cols))]


DimPlot(aa, label = FALSE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE) +
  NoLegend()


DimPlot(aa, label = TRUE, group.by =  'clusters', repel = TRUE, raster=FALSE, 
        cols = cols_cluster) + NoLegend()

aa$celltypes[which(aa$clusters == 2)] = 'Head'
aa$celltypes[which(aa$clusters == 3)] = 'pFP'
aa$celltypes[which(aa$clusters == 5)] = 'eNP'

DimPlot(aa, label = TRUE, group.by =  'celltypes', repel = TRUE, raster=FALSE) + NoLegend()

knownGenes =  c('Zfp42', 'Utf1', 'Pou5f1',  
                'Cyp26a1', 'Stra8','Hoxa1', 'Hoxa3', 'Hoxb4', 'Sox2',
                'Pax6', 'Sox1', 'Pax3',   'Irx3', 'Irx5',
                'Foxa2', 'Shh', 'Arx', "Spon1", 'Slit2', "Ntn1",
                'Nkx6-1', 'Nkx2-2',  'Olig2',  
                'Tubb3',  'Neurod4', 'Elavl4', 'Pou3f2')
genes.to.plot <- knownGenes

library(RColorBrewer)
feat_cols = c("#F0F0F0", "#EFFAB6", "#69C6BE", "#007BB7", "#121D60")
#to set the gradient of the colour scale from min to max
feat_values = c(-1, -0.5, 0.5, 1.5)
#colGEX = c("grey85", brewer.pal(7, "Reds"))
aa$clusters = factor(aa$clusters, levels = c(7, 5, 6, 4, 3, 2, 1))

DotPlot(aa, group.by = "clusters", features = genes.to.plot) + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  #scale_colour_viridis(option="magma") +
  #guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  scale_color_gradientn(colors = feat_cols) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(axis.text.x = element_text(angle = 90, size = 12, vjust = 0.),
        axis.text.y = element_text(angle = 0, size = 14)) +
  labs( x = '', y = 'Clusters' )


##########################################
# prepare the mouse gastrulation data  
##########################################
bb = readRDS(file = paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/scRNAseq_mNT/saved_seuratObj/',
                           'seuratObject_EmbryoAtlasData_all36sample.rds'))

celltype_colours = c("Epiblast" = "#635547",
                     "Primitive Streak" = "#DABE99",
                     "Caudal epiblast" = "#9e6762",
                     
                     "PGC" = "#FACB12",
                     
                     "Anterior Primitive Streak" = "#c19f70",
                     "Notochord" = "#0F4A9C",
                     "Def. endoderm" = "#F397C0",
                     "Gut" = "#EF5A9D",
                     
                     "Nascent mesoderm" = "#C594BF",
                     "Mixed mesoderm" = "#DFCDE4",
                     "Intermediate mesoderm" = "#139992",
                     "Caudal Mesoderm" = "#3F84AA",
                     "Paraxial mesoderm" = "#8DB5CE",
                     "Somitic mesoderm" = "#005579",
                     "Pharyngeal mesoderm" = "#C9EBFB",
                     "Cardiomyocytes" = "#B51D8D",
                     "Allantois" = "#532C8A",
                     "ExE mesoderm" = "#8870ad",
                     "Mesenchyme" = "#cc7818",
                     
                     "Haematoendothelial progenitors" = "#FBBE92",
                     "Endothelium" = "#ff891c",
                     "Blood progenitors 1" = "#f9decf",
                     "Blood progenitors 2" = "#c9a997",
                     "Erythroid1" = "#C72228",
                     "Erythroid2" = "#f79083",
                     "Erythroid3" = "#EF4E22",
                     
                     "NMP" = "#8EC792",
                     
                     "Rostral neurectoderm" = "#65A83E",
                     "Caudal neurectoderm" = "#354E23",
                     "Neural crest" = "#C3C388",
                     "Forebrain/Midbrain/Hindbrain" = "#647a4f",
                     "Spinal cord" = "#CDE088",
                     
                     "Surface ectoderm" = "#f7f79e",
                     
                     "Visceral endoderm" = "#F6BFCB",
                     "ExE endoderm" = "#7F6874",
                     "ExE ectoderm" = "#989898",
                     "Parietal endoderm" = "#1A1A1A"
                     
)

stage_colours = c("E6.5" = "#D53E4F",
                  "E6.75" = "#F46D43",
                  "E7.0" = "#FDAE61",
                  "E7.25" = "#FEE08B",
                  "E7.5" = "#FFFFBF",
                  "E7.75" = "#E6F598",
                  "E8.0" = "#ABDDA4",
                  "E8.25" = "#66C2A5",
                  "E8.5" = "#3288BD"
)

coexpcol = c("Not" = "grey", "CoExp" = "purple")

## subset cell types
Idents(bb) = as.factor(bb$celltype)

VlnPlot(bb, features = c('Foxa2', 'Pax6'), group.by = 'celltype', ncol = 1)

ggsave(filename = paste0(resDir, '/mouse_gastrulation_VlnPlot_Foxa2_Pax6.pdf'), 
       width = 12, height = 8)

## select the follow cell types specified by Hannah
'Spinal cord'
'Forebrain/Midbrain/Hindbrain'
'Rostral neurectoderm'
'Caudal neurectoderm'
"Forebrain/Midbrain/Hindbrain"

## select the cell types expressed either Pax6 or FoxA2
celltypes_sels = c('Anterior Primitive Streak',
                   
                   'Def. endoderm', 
                   'Epiblast', 
                   'Forebrain/Midbrain/Hindbrain',
                   'Gut', 
                   'NMP', 
                   'Notochord',
                   'Primitive Streak',
                   'Rostral neurectoderm',
                   'Spinal cord', 
                   'Surface ectoderm',
                   'Visceral endoderm',
                   
                   'Forebrain/Midbrain/Hindbrain',
                   'Rostral neurectoderm',
                   'Caudal neurectoderm')

celltypes_sels = unique(celltypes_sels)

sub.obj = subset(bb, idents = celltypes_sels)

VlnPlot(sub.obj, features = c('Foxa2', 'Pax6'), group.by = 'celltype', ncol = 1)

ggsave(filename = paste0(resDir, '/mouse_gastrulation_selectedCellTypes_VlnPlot_Foxa2_Pax6.pdf'), 
       width = 12, height = 8)

FeatureScatter(sub.obj, feature1 = 'Foxa2', feature2 = 'Pax6', 
               group.by = 'stage', cols = stage_colours, raster = FALSE, pt.size = 1.5) +
  theme(plot.title = element_blank()) + 
  labs( x = 'FoxA2', y = 'Pax6' ) +
  theme(axis.text.x = element_text(angle = 0, size = 12),
        axis.text.y = element_text(angle = 0, size = 12))

## or subset the FoxA2+ & Pax6+ cells
xx = subset(x = sub.obj, subset = Foxa2 > 0 & Pax6 > 0)

counts = table(xx$celltype)

pcts = counts
for(n in 1:length(pcts))
{
  pcts[n] = pcts[n]/length(which(sub.obj$celltype == names(pcts)[n]))
}

pcts = pcts[order(-pcts)]

# sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
# sub.obj = ScaleData(sub.obj)
# sub.obj <- RunPCA(sub.obj, verbose = FALSE, weight.by.var = TRUE)
# ElbowPlot(sub.obj, ndims = 50)

#Idents(sub.obj) = sub.obj$condition

sub.obj <- RunUMAP(sub.obj, reduction = 'pca.corrected', dims = 1:30, n.neighbors = 50, min.dist = 0.1)

p1 = DimPlot(sub.obj, group.by = 'celltype', label = TRUE, repel = TRUE, raster = FALSE)

p2 = FeaturePlot(sub.obj, features = c('Foxa2', 'Pax6'))

p1 / p2

ggsave(filename = paste0(resDir, '/mouse_gastrulation_selectedCellTypes_umap_Foxa2_Pax6.pdf'), 
       width = 12, height = 12)

saveRDS(sub.obj, file = paste0(RdataDir, 'seuratObject_EmbryoAtlasData_subset_celltypes_FoxA2.Pax6.rds'))


##########################################
# test the umap with  only TFs and signaling molecules
##########################################
tfs_sps = unique(c(tfs, sps))

mm = intersect(which(!is.na(match(tfs_sps, rownames(aa)))), which(!is.na(match(tfs_sps, rownames(sub.obj)))))

xx = ScaleData(sub.obj, features = tfs_sps[mm])
xx[['umap.correctd']] = xx[['umap']]

xx <- RunPCA(xx, verbose = FALSE,  features = rownames(xx), weight.by.var = TRUE)
ElbowPlot(xx, ndims = 50)
xx <- RunUMAP(xx, reduction = 'pca', dims = 1:30, n.neighbors = 50, min.dist = 0.1)

p1 = DimPlot(xx, group.by = 'celltype', label = TRUE, repel = TRUE, raster = FALSE)

p2 = FeaturePlot(xx, features = c('Foxa2', 'Pax6'))

p1 / p2


ggsave(filename = paste0(resDir, '/mouse_gastrulation_selectedCellTypes_tfs_umap_Foxa2_Pax6.pdf'), 
       width = 12, height = 12)


saveRDS(xx, file = paste0(RdataDir, 'seuratObject_EmbryoAtlasData_subset_celltypes_FoxA2.Pax6_umap_tfs.rds'))

##########################################
# processe the mNTO and embryo data  using only TFs and SPs
##########################################
ref = readRDS(file = paste0(RdataDir, 
                           'seuratObject_EmbryoAtlasData_subset_celltypes_FoxA2.Pax6_umap_tfs.rds'))

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam_clusters_4save.rds'))

aa$celltypes[which(aa$clusters == 2)] = 'Head'
aa$celltypes[which(aa$clusters == 3)] = 'pFP'
aa$celltypes[which(aa$clusters == 5)] = 'eNP'
aa$celltypes[which(aa$clusters == 1)] = 'c1'
aa$celltypes[which(aa$clusters == 7)] = 'NP'
aa$celltypes[which(aa$clusters == 6)] = 'FP'
aa$celltypes[which(aa$clusters == 4)] = 'c5'

aa$dataset = 'mNT'
aa$stage = aa$condition
aa$celltype = aa$celltypes
aa$celltype = paste0('mNT_', aa$celltype)

ref$dataset = 'ref'

tfs_sps = unique(c(tfs, sps))
mm = intersect(which(!is.na(match(tfs_sps, rownames(aa)))), 
               which(!is.na(match(tfs_sps, rownames(ref)))))
xx = subset(ref, features = tfs_sps[mm])
yy = subset(aa, features = tfs_sps[mm])

#DimPlot(yy, group.by = 'clusters', label = TRUE)
DimPlot(yy, group.by = 'celltype', label = TRUE)
DimPlot(xx, group.by = 'celltype', label = TRUE)

newxx = CreateSeuratObject(counts = xx@assays$originalexp@counts, assay = 'RNA', 
                           meta.data = xx@meta.data)

refs.merged = merge(yy, y = newxx, add.cell.ids = c("mNT", "mouseGastrulation"), 
                    project = "RA_competence", merge.data = TRUE)

saveRDS(refs.merged, file = paste0(RdataDir, 
                                   'mergedSeuratObject_EmbryoAtlasData_subsetCelltypes_mNTO_tfs.rds'))


##########################################
# calculate the similarity between the mNTO cell types and those mouse embryo 
##########################################
refs.merged = readRDS(file =  paste0(RdataDir, 
                                   'mergedSeuratObject_EmbryoAtlasData_subsetCelltypes_mNTO_tfs.rds'))

ref.list <- SplitObject(refs.merged, split.by = "dataset")

ref.list <- lapply(X = ref.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize")
  x = FindVariableFeatures(x, nfeatures = 500)
  #x <- ScaleData(x, features = features.common, verbose = TRUE)
  #x <- RunPCA(x, features = features.common, verbose = FALSE)
})

#tfs_sel <- intersect(VariableFeatures(ref.list[[1]]), VariableFeatures(ref.list[[2]]))

tfs_sel = VariableFeatures(ref.list[[2]])
cat(length(tfs_sel), ' tfs and sps selected for correlation calculation \n')


#an_orthologs$tfs = FALSE
#an_orthologs$tfs[which(!is.na(match(toupper(an_orthologs$query), c(tfs, sps))))] = TRUE
avg.ax = AggregateExpression(refs.merged, features = tfs_sel, slot = 'data', group.by = "celltype")
avg.ax = data.frame(avg.ax$RNA)

cat(nrow(avg.ax), 'TFs used for correlation analysis\n')

avg.nt = avg.ax[, grep('mNT_', colnames(avg.ax))]  
avg.me = avg.ax[, grep('mNT_', colnames(avg.ax), invert = TRUE)]  

cort = psych::corr.test(avg.nt, avg.me, method = "spearman", 
                        adjust = "fdr", alpha = 0.05, ci = F)

range(cort$r)

cort$maxrow = apply(cort$r, 1, which.max)
cort$maxcol = apply(cort$r, 2, which.max)

# cluster and order labels
hcr = hclust(dist(cort$r), method = "ward.D2")
hcc = hclust(dist(t(cort$r)), method = "ward.D2")
hcr = hcr$labels[hcr$order]
hcc = hcc$labels[hcc$order]

# reshaping the correlations
plot_df = reshape2::melt(cort$r)
plot_df$Var1 = factor(plot_df$Var1, levels = rev(hcr))
plot_df$Var2 = factor(plot_df$Var2, levels = hcc)

# add pvalue and max cor infor
plot_df$padj = -log10(reshape2::melt(cort$p.adj+min(cort$p.adj[cort$p.adj>0])/10)$value)
plot_df$rowmax = apply(Reduce(cbind, lapply(names(cort$maxrow), 
                                            function(n) plot_df$Var1==n &
                                              plot_df$Var2==colnames(cort$r)[cort$maxrow[n]])), 
                       1, any)
plot_df$colmax = apply(Reduce(cbind, lapply(names(cort$maxcol), 
                                            function(n) plot_df$Var2==n &
                                              plot_df$Var1==rownames(cort$r)[cort$maxcol[n]])), 
                       1, any)
plot_df$markcol = plot_df$value>quantile(plot_df$value, 0.98)


# getting a colourscale where 0 is white in the middle, and intensity leveled by max(abs(value))
cols = colorRampPalette(c(rev(RColorBrewer::brewer.pal(9, "Blues")),
                          RColorBrewer::brewer.pal(9, "Reds")))(101)
#cols = colorRampPalette(c(RColorBrewer::brewer.pal(9, "Reds")))(101)
br = seq((min(cort$r)), max(abs(cort$r)), length.out = 101)
cols = cols[!(br>max(cort$r) | br<min(cort$r))]

levels_y = c('mNT_c1', 'mNT_Head', 'mNT_pFP', 
             'mNT_eNP', 'mNT_c5', 'mNT_FP', 'mNT_NP')
levels_y = levels_y[c(length(levels_y):1)]

ggplot() +
  geom_point(data = plot_df, mapping = aes(x = Var2, 
                                           y = factor(Var1, levels = levels_y), 
                                           fill = value), 
             shape = 21, size = 10) +
  #geom_point(data = plot_df[plot_df$rowmax,], mapping = aes(x = Var2, y = Var1, size = 30), 
  #           shape = "—", show.legend = F, colour = "grey10")+
  #geom_point(data = plot_df[plot_df$colmax,], mapping = aes(x = Var2, y = Var1, size = 30), 
  #           shape = "|", show.legend = F, colour = "grey10")+
  scale_x_discrete(expand = c(0,0.7)) +
  scale_y_discrete(expand = c(0,0.7)) +
  scale_fill_gradientn(breaks = signif(c(min(cort$r)+0.005, 0, max(cort$r)-0.005),2), 
                       #values = scales::rescale(c(min(br), 0, max(br))),
                       colours = cols) +
  labs(x = 'mouse', y = 'axolotl', fill = "Spearman's\nrho", size = "-log10\nadj. p-value")+
  theme_classic()+
  theme(axis.title = element_text(colour = "black", face = "bold", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1.0, vjust = 1.0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))

ggsave(paste0(resDir, "mNTO_mEmbryo_correlationAnalysis.pdf"), width=8, height=6, dpi=300)


##########################################
# try reference mapping
##########################################
Test_Seurat_query_ref_mapping = FALSE 
if(Test_Seurat_query_ref_mapping){
  mapping_method = "seurat_query_ref_mapping"
  
  ref = readRDS(file = paste0(RdataDir, 
                              'seuratObject_EmbryoAtlasData_subset_celltypes_FoxA2.Pax6_umap_tfs.rds'))
  
  aa = readRDS(file = paste0(RdataDir, 
                             'seuratObj_clustersFiltered_umap_RAsamples_selectUMAPparam_clusters_4save.rds'))
  
  aa$celltypes[which(aa$clusters == 2)] = 'Head'
  aa$celltypes[which(aa$clusters == 3)] = 'pFP'
  aa$celltypes[which(aa$clusters == 5)] = 'eNP'
  aa$celltypes[which(aa$clusters == 1)] = 'c1'
  aa$celltypes[which(aa$clusters == 7)] = 'NP'
  aa$celltypes[which(aa$clusters == 6)] = 'FP'
  aa$celltypes[which(aa$clusters == 4)] = 'c5'
  
  aa$dataset = 'mNT'
  aa$stage = aa$condition
  aa$celltype = aa$celltypes
  aa$celltype = paste0('mNT_', aa$celltype)
  
  ref$dataset = 'ref'
  
  tfs_sps = unique(c(tfs, sps))
  mm = intersect(which(!is.na(match(tfs_sps, rownames(aa)))), 
                 which(!is.na(match(tfs_sps, rownames(ref)))))
  
  ref = subset(ref, features = tfs_sps[mm])
  aa = subset(aa, features = tfs_sps[mm])
  
  
  #ElbowPlot(ref, ndims = 50, reduction = 'pca')
  #ref = RunUMAP(ref, reduction = "pca", dims = 1:30, n.neighbors = 30, 
  #              min.dist = 0.1, return.model = TRUE) 
  
  DimPlot(ref, reduction = "umap", 
          group.by = "celltype", label = TRUE,
          repel = TRUE, raster=FALSE) 
  
  ref$labels = factor(paste0(ref$celltype, '_', ref$stage))
  
  # In data transfer, Seurat has an option (set by default) to project the PCA structure of a reference 
  # onto the query, instead of learning a joint structure with CCA. 
  # We generally suggest using this option when projecting data between scRNA-seq datasets.
  anchors <- FindTransferAnchors(reference = ref, 
                                 query = aa, 
                                 dims = 1:50,
                                 normalization.method = "LogNormalize",
                                 reference.assay = 'originalexp',
                                 query.assay = 'RNA', 
                                 features = rownames(aa),
                                 #reference.reduction = "pca.corrected",
                                 max.features = 200,
                                 k.anchor = 10,
                                 reduction = "cca" 
  )
  saveRDS(anchors, file = paste0(RdataDir, '/cca_anchors.rds'))
  
  #pancreas.anchors <- FindTransferAnchors(reference = pancreas.ref, query = pancreas.query, dims = 1:30,
  #                                        reference.reduction = "pca")
  anchors = readRDS( file = paste0(RdataDir, '/cca_anchors.rds'))
  predictions <- TransferData(anchorset = anchors, refdata = ref$labels, 
                              reference = ref,
                              query = aa,
                              weight.reduction = 'cca',
                              dims = 1:30)
  
  
  
  query = aa;
  query$predicted.id = predictions$predicted.id 
  query$predicted.score = predictions$predicted.id.score
  
  
  p1 = DimPlot(query, reduction = "umap", 
               group.by = "predicted.id", label = TRUE,
               repel = TRUE, raster=FALSE) 
  p2 = FeaturePlot(query, features = 'predicted.score')
  
  p1
  ggsave(paste0(outDir, '/transfer_learning_seurat_cca.ref_v3.pdf'), 
         width = 16, height = 8)
  
  
  
  ## visualize the query cells alongside our reference and didn't work well
  ## (MapQuery in Seurat didn't work well)
  # query <- Seurat::MapQuery(anchorset = anchors, 
  #                   reference = ref, 
  #                   query = aa,
  #                   refdata = list(labels = "labels"),
  #                   #refdata = list(celltype = "celltype"), 
  #                   transferdata.args = list(weight.reduction = 'rpca.ref'), 
  #                   reference.reduction = "pca", 
  #                   reduction.model = "umap")
  # 
  # p1 <- DimPlot(ref, reduction = "umap", group.by = "celltype", 
  #               label = TRUE, label.size = 3,
  #               repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
  # p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.labels", 
  #               label = TRUE,
  #               label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
  # p1 + p2
  
  
}



ref.list <- SplitObject(refs.merged, split.by = "dataset")

#ref.list <- SplitObject(refs.merged, split.by = "dataset")
#rm(list = c('refs.merged')) # remove big seurat objects to clear memory

# normalize and identify variable features for each dataset independently
features.common = rownames(refs.merged)

ref.list <- lapply(X = ref.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize")
  x <- ScaleData(x, features = features.common, verbose = TRUE)
  x <- RunPCA(x, features = features.common, verbose = FALSE)
  
})


ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                      anchor.features = features.common, 
                                      reference = c(2),
                                      reduction = "cca", 
                                      #reduction = 'rpca',
                                      k.anchor = 5,
                                      dims = 1:50)

rm(ref.list)

# this command creates an 'integrated' data assay
ref.combined <- IntegrateData(anchorset = ref.anchors, features.to.integrate = features.common, 
                              dims = 1:50) ## take ~80G memory

rm(ref.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(ref.combined) <- "integrated"

ref.combined <- ScaleData(ref.combined, verbose = FALSE)
ref.combined <- RunPCA(ref.combined, npcs = 50, verbose = FALSE)

ElbowPlot(ref.combined, ndims = 50)

#kk = which(ref.combined$dataset == 'mNT') 
#ref.combined$celltype[kk] = paste0('mNT_', ref.combined$condition[kk])
#names(cols_sel) = paste0('mNT_', names(cols_sel))

#cols_used = c(cols_mouse, cols_sel)
#ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
#ref.combined <- FindClusters(ref.combined, resolution = 0.2)

ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:50, n.neighbors = 50, 
                        min.dist = 0.2) 

DimPlot(ref.combined, reduction = "umap", group.by = "celltypes", label = TRUE,
        repel = TRUE, raster=FALSE)

saveRDS(ref.combined, file = paste0(RdataDir, '/integrated_mNT_mouseGastrulation_SeuratCCA.rds'))

ref.combined = readRDS(file = paste0(RdataDir, '/integrated_mNT_mouseGastrulation_SeuratCCA.rds'))


saveRDS(cols_used, file = paste0(outDir, '/integrated_mNT_mouseGastrulation_colorsUsed.rds'))

# Visualization
DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
        repel = TRUE, raster=FALSE, cols = c(cols_mouse, cols_sel)) 

ggsave(paste0(outDir, '/Integration_mNT_celltypes_noref.batchcorrection.pdf'), 
       width = 16, height = 8)

#DimPlot(ref.combined, reduction = "umap")



