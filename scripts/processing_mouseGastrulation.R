##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: process the mouse gastrulation data from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3457441
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Oct 30 14:55:46 2023
##########################################################################
##########################################################################
rm(list = ls())

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')
source(paste0(functionDir, 'functions_dataIntegration.R'))

library(pryr) # monitor the memory usage
require(ggplot2)

require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(DropletUtils)
library(edgeR)
library(future)
options(future.globals.maxSize = 160000 * 1024^2)
mem_used()

species = 'mm10'
version.analysis = '_mouse_gastrulation_Chan.et.al'
resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../mouse_gastrulation/Chan_et_al/'

########################################################
########################################################
# Section I : process the chan et al. gastrulation data
# 
########################################################
########################################################

##########################################
# import the raw data 
##########################################
dir_list = list.dirs(path = dataDir, full.names = TRUE)
dir_list = dir_list[grep('GSE122187_WT', dir_list)]
conds = basename(dir_list)

design = data.frame(conds, gsub('GSE122187_WT_', '', conds))
colnames(design) = c('sampleID', 'condition')

rm(list = c("conds", "dir_list"))

# import data from cellranger output
for(n in 1:nrow(design))
{
  # n = 1
  cat(n, ' : ', design$condition[n], '\n')
  
  topdir = paste0(dataDir, design$sampleID[n])
  exp = Matrix::readMM(paste0(topdir, "/matrix.mtx")) #read matrix
  bc = read.csv(paste0(topdir, "/barcodes.tsv"), header = F, stringsAsFactors = F)
  g = read.csv(paste0(topdir, "/genes.tsv"), header = F, stringsAsFactors = F, sep = '\t')
  
  ## make unique gene names
  g$name = g$V2
  gg.counts = table(g$V2)
  gg.dup = names(gg.counts)[which(gg.counts>1)]
  index.dup = which(!is.na(match(g$V2, gg.dup)))
  g$name[index.dup] = paste0(g$V2[index.dup], '_', g$V1[index.dup])
  
  colnames(exp) = bc$V1
  rownames(exp) = g$name
  
  count.data = exp
  rm(exp);
  
  cat('get empty drops with UMI rank \n')
  
  # get emptyDrops and default cutoff cell estimates
  #iscell_dd = defaultDrops(count.data, expected = 8000) # default cell estimate, similar to 10x cellranger
  #sum(iscell_dd, na.rm=TRUE)
  
  ## not used the emptyDrops too slow 
  # eout = emptyDrops(count.data, lower = 200)
  # eout$FDR[is.na(eout$FDR)] = 1
  # iscell_ed = eout$FDR<=0.01
  # sum(iscell_ed, na.rm=TRUE)
  meta = data.frame(row.names = colnames(count.data), 
                    cellId = colnames(count.data),
                    condition = design$condition[n])
  
  # plot rankings for number of UMI
  br.out <- barcodeRanks(count.data)
  
  # pdf(paste0(resDir, "/UMIrank_emptyDrop_", design$condition[n], "_", design$sampleID[n],  ".pdf"), 
  #     height = 6, width =10, useDingbats = FALSE)
  # 
  # plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  # 
  # o <- order(br.out$rank)
  # lines(br.out$rank[o], br.out$fitted[o], col="red")
  # abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  # abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  # abline(v = sum(iscell_dd), col = 'darkgreen', lwd = 2.0)
  # abline(v = c(3000, 5000, 8000, 10000, 12000), col = 'gray')
  # text(x = c(3000, 5000, 8000, 10000, 12000), y =10000, labels = c(3000, 5000, 8000, 10000, 12000), 
  #      col = 'red')
  # legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
  #        legend=c("knee", "inflection"))
  # 
  # dev.off()
  
  # use defaultDrop to select cells.
  aa = CreateSeuratObject(counts = count.data,
                          meta.data = meta, 
                          min.cells = 20, min.features = 100)
  
  aa$cell.id = paste0(colnames(aa), '_', design$condition[n], '_', design$sampleID[n])
  
  if(n == 1) {
    mnt = aa
  }else{
    mnt = merge(mnt, aa)
  }
  
  rm(aa)
  
}

mnt[["percent.mt"]] <- PercentageFeatureSet(mnt, pattern = "^mt-")

save(design, mnt, 
     file = paste0(RdataDir, 'seuratObject_design_variableGenes', version.analysis, '.Rdata'))

rm(count.data); rm(br.out); rm(bc);rm(meta)

##########################################
# QCs and cell filtering  
##########################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes', version.analysis, '.Rdata'))
aa = mnt
rm(mnt)

levels = design$condition

pdfname = paste0(resDir, '/QCs_nCounts_nFeatures_percentMT.pdf')
pdf(pdfname, width=16, height = 8)

table(aa$condition) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(condition = factor(Var1)) %>%
  #mutate(cellNbs = integer(Freq))
  ggplot(aes(x=condition, y=Freq)) +
  geom_bar(stat="identity", width=0.5) +
  theme_classic() +
  labs( x = '', y = 'detected cell # from cellRanger barcodes' )  +
  theme(axis.text.x = element_text(angle = 90, size = 10)) + 
  geom_hline(yintercept = c(3000, 5000, 7000), col = 'red')

Idents(aa) = factor(aa$condition, levels = levels)

VlnPlot(aa, features = 'nFeature_RNA', y.max = 10000)
VlnPlot(aa, features = 'nCount_RNA', y.max = 100000)
VlnPlot(aa, features = 'percent.mt', y.max =50)

FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "percent.mt")

#ggs = sapply(rownames(aa), function(x) {x = unlist(strsplit(x, '-')); x = x[grep('ENSMUSG', x, invert = TRUE)];
#                                                                               paste0(x, collapse = '-')})

Check_geneExamples = FALSE
if(Check_geneExamples){
  features = c('Pou5f1', 'Sox2', 'Lef1', 'Otx2', 'Zfp703', 'Pax6', 'Foxa2', 'Shh', 'Nkx6-1', 'Nkx2-2', 'Olig2', 
               'Sox1', 'Tubb3')
  features = rownames(aa)[!is.na(match(rownames(aa), features))]
  
  Idents(aa) = factor(aa$condition, levels = levels)
  
  for(n in 1:length(features))
  {
    p1 = VlnPlot(aa, features = features[n])
    plot(p1)
  }
  
}

dev.off()


##########################################
## filter cells here
##########################################
VlnPlot(aa, features = 'nFeature_RNA', y.max = 10000) +
  geom_hline(yintercept = c(1000, 8000), col = 'red')

VlnPlot(aa, features = 'percent.mt', y.max =20) + 
  geom_hline(yintercept = c(5), col = 'red')

aa <- subset(aa, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_', species, version.analysis, '.rds'))


##########################################
# normalization, umap and clustering
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_', species, version.analysis, '.rds'))

Idents(aa) = factor(aa$condition, levels = levels)

aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)

aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)

ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.3)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 1.0)

DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)

ggsave(filename = paste0(resDir, '/first_test_umap.pdf'), width = 10, height = 8)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)

p1 + p2

ggsave(filename = paste0(resDir, '/umap_timepoints_clusters.pdf'), width = 16, height = 6)


features = c('Pou5f1', 'Sox2', 'Lef1', 'Otx2', 'Zfp703', 'Pax6', 'Foxa2', 
             'Shh', 'Nkx6-1', 'Nkx2-2', 'Olig2', 'Sox1', 'Tubb3',
             'Otx2', 'Zic2', 'Nodal', 'Fgf5', 'Fgf4', 'Utf1'
             )

features = rownames(aa)[!is.na(match(rownames(aa), features))]
Idents(aa) = factor(aa$condition, levels = levels)

pdfname = paste0(resDir, '/Gene_Examples.pdf')
pdf(pdfname, width=16, height = 8)

for(n in 1:length(features))
{
  cat(n, ' -- ',  features[n], '\n')
  p1 = FeaturePlot(aa, features = features[n])
  plot(p1)
}

dev.off()

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                          '_lognormamlized_pca_umap_clustered.rds'))


##########################################
# reprocess the data according to the original paper
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                           '_lognormamlized_pca_umap_clustered.rds'))

VlnPlot(aa, features = 'nCount_RNA', y.max = 100000) +
  geom_hline(yintercept = c(1000, 8000), col = 'red')

# reproduce the processing from the original paper
aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
aa <- FindVariableFeatures(aa, selection.method = "dispersion", 
                           mean.function = ExpMean,
                           dispersion.function = LogVMR, 
                           mean.cutoff = c(0.0125, 3),
                           dispersion.cutoff = c(0.5, Inf)
                           )

aa <- ScaleData(aa, vars.to.regress = 'nCount_RNA')
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)

ElbowPlot(aa, ndims = 50)

aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 1.0)

Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 50, min.dist = 0.3)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)

p1 + p2

ggsave(filename = paste0(resDir, '/umap_timepoints_clusters_regress.nCounts.pdf'), width = 16, height = 6)


saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                          '_lognormamlized_var.to.regress.nCount.RNA_pca_clustering_umap.rds'))


##########################################
# assign cell type lables according to the correlation of marker genes defined by the paper 
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                           '_lognormamlized_var.to.regress.nCount.RNA_pca_clustering_umap.rds'))

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)

p1 + p2

## over clustering to assign cell labels with correlations
aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 10)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
ggsave(filename = paste0(resDir, '/umap_overclustered_161clusters.pdf'), width = 12, height = 6)

# pseudo-bulk by per donor per cell type
pb <- AverageExpression(aa, return.seurat = TRUE, slot = 'counts', group.by = c("seurat_clusters"))
clusters = pb@assays$RNA@data

## import the vectors of cell states used by original paper
refs = readxl::read_xls(paste0(dataDir, 'GSE122187_CellStateKernels.xls'), col_names = TRUE, skip = 1)
refs = data.frame(refs)
rownames(refs) = refs[,1]
refs = refs[, -1]

mm = match(rownames(refs), rownames(clusters))
refs = refs[which(!is.na(mm)), ]
clusters = clusters[mm[which(!is.na(mm))], ]

keep = matrix(NA, nrow = ncol(clusters), ncol = 6)
rownames(keep) = colnames(clusters)
colnames(keep) = c('cluster.pearson', 'cluster.spearman', 'dist', 
                   'cor.pearson', 'cor.spearman', 'cluster.dist')

library(lsa)
euclidean <- function(a, b) sqrt(sum((a - b)^2))

for(n in 1:ncol(clusters))
{
  # n = 1
  cat(n, ' -- ', colnames(clusters)[n], '\n')
  cors = cor(clusters[,n], refs, use = 'everything', method = 'pearson')
  keep[n, 1] = as.numeric(gsub('X', '', colnames(refs)[which.max(cors)]))
  keep[n, 4] = cors[which.max(cors)]
  
  cors = cor(clusters[,n], refs, use = 'everything', method = 'spearman')
  keep[n, 2] = as.numeric(gsub('X', '', colnames(refs)[which.max(cors)]))
  keep[n, 5] = cors[which.max(cors)]
  
  cs = rep(NA, ncol(refs))
  for(m in 1:ncol(refs)){
    #cs[m] = cosine(clusters[,n],  refs[,m])
    cs[m] = euclidean(clusters[,n],  refs[,m])
  }
  
  keep[n, 3] = as.numeric(gsub('X', '', colnames(refs)[which.min(cs)]))
  keep[n, 6] = cs[which.min(cs)]
  
}

keep = data.frame(keep)
keep$cluster = NA
aa$clusterID = NA

for(n in 1:nrow(keep))
{
  # n = 1
  test = table(as.numeric(keep[n, 1:3]))
  kk = which(test>1)
  if(length(kk)>0){
    id = names(test[kk])
    cat(n, ' -- ', id, '\n' )
    keep$cluster[n] = id
    aa$clusterID[which(aa$seurat_clusters == rownames(keep)[n])] = id
  }else{
    cat(n, ' -- no id found \n' )
  }
  
}

annots = openxlsx::read.xlsx(paste0(dataDir, 'clusterID_cellstates.xlsx'), colNames = FALSE, rowNames = FALSE)
annots$clusterID = annots[,1]
annots$celltype = annots[,1]
annots$clusterID = sapply(annots$clusterID, function(x){
  unlist(strsplit(x, ' - '))[1]
})
annots$celltype = sapply(annots$celltype, function(x){
  unlist(strsplit(x, ' - '))[2]
})

annots$celltype = gsub(' $','', annots$celltype)
annots$celltype = gsub(' ','_', annots$celltype)

aa$celltype = NA
for(n in 1:nrow(annots))
{
  aa$celltype[which(aa$clusterID == annots$clusterID[n])] = annots$celltype[n]
}

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE) + NoLegend()

p2

ggsave(filename = paste0(resDir, '/umap_tranferred_celltypes_scmap.method.pdf'), width = 12, height = 10)


saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis,
                          '_lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltypes.rds'))


Test_batchCorrection_fastMNN = FALSE
if(Test_batchCorrection_fastMNN){
  
  aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis,
                '_lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltypes.rds'))
  
  aa$condition = factor(aa$condition)
 
  source(paste0(functionDir, 'functions_dataIntegration.R'))
  xx = IntegrateData_runFastMNN(aa, group.by = 'condition', nfeatures = 3000,
                                ndims = c(1:50), 
                                merge.order = c("E8.5_1ab", "E8.0_1ab","E7.5_2", "E7.5_1", "E7.0_1", "E6.5_1"),
                                correct.all = TRUE,
                                reference = NULL)
  
  xx[['pca']] = aa[['pca']]
  xx[['umap.mnn']] = xx[['umap']]
  xx[['umap']] = aa[['umap']]
  
  aa = xx
  rm(xx)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE) + NoLegend()
  
  p1 / p2  
  
  ggsave(filename = paste0(resDir, '/umap_tranferred_celltypes_fastMNN.batchCorrected.pdf'), width = 14, 
         height = 22)
  
  saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis,
                            '_lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltypes_fastmnn.rds'))
  
  
  
}

##########################################
# subset the Chan2019 data according to Hannah's manual selection of cell types
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis,
                           '_lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltypes.rds'))

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE) + NoLegend()

cluster_sels = c(1,2,3,4,6,7,8,10,
                 11,16,17,19,
                 23,24,26,27,29,
                 31,33,35,37,38,39,40)

mm = which(!is.na(match(aa$clusterID, cluster_sels)))
aa = subset(aa, cells = colnames(aa)[mm])

aa <- FindVariableFeatures(aa, selection.method = "dispersion", 
                           mean.function = ExpMean,
                           dispersion.function = LogVMR, 
                           mean.cutoff = c(0.0125, 3),
                           dispersion.cutoff = c(0.5, Inf)
)

aa = FindVariableFeatures(aa, selection.method = 'vst', nfeatures = 2000)
aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 10)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.2)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE) + NoLegend()

ggsave(filename = paste0(resDir, '/umap_Chan2019_selectedCelltypes.pdf'), width = 12, height = 8)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis,
                          '_lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltype.subsets.rds'))


##########################################
# annotate the floorplate  
##########################################
aa = readRDS(paste0(RdataDir, 'seuratObject_', species, version.analysis,
                    '_lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltype.subsets.rds'))

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE)

ggsave(filename = paste0(resDir, '/Ref_Chan2019_selectedCelltypes.pdf'), width = 16, height = 8)

aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.6)

DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)

FeaturePlot(aa, features = c('Shh', 'Foxa2', 'Arx', 'Nkx6-1', 'Ntn'))

cells = which(aa$celltype == 'anterior_primitive_streak')
DimPlot(aa, reduction = "umap", group.by = "celltype", label = TRUE,
        repel = TRUE, raster=FALSE, 
        #cols = c(cols_sel, cols_mouse),
        #order = c('day3_RA.rep1'),
        cells.highlight = cells,
        #shuffle = TRUE,
        cols.highlight = 'red'
) + NoLegend()

ggsave(filename = paste0(resDir, '/Ref_Chan2019_selectedCelltypes_highlightedCelltypes_closeToNode.pdf'), 
       width = 16, height = 8)

celltype_sels = c('anterior_primitive_streak', 'node', 'notochord', 'future_spinal_cord', 'fore/midbrain', 
                  'primitive/definitive_endoderm', 'gut_endoderm')

mm = match(aa$celltype, celltype_sels)
aa = subset(aa, cells = colnames(aa)[which(!is.na(mm))])

aa = FindVariableFeatures(aa, selection.method = 'vst', nfeatures = 2000)
aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.2)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p1 / p2
ggsave(filename = paste0(resDir, '/Ref_Chan2019_selectedCelltypes_subsetting_for_FP.pdf'), 
       width = 16, height = 16)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis,
                          '_lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltype.subsets_APS_.rds'))


## subset again the forbrain
aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis,
                           '_lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltype.subsets_APS_.rds'))


aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 10)


ggs = c('Shh', 'Foxa2', 'Arx', 'Nkx6-1', 'Sox17', 'Pax6', 
        'Gsc', 'Lhx1', 'Cer1', 'Eomes', 'Mixl1', 'Tdgf1', 'T', 
        'Afp', 'Acvr2b', 'Acvr1b', 'Aldh1a2', 'Apc', 'Bmp4', 'Smad2',
        'Aifm2', 'Nodal', 'Foxa3',
        "Hhex", 'Cdx1', 'Foxa1', 'Ctnnb1', 'Gata4', 'Gata6', 'Cdx1', 'Ihh', 'Hesx1')

ggs = ggs[!is.na(match(ggs, rownames(aa)))]

pdf(paste0(resDir, '/FeaturePlot_Markers_celltype.subsets_APS_gutEndoderm.pdf'),
    width =10, height = 8, useDingbats = FALSE)
for(n in 1:length(ggs))
{
  cat(n, '--', ggs[n], '\n')
  p1 = FeaturePlot(aa, features = ggs[n], min.cutoff = 'q5')
  #FeaturePlot(ref.combined, features = 'Foxa2', min.cutoff = 'q5')
  #FeaturePlot(ref.combined, features = 'Sox17', min.cutoff = 'q5')
  plot(p1)
  
}

dev.off()


FeaturePlot(aa, features = c('Shh', 'Foxa2', 'Arx', 'Nkx6-1', 'Sox17', 'Pax6'))

celltype_sels = c('future_spinal_cord', 'fore/midbrain')

mm = match(aa$celltype, celltype_sels)
aa = subset(aa, cells = colnames(aa)[which(!is.na(mm))])

aa = FindVariableFeatures(aa, selection.method = 'vst', nfeatures = 2000)
aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.2)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE)

aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 1)

DimPlot(aa, label = TRUE, repel = TRUE,  raster=FALSE)
ggsave(filename = paste0(resDir, '/Ref_Chan2019_selectedCelltypes_subsetting_for_FP_reclustering.pdf'), 
       width = 16, height = 12)

FeaturePlot(aa, features = c('Shh', 'Foxa2', 'Arx', 'Nkx6-1'))
ggsave(filename = paste0(resDir, '/Ref_Chan2019_selectedCelltypes_subsetting_for_FP_markerGenes.pdf'), 
       width = 16, height = 12)

VlnPlot(aa, features = c('Shh', 'Foxa2', "Arx"))

fp_cells = colnames(aa)[which(aa$seurat_clusters == '15')]

saveRDS(fp_cells, file = paste0(RdataDir, 'mouseGastrulation_Chan2019_FPcells.annotated.rds'))


########################################################
########################################################
# Section II : Process mouse gastrulation data Marioni2019
# # the data from the following paper:
# https://www.nature.com/articles/s41586-019-0933-9
# R pacakge of the data were found: 
# https://bioconductor.org/packages/release/data/experiment/html/MouseGastrulationData.html
########################################################
########################################################
Load_process_MouseGastrulation_Marioni2019 = FALSE
if(Load_process_MouseGastrulation_Marioni2019){
  library(MouseGastrulationData)
  head(AtlasSampleMetadata, n = 3)
  
  sce <- EmbryoAtlasData(type = 'processed', samples = NULL)
  sce
  
  saveRDS(sce, file = paste0(RdataDir, 'EmbryoAtlasData_all36sample.rds'))
  
  sce = readRDS(paste0('../results/dataset_scRNAseq_MouseGastrulationData/Rdata/', 
                       'EmbryoAtlasData_all36sample.rds'))
  
  head(rowData(sce))
  rownames(sce) = rowData(sce)$SYMBOL
  
  head(colData(sce))
  
  #exclude technical artefacts
  singlets <- which(!(colData(sce)$doublet | colData(sce)$stripped))
  
  plot(
    x = reducedDim(sce, "umap")[singlets, 1],
    y = reducedDim(sce, "umap")[singlets, 2],
    col = sce$colour[singlets],
    pch = 19,
    xaxt = "n", yaxt = "n",
    xlab = "UMAP1", ylab = "UMAP2"
  )
  
  sce = sce[, singlets]
  EmbryoCelltypeColours = EmbryoCelltypeColours[colData(sce)$celltype[singlets]]
  
  sce = scuttle::logNormCounts(sce)
  rownames(sce) = make.unique(rownames(sce))
  srat = as.Seurat(sce, counts = "counts",  assay = NULL)
  
  rm(sce)
  
  DimPlot(srat, reduction = 'umap', cols = EmbryoCelltypeColours, group.by = 'celltype', 
          label = TRUE, repel = TRUE,
          raster=FALSE)
  
  ggsave(filename = paste0(resDir, '/umap_celltypes.pdf'), width = 20, height = 10)
  
  saveRDS(srat, file = paste0(RdataDir, 'seuratObject_EmbryoAtlasData_all36sample.rds'))
  
}

Modify.Assay.Name = FALSE
if(Modify.Assay.Name){
  srat = readRDS(file = paste0('../results/dataset_scRNAseq_MouseGastrulationData/Rdata/',
                               'seuratObject_EmbryoAtlasData_all36sample.rds'))
  
  ## change assay name
  adt.data <- GetAssayData(object =  srat[['originalexp']], slot = 'counts')
  srat[["RNA"]] <- CreateAssayObject(counts = adt.data )
  DefaultAssay(srat) <- "RNA"
  srat[['originalexp']] = NULL
  
  srat <- NormalizeData(srat, normalization.method = "LogNormalize")
  srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 3000)
  srat <- ScaleData(srat, verbose = FALSE)
  srat <- RunPCA(srat, verbose = FALSE)
  
  saveRDS(srat, file = paste0(RdataDir,  'seuratObject_EmbryoAtlasData_all36sample_RNAassay.rds'))
  
  
}

##########################################
# keep only relevant cell types 
##########################################
Filter_unrelevant_celltype_Marioni2019 = FALSE
if(Filter_unrelevant_celltype_Marioni2019){
  srat = readRDS(file = paste0(RdataDir, 
                               'seuratObject_EmbryoAtlasData_all36sample_RNAassay.rds'))
  xx = readRDS(file = paste0('../results/dataset_scRNAseq_MouseGastrulationData/Rdata/',
                             'seuratObject_EmbryoAtlasData_all36sample.rds'))
  
  umap.embedding = xx@reductions$umap@cell.embeddings
  umap.embedding = umap.embedding[match(colnames(srat), rownames(umap.embedding)), ]
  srat[['umap']] = Seurat::CreateDimReducObject(embeddings=umap.embedding,
                                                key='UMAP_',
                                                assay='RNA')
  
  pca.embedding = xx@reductions$pca.corrected@cell.embeddings
  pca.embedding = pca.embedding[match(colnames(srat), rownames(pca.embedding)), ]
  srat[['pca.corrected']] = Seurat::CreateDimReducObject(embeddings=pca.embedding,
                                                key='PCAcorrected_',
                                                assay='RNA')
  
  rm(xx)
  rm(umap.embedding)
  rm(pca.embedding)
  
  p1 = DimPlot(srat, reduction = 'umap', 
               #cols = EmbryoCelltypeColours, 
               group.by = 'celltype', 
               label = TRUE, repel = TRUE,
               raster=FALSE) 
  
  p2 = DimPlot(srat, reduction = 'umap', 
               #cols = EmbryoCelltypeColours, 
               group.by = 'stage', 
               label = TRUE, repel = TRUE,
               raster=FALSE) 
  
  p1 /p2
  
  ggsave(filename = paste0(resDir, '/MouseGastrulation_Marioni2019_all.celltypes.pdf'), width = 18, height = 20)
  
  saveRDS(srat, file = paste0(RdataDir,
                              'seuratObject_MouseGastrulation_Marioni2019_all36sample_RNAassay_all.celltypes.rds'))
  
  
  ## filter unlikely celltypes in the reference
  sels = grep('Erythroid|Blood|Allantois|mesoderm|Haemato|Cardiomy|Endothelium|Mesenchyme|ExE', srat$celltype, 
              invert = TRUE)
  srat = subset(srat, cells = colnames(srat)[sels])
  
  #saveRDS(srat, file = paste0(RdataDir,  
  #                            'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes.rds'))
  
  #srat = readRDS(file = paste0(RdataDir,  
  #                             'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes.rds'))
  
  sels = grep('Parietal', srat$celltype, 
              invert = TRUE)
  srat = subset(srat, cells = colnames(srat)[sels])
  
  saveRDS(srat, file = paste0(RdataDir,  
                              'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v3.rds'))
  
  p1 = DimPlot(srat, reduction = 'umap', 
               #cols = EmbryoCelltypeColours, 
               group.by = 'celltype', 
               label = TRUE, repel = TRUE,
               raster=FALSE) 
  
  p2 = DimPlot(srat, reduction = 'umap', 
               #cols = EmbryoCelltypeColours, 
               group.by = 'stage', 
               label = TRUE, repel = TRUE,
               raster=FALSE) 
  
  p1 /p2
  
  ggsave(filename = paste0(resDir, '/MouseGastrulation_celltypes_stage.pdf'), width = 18, height = 20)
  
  
}

##########################################
# annot FP in the Marioni2019 dataset 
##########################################
Annot_FP_Marioni219 = FALSE
if(Annot_FP_Marioni219){
  aa = readRDS(file = paste0('../results/scRNAseq_R13547_10x_mNT_20220813/mapping_to_MouseGastrulationData/Rdata/',  
                             'seuratObject_EmbryoAtlasData_all36sample_RNAassay_keep.relevant.celltypes_v2.rds'))
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE)
  
  ggsave(filename = paste0(resDir, '/Ref_Marioni2019_originalUMAP_selectedCelltypes.pdf'), width = 16, height = 8)
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  aa <- ScaleData(aa, verbose = FALSE)
  aa <- RunPCA(aa, verbose = FALSE)
  
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.2)
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE)
  
  aa <- FindNeighbors(aa, dims = 1:20)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 1.0)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)
  
  p1 /p2
  
  ggsave(filename = paste0(resDir, '/Ref_Marioni2019_selectedCelltypes_redoUMAP.clustering.pdf'), 
         width = 16, height = 20)
  
  
  FeaturePlot(aa, features = c('Shh', 'Foxa2', 'Arx', 'Nkx6-1', 'Ntn'))
  
  ggsave(filename = paste0(resDir, '/Ref_Marioni2019_selectedCelltypes_FPmarkers.pdf'), 
         width = 16, height = 8)
  
  celltype_sels = c('Anterior Primitive Streak', 'Def. endoderm', 
                    'Forebrain/Midbrain/Hindbrain', 'Gut', 'Notochord', 
                    'Primitive Streak', 'Spinal cord')
  
  mm = match(aa$celltype, celltype_sels)
  aa = subset(aa, cells = colnames(aa)[which(!is.na(mm))])
  
  aa = FindVariableFeatures(aa, selection.method = 'vst', nfeatures = 2000)
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.2)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'stage', raster=FALSE)
  p1 / p2
  
  ggsave(filename = paste0(resDir, '/Ref_Marioni2019_selectedCelltypes_subsetting_for_FP.pdf'), 
         width = 16, height = 16)
  
  saveRDS(aa, file = paste0(RdataDir, 'seuratObject_mouseGastrulation_Marioni2019',
                            '_celltype.subsets_APS.rds'))
  
  
  ## subset again the forbrain
  aa = readRDS(file = paste0(RdataDir, 'seuratObject_mouseGastrulation_Marioni2019',
                             '_celltype.subsets_APS.rds'))
  
  
  aa <- FindNeighbors(aa, dims = 1:20)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 10)
  
  
  ggs = c('Shh', 'Foxa2', 'Arx', 'Nkx6-1', 'Sox17', 'Pax6', 
          'Gsc', 'Lhx1', 'Cer1', 'Eomes', 'Mixl1', 'Tdgf1', 'T', 
          'Afp', 'Acvr2b', 'Acvr1b', 'Aldh1a2', 'Apc', 'Bmp4', 'Smad2',
          'Aifm2', 'Nodal', 'Foxa3',
          "Hhex", 'Cdx1', 'Foxa1', 'Ctnnb1', 'Gata4', 'Gata6', 'Cdx1', 'Ihh', 'Hesx1')
  
  ggs = ggs[!is.na(match(ggs, rownames(aa)))]
  
  
  pdf(paste0(resDir, '/Ref_Marioni2019_FeaturePlot_Markers_celltype.subsets_APS_gutEndoderm.pdf'),
      width =10, height = 8, useDingbats = FALSE)
  for(n in 1:length(ggs))
  {
    cat(n, '--', ggs[n], '\n')
    p1 = FeaturePlot(aa, features = ggs[n], min.cutoff = 'q5')
    #FeaturePlot(ref.combined, features = 'Foxa2', min.cutoff = 'q5')
    #FeaturePlot(ref.combined, features = 'Sox17', min.cutoff = 'q5')
    plot(p1)
    
  }
  
  dev.off()
  
  
  FeaturePlot(aa, features = c('Shh', 'Foxa2', 'Arx', 'Nkx6-1', 'Sox17', 'Pax6'))
  
  celltype_sels = c('Forebrain/Midbrain/Hindbrain', 'Spinal cord')
  
  mm = match(aa$celltype, celltype_sels)
  aa = subset(aa, cells = colnames(aa)[which(!is.na(mm))])
  
  aa = FindVariableFeatures(aa, selection.method = 'vst', nfeatures = 2000)
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.2)
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE)
  
  
  aa <- FindNeighbors(aa, dims = 1:20)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 1)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE,  raster=FALSE)
  p1 /p2
  
  ggsave(filename = paste0(resDir, '/Ref_Marioni2019_selectedCelltypes_subsetting_for_FP_reclustering.pdf'), 
         width = 16, height = 20)
  
  FeaturePlot(aa, features = c('Shh', 'Foxa2', 'Arx', 'Nkx6-1'))
  ggsave(filename = paste0(resDir, '/Ref_Chan2019_selectedCelltypes_subsetting_for_FP_markerGenes.pdf'), 
         width = 16, height = 12)
  
  VlnPlot(aa, features = c('Shh', 'Foxa2', "Arx"))
  
  fp_cells = colnames(aa)[which(aa$seurat_clusters == '8')]
  
  saveRDS(fp_cells, file = paste0(RdataDir, 'mouseGastrulation_Marioni2019_FPcells.annotated.rds'))
  
}


########################################################
########################################################
# Section III:
# compare those two reference datasets and search for Foxa2-Pax6 double positive clusters
########################################################
########################################################
Compare_Refs_Marioni2019_Chan2019 = FALSE
if(Compare_Refs_Marioni2019_Chan2019){
  
  srat = readRDS(file = paste0(RdataDir,
                              'seuratObject_MouseGastrulation_Marioni2019_all36sample_RNAassay_all.celltypes.rds'))
    
  srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-")
  
  pdfname = paste0(resDir, '/QCs_nCounts_nFeatures_percentMT_Marioni.dataset.pdf')
  pdf(pdfname, width=16, height = 8)
  
  #Idents(srat) = factor(srat$stage)
  
  VlnPlot(srat, features = 'nFeature_RNA', group.by = 'sample', y.max = 6000, raster = FALSE) + NoLegend() +
    geom_hline(yintercept = c(1000, 2000, 3000))
  VlnPlot(srat, features = 'nCount_RNA', group.by = 'sample', y.max = 50000, raster = FALSE) + NoLegend() +
    geom_hline(yintercept = c(10000, 20000, 30000))
  VlnPlot(srat, features = 'percent.mt',  group.by = 'sample', y.max =3, raster = FALSE) + NoLegend()
  
  dev.off()
  
  
  aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis,
                             '_lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltypes.rds'))
  
  pdfname = paste0(resDir, '/QCs_nCounts_nFeatures_percentMT_mouseGastrulation_Chan2019.pdf')
  pdf(pdfname, width=16, height = 8)
  
  VlnPlot(aa, features = 'nFeature_RNA', group.by = 'condition', y.max = 10000, raster = FALSE) + NoLegend() +
    geom_hline(yintercept = c(1000, 2000, 3000, 4000))
  VlnPlot(aa, features = 'nCount_RNA', group.by = 'condition', y.max = 60000, raster = FALSE) + NoLegend() +
    geom_hline(yintercept = c(10000, 20000, 30000))
  VlnPlot(aa, features = 'percent.mt',  group.by = 'condition', y.max =10, raster = FALSE) + NoLegend()
  
  dev.off()
  
}

Search_FoxA2_Pax6_doublePositive = FALSE
if(Search_FoxA2_Pax6_doublePositive){
 
  srat = readRDS(file = paste0(RdataDir,
                               'seuratObject_MouseGastrulation_Marioni2019_all36sample_RNAassay_all.celltypes.rds'))
  
  aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis,
                             '_lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltypes.rds'))
  
  
  ##########################################
  # check the Pax6 and FoxA2 co-expression
  ##########################################
  p1 = VlnPlot(srat, features = c('Pax6'), group.by = 'celltype') + NoLegend()
  p2 = VlnPlot(srat, features = 'Foxa2', group.by = 'celltype') + NoLegend()
  
  p1 / p2
  
  ggsave(filename = paste0(resDir, '/Vlnplot_Pax6_Foxa2_celltypes.pdf'), width = 14, height = 10)
  
  
  FeaturePlot(srat, features = c('Pax6', 'Foxa2'), blend = TRUE, raster = TRUE, order = TRUE)
  
  ggsave(filename = paste0(resDir, '/Featureplots_Pax6_Foxa2_blended_ordered.pdf'), width = 20, height = 8)
  
  
  ## subset potential clusters
  Idents(srat) = as.factor(srat$celltype)
  
  celltypes_sels = c('Surface ectoderm', 'Spinal cord', 'Rostral neurectoderm',
                     'NMP', 'Forebrain/Midbrain/Hindbrain', 
                     'Caudal neurectoderm', 'Caudal epiblast')
  
  sub.obj = subset(srat, idents = celltypes_sels)
  
  Idents(sub.obj) = as.factor(sub.obj$celltype)
  
  DimPlot(sub.obj, reduction = 'umap', group.by = 'celltype', 
          label = TRUE, repel = TRUE,
          raster=FALSE)
  
  
  sub.obj <- RunUMAP(sub.obj, reduction = "pca.corrected", dims = 1:50, 
                     n.neighbors = 30, min.dist = 0.1)
  
  p1 = DimPlot(sub.obj, reduction = 'umap', group.by = 'celltype', 
               label = TRUE, repel = TRUE,
               raster=FALSE)
  p2 = DimPlot(sub.obj, reduction = 'umap', group.by = 'stage', 
               label = TRUE, repel = TRUE,
               raster=FALSE)
  
  p1 + p2
  
  ggsave(filename = paste0(resDir, '/umap_subsetting_celltypes_stages.pdf'), width = 25, height = 8)
  
  FeaturePlot(sub.obj, features = c('Pax6', 'Foxa2'), blend = TRUE, raster = TRUE, order = TRUE)
  ggsave(filename = paste0(resDir, '/Featureplots_Pax6_Foxa2_blended_ordered_subsetting.pdf'), 
         width = 20, height = 8)
  
  
  
   
}

########################################################
########################################################
# Section IV : harmonize the two dataset as unified reference
# 
########################################################
########################################################
##########################################
# merge first two atlas
##########################################
source(paste0(functionDir, 'functions_dataIntegration.R'))
aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis,
              '_lognormamlized_var.to.regress.nCount.RNA_pca_clusterIDs_celltypes_fastmnn.rds'))

aa$celltype[which(is.na(aa$celltype))] = 'unknown'

srat = readRDS(file = paste0(RdataDir, 
                             'seuratObject_EmbryoAtlasData_all36sample_Marioni_pca.corrected.umap.rds'))

srat[['mnt']] = srat[['pca.corrected']]

## merge two references
aa$dataset = 'Chan2019'
srat$dataset = 'Marioni2019'
features.common = intersect(rownames(aa), rownames(srat))

xx = subset(aa, features = features.common)
yy = subset(srat, features = features.common)

yy[['mnn']] = yy[['pca.corrected']]
yy[['umap.mnn']] = yy[['umap.orig']]

rm(list = c('aa', 'srat'))

refs = merge(yy, y = xx, add.cell.ids = c("Marioni2019", "Chan2019"), project = "mouseGastrulation",
             merge.dr = c('pca', 'mnn', 'umap.mnn'))

rm(list = c('xx', 'yy'))
#refs = subset(refs, features = features.common)

saveRDS(refs, file = paste0(RdataDir, 
                            'seuratObject_mouseGastrulationAtlasData_Marioni2019_Chan2019_mnn.rds'))


##########################################
# benchmark different integration methods
##########################################

