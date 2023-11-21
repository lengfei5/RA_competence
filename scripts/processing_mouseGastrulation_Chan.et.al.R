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

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_dataIntegration.R')

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
  
  aa$condition = factor(aa$condition)
  source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_dataIntegration.R')
  xx = IntegrateData_runFastMNN(aa, group.by = 'condition', nfeatures = 3000,
                                ndims = c(1:50), 
                                merge.order = c("E8.5_1ab", "E8.0_1ab","E7.5_2", "E7.5_1", "E7.0_1", "E6.5_1"),
                                correct.all = TRUE,
                                reference = NULL)
  
  aa = xx
  rm(xx)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltype', raster=FALSE) + NoLegend()
  
  p1 / p2  
  
  ggsave(filename = paste0(resDir, '/umap_tranferred_celltypes_fastMNN.batchCorrected.pdf'), width = 14, 
         height = 22)
  
}

########################################################
########################################################
# Section II : harmonize the two dataset as unified reference
# 
########################################################
########################################################
srat = readRDS(file = paste0('../results/scRNAseq_R13547_10x_mNT_20220813/mapping_to_MouseGastrulationData/Rdata',  
                             '/seuratObject_EmbryoAtlasData_all36sample_RNAassay.rds'))
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
                                              key='PCACorrected_',
                                              assay='RNA')

rm(list = c("xx", "pca.embedding", "umap.embedding"))

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


srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-")

pdfname = paste0(resDir, '/QCs_nCounts_nFeatures_percentMT_Marioni.dataset.pdf')
pdf(pdfname, width=16, height = 8)

Idents(srat) = factor(srat$stage)

VlnPlot(srat, features = 'nFeature_RNA', y.max = 10000, raster = FALSE)
VlnPlot(srat, features = 'nCount_RNA', y.max = 100000, raster = FALSE)
VlnPlot(srat, features = 'percent.mt', y.max =50, raster = FALSE)

dev.off()



saveRDS(srat, file = paste0(RdataDir, 'seuratObject_EmbryoAtlasData_all36sample_Marioni.rds'))


aa = readRDS(paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                    '_lognormamlized_pca_umap_clustered.rds'))
