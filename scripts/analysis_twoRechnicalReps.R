##########################################
# compare the two replicates of day3_RA
# there are some differences between those two replicates
# the conclusion is that the difference is due to diverse sequence depth and feature number.
# normalization of Seurat or scran can not resolve this issue.
# regress out nCount_RNA works better for now
##########################################
# doulet removed already here
aa = readRDS(file = paste0(RdataDir, 
                          'seuratObject_merged_cellFiltered_doublet.rm_', 
                          species, version.analysis, '_v3.rds'))

Idents(aa) = factor(aa$condition, levels = levels)

outDir = paste0(resDir, '/twoReps_issue_v2')
system(paste0('mkdir -p ', outDir))

FeatureScatter(aa, feature1="percent.rb", feature2="nFeature_RNA")
FeatureScatter(aa, feature1="percent.rb", feature2="percent.mt")
FeatureScatter(aa, feature1 = "percent.rb", feature2="nCount_RNA")
FeatureScatter(aa, feature1="nCount_RNA", feature2="nFeature_RNA")
VlnPlot(aa, features = 'percent.rb', y.max =50)

## subset the two replicates and the time point before and after
subs = subset(aa, condition == 'day3_RA.rep1'|condition == 'day3_RA.rep2'| 
                condition == 'day3.5_RA'| condition == 'day2.5_RA')

DimPlot(subs, label = TRUE, repel = TRUE, raster=FALSE,
        reduction = 'umap')
ggsave(filename = paste0(outDir, '/UMAP_RA_day2.5_day3.twoReps_day3.5_beforeRegression.pdf'), 
       width = 10, height = 8)

p1 = VlnPlot(subs, features = 'nFeature_RNA', y.max =10000)
p2 = VlnPlot(subs, features = 'nCount_RNA', y.max = 50000)
p3 = VlnPlot(subs, features = 'percent.mt', y.max =5)
p4 = VlnPlot(subs, features = 'percent.rb', y.max =50)

(p1 + p2)|(p3 + p4) 

ggsave(filename = paste0(outDir, '/QC_features_beforeRegression.pdf'), 
       width = 16, height = 10)

subs_orgi = subset(aa, condition == 'day3_RA.rep1'|condition == 'day3_RA.rep2'| 
                     condition == 'day3.5_RA'| condition == 'day2.5_RA')
Idents(subs_orgi) = subs_orgi$condition

subs_orgi = subset(subs_orgi, downsample = 2000)

##########################################
# Test Seurat regression nCount_RNA  without gene filtering
##########################################
subs = subs_orgi

subs <- NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000)
#subs <- ScaleData(subs)
#subs <- ScaleData(subs, vars.to.regress = 'nFeature_RNA')
subs <- ScaleData(subs, vars.to.regress = 'nCount_RNA')

subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE)
ElbowPlot(subs, ndims = 30)

Idents(subs) = subs$condition
subs <- RunUMAP(subs, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(subs, label = TRUE, repel = TRUE, raster=FALSE,
        reduction = 'umap')

ggsave(filename = paste0(outDir, '/Day3_RA_twoReps.subsample_regressed.nCount.RNA_noGeneFiltering_v3.pdf'), 
       width = 10, height = 8)


##########################################
# Test Seurat regress nFeature_RNA  without gene filtering
##########################################
subs = subs_orgi

subs <- NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000)
#subs <- ScaleData(subs)
subs <- ScaleData(subs, vars.to.regress = 'nFeature_RNA')
#subs <- ScaleData(subs, vars.to.regress = 'nCount_RNA')

subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE)
ElbowPlot(subs, ndims = 30)

Idents(subs) = subs$condition
subs <- RunUMAP(subs, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(subs, label = TRUE, repel = TRUE, raster=FALSE,
        reduction = 'umap')

ggsave(filename = paste0(outDir, '/Day3_RA_twoReps.subsample_regressed.nFeature.RNA_noGeneFiltering_v3.pdf'), 
       width = 10, height = 8)

##########################################
# test mt rb gene filtering + NO regress out nCounts_RNA
##########################################
source(paste0(functionDir, '/functions_scRNAseq.R'))
subs = subs_orgi

subs = subset(subs, features = rownames(subs)[grep('^Rp[sl]|^mt-', rownames(subs), invert = TRUE)])

# gg.rb = rownames(aa)[grep('^Rp[sl]', rownames(aa))]
# gg.mt = rownames(aa)[grep('^mt-', rownames(aa))]
# aa = geneFiltering.scran(aa)

subs <- NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000)
subs <- ScaleData(subs)
#subs <- ScaleData(subs, vars.to.regress = 'nFeature_RNA')
#subs <- ScaleData(subs, vars.to.regress = 'nCount_RNA')

subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE)
ElbowPlot(subs, ndims = 30)

Idents(subs) = subs$condition
subs <- RunUMAP(subs, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(subs, label = TRUE, repel = TRUE, raster=FALSE,
        reduction = 'umap')

ggsave(filename = paste0(outDir, '/Day3_RA_twoReps.subsample_noRegressed.nCount_mt.ribo.GeneFiltering_v3.pdf'), 
       width = 10, height = 8)


##########################################
# filtering mt rb gene and highly expressed genes + NO regress out nCounts_RNA
##########################################
require(SingleCellExperiment)
library(scran)
library(scater)
library(scuttle)
library(Seurat)
library(SeuratObject)

source(paste0(functionDir, '/functions_scRNAseq.R'))
subs = subs_orgi

# subs = subset(subs, features = rownames(subs)[grep('^Rp[sl]|^mt-', rownames(subs), invert = TRUE)])
gg.rb = rownames(subs)[grep('^Rp[sl]', rownames(subs))]
gg.mt = rownames(subs)[grep('^mt-', rownames(subs))]


sce <- as.SingleCellExperiment(subs)

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotHighestExprs(sce, n=50, exprs_values = "counts")
ggsave(filename = paste0(outDir, 
                         '/Day3_RA_twoReps.subsample_mt.ribo.Gene_countPercentage.pdf'), 
       width = 14, height = 6)

#fontsize <- theme(axis.text=element_text(size=16), axis.title=element_text(size=16))
#plotHighestExprs(sce, n=30) + fontsize

ave.counts <- calculateAverage(sce, assay.type = "counts")

hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))

num.cells <- nexprs(sce, byrow=TRUE)

smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
              xlab=expression(Log[10]~"average count"))

# detected in >= 5 cells, ave.counts >=5 but not too high
genes.to.keep <- num.cells > 5 & ave.counts >= 10^-4  & ave.counts <10^2  
summary(genes.to.keep)

# remove mt and ribo genes
genes.to.keep = genes.to.keep & ! rownames(sce) %in% gg.mt & ! rownames(sce) %in% gg.rb
summary(genes.to.keep)

sce <- sce[genes.to.keep, ]

subs = Seurat::as.Seurat(sce)


subs <- NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000)
subs <- ScaleData(subs)
#subs <- ScaleData(subs, vars.to.regress = 'nFeature_RNA')
#subs <- ScaleData(subs, vars.to.regress = 'nCount_RNA')

subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE)
ElbowPlot(subs, ndims = 30)

Idents(subs) = subs$condition
subs <- RunUMAP(subs, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(subs, label = TRUE, repel = TRUE, raster=FALSE,
        reduction = 'umap')

ggsave(filename = 
      paste0(outDir, 
      '/Day3_RA_twoReps.subsample_noRegressed.nCount_mt.ribo.lowly.hihglyExpr.GeneFiltering_v3.pdf'), 
      width = 10, height = 8)

##########################################
# test mt rb gene filtering + regress out nCounts_RNA
##########################################
source(paste0(functionDir, '/functions_scRNAseq.R'))
subs = subs_orgi

FeatureScatter(subs, feature2="percent.rb", feature1="nFeature_RNA")
FeatureScatter(subs, feature1="nFeature_RNA", feature2="percent.mt")
FeatureScatter(subs, feature1 = "percent.rb", feature2="nCount_RNA")

subs = subset(subs, features = rownames(subs)[grep('^Rp[sl]|^mt-', rownames(subs), invert = TRUE)])

# gg.rb = rownames(aa)[grep('^Rp[sl]', rownames(aa))]
# gg.mt = rownames(aa)[grep('^mt-', rownames(aa))]
# aa = geneFiltering.scran(aa)

subs <- NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000)
#subs <- ScaleData(subs)
#subs <- ScaleData(subs, vars.to.regress = 'nFeature_RNA')
subs <- ScaleData(subs, vars.to.regress = 'nCount_RNA')

subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE)
ElbowPlot(subs, ndims = 30)

Idents(subs) = subs$condition
subs <- RunUMAP(subs, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(subs, label = TRUE, repel = TRUE, raster=FALSE,
        reduction = 'umap')

ggsave(filename = paste0(outDir, '/Day3_RA_twoReps.subsample_regressed.nCount_mt.ribo.GeneFiltering_v3.pdf'), 
       width = 10, height = 8)


##########################################
# test mt rb gene filtering + regress out nCounts_RNA for the full sample of day3.5_RA
##########################################
subs = subset(aa, condition == 'day3_RA.rep1'|condition == 'day3_RA.rep2'| 
                condition == 'day3.5_RA'| condition == 'day2.5_RA')

subs = subset(subs, features = rownames(subs)[grep('^Rp[sl]|^mt-', rownames(subs), invert = TRUE)])

# gg.rb = rownames(aa)[grep('^Rp[sl]', rownames(aa))]
# gg.mt = rownames(aa)[grep('^mt-', rownames(aa))]
# aa = geneFiltering.scran(aa)

subs <- NormalizeData(subs, normalization.method = "LogNormalize", scale.factor = 10000)
subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000)
subs <- ScaleData(subs)
#subs <- ScaleData(subs, vars.to.regress = 'nFeature_RNA')
#subs <- ScaleData(subs, vars.to.regress = 'nCount_RNA')

subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE)
ElbowPlot(subs, ndims = 30)

Idents(subs) = subs$condition
subs <- RunUMAP(subs, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(subs, label = TRUE, repel = TRUE, raster=FALSE,
        reduction = 'umap')

ggsave(filename = paste0(outDir, '/Day3_RA_twoReps.fullSample_Regressed.nCount_Nomt.ribo.GeneFiltering_v4.pdf'), 
       width = 10, height = 8)

##########################################
# Test different normalization to solve the replicate issue
##########################################
Test.differentNormalization.for.twoReps = FALSE
if(Test.differentNormalization.for.twoReps){
  subs = subset(aa, condition == 'day3_RA.rep1'|condition == 'day3_RA.rep2'| 
                  condition == 'day3.5_RA'| condition == 'day2.5_RA')
  
  subs = subset(subs, features = rownames(subs)[grep('^Rp[sl]|^mt-', rownames(subs), invert = TRUE)])
  
  ## test scran normalization
  Test.scran.normalization = FALSE
  if(Test.scran.normalization){
    DefaultAssay(subs) = 'RNA'
    tic()
    subs = Normalize_with_scran(subs)
    toc()
    
    subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 3000)
    subs <- ScaleData(subs)
    subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = FALSE, reduction.key = "PC_")
    
    subs <- RunUMAP(subs, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
    
    DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'condition', 
            reduction = 'umap')
    ggsave(filename = paste0(outDir, '/Day3_RA_twoReps.fullSample_mt.ribo.GeneFiltering_scranNorm_v4.pdf'), 
           width = 10, height = 8)
    
  }
  
  Test.GLMPCA = FALSE
  if(Test.GLMPCA){
    require(glmpca)
    
    set.seed(202)
    Y = GetAssayData(object = subs, slot = "counts")
    Y <- Y[rowSums(Y) > 0, ]
    system.time(res1<-glmpca(Y,50,fam="nb", minibatch = 'stochastic')) #about 9 seconds
    
    
  }
  
  ## test downsample method for normalization
  Test.downsample = FALSE
  if(Test.downsample){
    xx1 = subset(aa, condition == 'day3_RA.rep2')
    xx2 = subset(aa, condition != 'day3_RA.rep2')
    
    counts = as.matrix(x = GetAssayData(object = xx1, assay = "RNA", slot = "counts"))
    downsampled = SampleUMI(data = counts, max.umi = median(aa$nCount_RNA[which(aa$condition == 'day3_RA.rep1')]))
    #head(x = downsampled)
    
    xx0 = CreateSeuratObject(counts = downsampled,
                             meta.data = xx1@meta.data)
    
    VlnPlot(xx0, features = 'nCount_RNA')
    
    rm(xx1);
    
    subs = merge(xx0, xx2)
    rm(xx0);  rm(xx2)
    metadata = subs@meta.data
    subs = CreateSeuratObject(counts = GetAssayData(object = subs, assay = "RNA", slot = "counts"),
                              meta.data = metadata[, c(6:7)])
    
    
    Idents(subs) = factor(subs$condition, levels = levels)
    VlnPlot(subs, features = 'nFeature_RNA', y.max = 10000)
    VlnPlot(subs, features = 'nCount_RNA', y.max = 100000)
    
    FeaturePlot(subs, features = 'nCount_RNA')
    
    FeaturePlot(subs, features = 'Pax6')
    VlnPlot(subs, features = 'Pax6')
    VlnPlot(subs, features = 'Foxa2')
    VlnPlot(subs, features = 'Sox2')
    VlnPlot(subs, features = 'Gm28438')
    
    VlnPlot(aa, features = 'Rpl31')
    VlnPlot(aa, features = 'Sox2')
    
    markers = FindMarkers(subs, ident.1 = 'day3_RA.rep2', ident.2 = 'day3_RA.rep1', 
                          min.pct = 0.25, only.pos = FALSE)
    head(markers, n = 5)
    
  }
 
  ## test SCTransform normalization
  Test.SCT.normalization = FALSE
  if(Test.SCT.normalization){
    subs = SCTransform(subs, variable.features.n = 3000, vst.flavor = "v2")
    
    subs = SCTransform(subs, variable.features.n = 3000, vst.flavor = "v2", vars.to.regress = 'nCount_RNA')
    
    subs <- RunPCA(subs, verbose = FALSE, weight.by.var = TRUE)
    ElbowPlot(subs, ndims = 30)
    
    subs <- RunUMAP(subs, dims = 1:30, n.neighbors = 30, min.dist = 0.3)
    DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'condition')
    
    ggsave(filename = paste0(resDir, '/RA_day3.5_twoReps_geneFiltering_SCT.regressed.nCount_RNA.pdf'), 
           width = 10, height = 8)
    
    
    DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'condition', reduction = 'pca')
    DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'condition', dims = c(3, 4), reduction = 'pca')
    
    DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'condition', dims = c(5, 6), reduction = 'pca')
    DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'condition', dims = c(7, 8), reduction = 'pca')
    
    
    
    xx = subs@reductions$pca@feature.loadings
    write.csv(xx, file = paste0(resDir, '/PCA_featureLoading_matrix_day3_RA.rep1_2_day3.5_RA.csv'), 
              row.names = TRUE)
  }
  
  
}
