##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: analyze the scRNA-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Jul 22 08:12:08 2022
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R13547_10x_mNT_20220813'

resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R13547_10x'
functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

library(pryr) # monitor the memory usage
require(ggplot2)

require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(DropletUtils)
library(edgeR)
library(future)
options(future.globals.maxSize = 80000 * 1024^2)
mem_used()

species = 'mNT_scRNAseq'

### Parameters for this project
levels = c("day0_beforeRA", "day1_beforeRA", 
           "day2_beforeRA",
           "day2.5_RA", "day3_RA.rep1", "day3_RA.rep2", 'day3.5_RA',
           "day4_RA", "day5_RA", "day6_RA",
           "day2.5_noRA", "day3_noRA", 'day3.5_noRA', "day4_noRA", "day5_noRA", "day6_noRA")

# manually set colors by Hannah
library(RColorBrewer)
library("viridis")
cols = rep(NA, length = 16)
names(cols) = levels
cols[grep('_beforeRA', names(cols))] = colorRampPalette((brewer.pal(n = 3, name ="Greys")))(3)
#cols[1:3] = viridis(3)
cols[grep('_noRA', names(cols))] = colorRampPalette((brewer.pal(n = 6, name ="Blues")))(6)
cols[grep('_RA', names(cols))] = colorRampPalette((brewer.pal(n = 7, name ="OrRd")))(7)


########################################################
########################################################
# Section I : import data and quick check QCs
# 
########################################################
########################################################
design = read.csv(file = paste0(dataDir, '/sampleInfos.csv'))
colnames(design) = c('sampleID', 'condition')

design = design[, c(1:2)]

# design = design[which(design$sampleID != '196323'), ]

# import data from cellranger output
for(n in 1:nrow(design))
{
  # n = 1
  cat(n, ' : ', design$condition[n], '\n')
  
  topdir = paste0(dataDir, '/cellRanger_outs/',  design$sampleID[n], '/outs/raw_feature_bc_matrix/')
  exp = Matrix::readMM(paste0(topdir, "matrix.mtx.gz")) #read matrix
  bc = read.csv(paste0(topdir, "/barcodes.tsv.gz"), header = F, stringsAsFactors = F)
  g = read.csv(paste0(topdir, "/features.tsv.gz"), header = F, stringsAsFactors = F, sep = '\t')
  
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
  iscell_dd = defaultDrops(count.data, expected = 8000) # default cell estimate, similar to 10x cellranger
  sum(iscell_dd, na.rm=TRUE)
  
  ## not used the emptyDrops too slow 
  # eout = emptyDrops(count.data, lower = 200)
  # eout$FDR[is.na(eout$FDR)] = 1
  # iscell_ed = eout$FDR<=0.01
  # sum(iscell_ed, na.rm=TRUE)
  
  meta = data.frame(row.names = colnames(count.data), condition = design$condition[n],
                    iscell_dd = iscell_dd)
  
  # plot rankings for number of UMI
  br.out <- barcodeRanks(count.data)
  
  pdf(paste0(resDir, "/UMIrank_emptyDrop_", design$condition[n], "_", design$sampleID[n],  ".pdf"), 
      height = 6, width =10, useDingbats = FALSE)
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")

  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  abline(v = sum(iscell_dd), col = 'darkgreen', lwd = 2.0)
  abline(v = c(3000, 5000, 8000, 10000, 12000), col = 'gray')
  text(x = c(3000, 5000, 8000, 10000, 12000), y =10000, labels = c(3000, 5000, 8000, 10000, 12000), 
       col = 'red')
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
         legend=c("knee", "inflection"))
  
  dev.off()
  
  # use defaultDrop to select cells.
  aa = CreateSeuratObject(counts = count.data[, iscell_dd],
                            meta.data = meta[iscell_dd, ], 
                            min.cells = 20, min.features = 100)
  aa$cell.id = paste0(colnames(aa), '_', design$condition[n], '_', design$sampleID[n])
  
  if(n == 1) {
    mnt = aa
  }else{
    mnt = merge(mnt, aa)
  }
}

mnt[["percent.mt"]] <- PercentageFeatureSet(mnt, pattern = "^mt-")

save(design, mnt, 
     file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, '.Rdata'))

##########################################
# QCs and cell filtering  
##########################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, '.Rdata'))
aa = mnt
rm(mnt)

pdfname = paste0(resDir, '/QCs_nCounts_nFeatures_percentMT.pdf')
pdf(pdfname, width=16, height = 8)

table(aa$condition) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(condition = factor(Var1, levels=levels)) %>%
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
VlnPlot(aa, features = 'percent.mt', y.max =20)


FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "percent.mt")

#ggs = sapply(rownames(aa), function(x) {x = unlist(strsplit(x, '-')); x = x[grep('ENSMUSG', x, invert = TRUE)];
#                                                                               paste0(x, collapse = '-')})

features = c('Pou5f1', 'Sox2', 'Lef1', 'Otx2', 'Zfp703', 'Pax6', 'Foxa2', 'Shh', 'Nkx6-1', 'Nkx2-2', 'Olig2', 
             'Sox1', 'Tubb3')
features = rownames(aa)[!is.na(match(rownames(aa), features))]

Idents(aa) = factor(aa$condition, levels = levels)

for(n in 1:length(features))
{
  p1 = VlnPlot(aa, features = features[n])
  plot(p1)
}

dev.off()

rm(count.data); rm(br.out); rm(bc);rm(meta)

##########################################
## filter cells here
##########################################
VlnPlot(aa, features = 'nFeature_RNA', y.max = 10000) +
  geom_hline(yintercept = c(1000, 8000), col = 'red')

VlnPlot(aa, features = 'percent.mt', y.max =20) + 
  geom_hline(yintercept = c(5), col = 'red')

aa <- subset(aa, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 5)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_', species, version.analysis, '.rds'))

########################################################
########################################################
# Section II : PCA and UMAP exploration for cleaned cells
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_', species, version.analysis, '.rds'))

aa$condition[which(aa$condition == 'day0')] = 'day0_beforeRA'

Idents(aa) = factor(aa$condition, levels = levels)

aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

plan()
# change the current plan to access parallelization
plan("multicore", workers = 16)
plan()

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(aa)

aa <- ScaleData(aa, features = all.genes)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)

ElbowPlot(aa, ndims = 30)

aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

plan("multiprocess", workers = 1)
plan()

Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                          '_lognormamlized_pca_umap_v3_regressed.pct.mt.rds'))


ggsave(filename = paste0(resDir, '/first_test_umap_v3.pdf'), width = 10, height = 8)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_lognormamlized_pca_umap_v3.rds'))

Explore.umap.parameters = FALSE
if(Explore.umap.parameters){
  source(paste0(functionDir, '/functions_scRNAseq.R'))
  explore.umap.params.combination(sub.obj = aa, resDir = resDir, 
                                  pdfname = 'axolotl_spliced_unspliced_umap_test.pdf',
                                  use.parallelization = TRUE,
                                  group.by = 'condition',
                                  cols = cols, 
                                  nfeatures.sampling = c(3000, 5000, 8000),
                                  nb.pcs.sampling = c(30, 50, 100), 
                                  n.neighbors.sampling = c(30, 100, 200),
                                  min.dist.sampling = c(0.1, 0.3)
                                  
                                  )
  
}

########################################################
########################################################
# Section III : Data exploration after cleaning steps
# after cell filtering and doublet detection
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_doubletFinderOut.v2_', 
                           species, version.analysis, '.rds'))

Idents(aa) = factor(aa$condition, levels = levels)

aa[["percent.rb"]] <- PercentageFeatureSet(aa, pattern = "^Rp[sl]")


# discard the doublet from doubletFinder
aa = subset(aa, DF_out == 'Singlet')

aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(aa)

aa <- ScaleData(aa, features = all.genes)

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
ElbowPlot(aa, ndims = 30)

Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)

saveRDS(aa, file = paste0(RdataDir, 
                    'seuratObject_merged_cellFiltered_doublet.rm_', species, version.analysis, '_v3.rds'))


##########################################
# there are some discrepency between those two replicates: day3.5_RA
# and they are not due to normalization after testing multiple normalization
# finally solve the issue by regressing out the nCounts
##########################################
# source('analysis_twoRechnicalReps.R') # test how to remove the discrepency of two technical replicates
# source('script_regressOut.nCount_RNA.R') # rm ribo and mt gene and regress out nCounts

# load the result
aa = readRDS(file = paste0(RdataDir, 
              'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_', 
              species, version.analysis, '_cbe_v4.rds'))

Idents(aa) = factor(aa$condition, levels = levels)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)

head(grep('^Rp[sl]|^mt-', rownames(aa))) # double check if ribo and mito genes 
     
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.2)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)

ggsave(filename = paste0(resDir, 
                         '/umap_rmDoublet_rmRiboMT_regressed.nCounts_5000features_30pcs_100neighbor_dist0.2.pdf'), 
       width = 10, height = 8)


##########################################
# calculate the cell cycle scores for the downstream analysis
# for the moment, don't regress out S.Score or G2M.Score
##########################################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- firstup(cc.genes$s.genes)
s.genes[which(s.genes == 'Mlf1ip')] = 'Cenpu'
g2m.genes <- firstup(cc.genes$g2m.genes)

aa <- CellCycleScoring(aa, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# view cell cycle scores and phase assignments
head(aa[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(aa, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
aa <- RunPCA(aa, features = c(s.genes, g2m.genes))
DimPlot(aa, reduction = 'pca', group.by = 'Phase')


saveRDS(aa, file = paste0(RdataDir, 
                      'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                      'cellCycleScoring_', species, version.analysis, '.rds'))

##########################################
# test to remove the cell cycle effect 
##########################################

##########################################
# explore the parameters for umap of all data 
##########################################
aa =  readRDS(file = paste0(RdataDir, 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', species, version.analysis, '.rds'))

Explore.umap.parameters = FALSE
if(Explore.umap.parameters){
  source(paste0(functionDir, '/functions_scRNAseq.R'))
  explore.umap.params.combination(sub.obj = aa, resDir = resDir, 
                                  pdfname = 'axolotl_umap_test_allSamples_cellCycleScoring_annot.v1.pdf',
                                  use.parallelization = FALSE,
                                  group.by = 'condition',
                                  cols = cols, 
                                  weight.by.var = TRUE,
                                  nfeatures.sampling = c(3000, 5000, 8000),
                                  nb.pcs.sampling = c(30, 50), 
                                  n.neighbors.sampling = c(30, 50, 100),
                                  min.dist.sampling = c(0.1, 0.3)
                                  
  )
  
  
}

########################################################
########################################################
# Section IV : cell annotations
# 
########################################################
########################################################
##########################################
# check the cluster-specific markers to figure out the cell populations
##########################################
### cluster markers
Run_FindAllMarkers = FALSE
if(Run_FindAllMarkers){
    
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
  all.genes <- rownames(aa)
  
  aa <- ScaleData(aa, features = all.genes)
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
  
  ElbowPlot(aa, ndims = 30)
  
  aa <- FindNeighbors(aa, dims = 1:30)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 1.0)
  
  markers = FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  saveRDS(markers, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_markers_v1.rds'))
  
}

markers = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_markers_v2.rds')) 

markers %>%
  filter(!str_detect(gene, '^(AMEX|LOC)')) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) -> top10

#saveRDS(top10, file = paste0(RdataDir, 'top10_markerGenes_coarseCluster.rds'))

xx = subset(aa, downsample = 500)
DoHeatmap(xx, features = top10$gene) + NoLegend()

ggsave(filename = paste0(resDir, '/first_test_clusterMarkers_v2.pdf'), width = 45, height = 40)

##########################################
# manual annotation and systematic annotation 
##########################################
# source('analysis_celltype_manualAnnotation.R')
# source('analysis_cellAnnot_BriscoeModule.R')

########################################################
########################################################
# Section V : what does RA do ? 
# by compare RA vs noRA at day2.before RA, day2.5 with/without RA day3, day3.5 and day4
# 
########################################################
########################################################
source('analysis_1st_bifurcation_RA.vs.noRA.R')


########################################################
########################################################
# Section VI : heterogeneity quantification by RA
# 
########################################################
########################################################



########################################################
########################################################
# Section VII: symmetry breaking by RA 
# after cell annotation, we will focus on the RA condition, in particular d2-d4 
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 
                          'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene',
                          '_regressed.nCounts_clusterd0.7_rmBloodCells', 
                           species, version.analysis, '.rds'))
Idents(aa) = factor(aa$condition, levels = levels)


sels = grep('_RA|day2_beforeRA', aa$condition)

aa = subset(aa, cells = colnames(aa)[sels])
aa$condition = droplevels(aa$condition)
cols_sel = cols[match(levels(aa$condition), names(cols))]

aa = FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE, npcs = 100)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.2)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', pt.size = 2, cols = cols_sel, raster=FALSE)

ggsave(filename = paste0(resDir, '/RA_umap_rmDoublet_5000features_30pcs_100neighbor_dist0.2.pdf'), 
       width = 10, height = 8)

saveRDS(aa, file = paste0(RdataDir, 
                           'savedUMAP_seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene',
                           '_regressed.nCounts_clusterd0.7_rmBloodCells_', 
                           species, version.analysis, '.rds'))

FeaturePlot(aa, features = 'Dhrs3') + scale_y_reverse()

FeaturePlot(aa, features = 'Rarg') + scale_y_reverse()

Explore.umap.parameters = FALSE
if(Explore.umap.parameters){
  source(paste0(functionDir, '/functions_scRNAseq.R'))
  explore.umap.params.combination(sub.obj = aa, resDir = resDir, 
                                  pdfname = 'RAtreated_umap_parameter_test_v2.pdf',
                                  use.parallelization = FALSE,
                                  group.by = 'condition',
                                  cols = cols_sel, 
                                  nfeatures.sampling = c(3000, 5000, 8000),
                                  nb.pcs.sampling = c(30, 50, 100), 
                                  n.neighbors.sampling = c(30, 50, 100),
                                  min.dist.sampling = c(0.2, 0.3, 0.5))
  
}

##########################################
# highlight mark genes in RA treated umap
##########################################
outDir = paste0(resDir, '/geneExamples_in_savedUMAP_RAtreated/')
system(paste0('mkdir -p ', outDir))


features_1 = unique(c('Sox2', 'Sox1', 'Tubb3', 'Elavl3', 
                    'Irx3', 'Irx5', 'Pax3', 'Pax7',
                    'Pax6', 'Olig2', 'Nkx2-9', 'Nkx2-2', 
                    'Nkx6-1', 'Foxa2', 'Arx', 'Shh'
)) # DV overview

features_2 = c('Dhrs3', 'Rarg', 'Cyp26a1',
               'Pou3f1', 'Hoxa1', 'Gas1', 'Spry4', 'Sox11', 
               'Cdh1', 'Cdh2', 'Shh', 'Rfx4', 'Zfp42', 'Tcf15', 'Prrx2', 'Gdf3',
               'Etv5', 'Fgf4', 'Otx2', 'Zscan10', 'Apoe', 'Peg10', 'Klf9', 'Tshz1', 'Skil', 'Zfp703')
features = unique(c(c('Pax6', 'Foxa2', 'Sox1', 'Sox2', 'Tubb3', 'Shh', 'Arx',
                    'Zfp703', 'Lef1', 'Irx5', 'Pou5f1', 'Otx2', 'Adgra2', 'Hoxb4', 
                    'Nkx2-2', 'Nkx2-9', 'Nkx6-1', 'Olig2', 'Pax3', 'Pax7'), 
                    features_1,
                    features_2)) # marker wanted by Hannah


for(n in 1:length(features))
{
  cat(n, '--', features[n], '\n')
  if(length(which(rownames(aa) == features[n])) == 1){
    p1 = FeaturePlot(aa, features = features[n], cols = c('gray', 'red'))  
    plot(p1)
    ggsave(filename = paste0(outDir, '/RAtreated_savedUMAP_', features[n], '.pdf'), 
           width = 10, height = 8)
  }
}
