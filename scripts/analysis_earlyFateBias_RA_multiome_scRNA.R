##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: search for the early bias/heterogeneity before RA, RA and/or RA withdraw
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jan 27 11:53:12 2026
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R16597_mNT_10xmultiome_reseq_20240517'

resDir = paste0("../results/scRNAseq", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/'

source(paste0(functionDir,  'functions_scATAC.R'))
source(paste0(functionDir, 'functions_scRNAseq.R'))
source(paste0(functionDir, 'functions_Visium.R'))

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
require(SeuratObject)
library(data.table)
library(DropletUtils)
library(edgeR)
library(future)
library(tictoc)

library(RColorBrewer)
require(Matrix)
library(RaceID)
require(tictoc)
require(gridExtra)
require(cowplot)

options(future.globals.maxSize = 200 * 1024^3)
set.seed(1234)
mem_used()

species = 'mNT_multiome'

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


## subset the mutliome conditions
levels_sels = c("day2_beforeRA",
                "day2.5_RA", "day3_RA", "day3.5_RA",  "day4_RA", "day5_RA", 
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA")

cc = names(cols)
cc[which(cc == 'day3_RA.rep1')] = 'day3_RA'
names(cols) = cc

cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/mNTs_multiome_RA_',
                'searching_earlyBias')
if(!dir.exists(outDir)) dir.create(outDir)


##########################################
# features  
##########################################
sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

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
                      'Nkx2-2', 'Nkx2-9', 'Nkx6-1', 'Olig2', 'Pax3', 'Pax7', 'Cyp26a1', 'Dhrs3'), 
                    features_1,
                    features_2)) # marker wanted by Hannah

gene_examples = unique(c('Foxa2', 'Pax6', c('Zfp42', 'Tcf15', 'Skil', 'Lef1',
                                            'Sox2', 'Pou5f1', 'Sall4', 'Tdgf1', # pluripotency markers
                                            'Nanog', 'Nr5a2', #'Prdm14', 
                                            'Klf4', 'Fgf4', 'Esrrb', 'Tcf3', 'Tbx3'), # naive pluripotency
                         c('Zfp42', 'Tcf15', 'Skil',
                           'Fgf5', 'Otx2', 'Pou3f1', 'Lef1', 'Dnmt3b', 'Dnmt3a',	
                           'Foxd3', 'Utf1', 'Tcf15', 'Zic3', 'Rhox5', 'Etv5', 'Etv4',	
                           'Lin28b', 'Sox4', 'Sox3', 'Sox11'
                         ),
                         c('Lhx1','Eomes', 'Sox2', 'Hoxb4', 'Hoxb5', 'Hoxb6','Zfp703'),
                         c('Zfp42', 'Tcf15', 'Skil', 'Lef1', 'Dhrs3', 'Rarg', 'Cyp26a1'),
                         features,
                         features_1,
                         features_2
                         
))


########################################################
########################################################
# Section I : check gene examples of Pluripotent, RA response, AP identity genes
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir,
                           'seuratObj_RAmultiome_snRNA.normalized.umap_scATAC.normalized',
                           'noMatureNeurons.smallClusters_wnn.umap.clusters.rds'))

sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
# aa <- RunUMAP(aa, nn.name = "weighted.nn", 
#               reduction.name = "wnn.umap", 
#               reduction.key = "wnnUMAP_", min.dist = 0.1)

FeaturePlot(aa, features = c('Pax6', 'Foxa2', 'Cdh1', 'Cdh2'), reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_Pax6_Foxa2_Cdh1.2.pdf'), 
       height = 12, width = 17)

## # naive pluripotency
genes = c('Nanog', 'Nr5a2', 'Prdm14', 'Klf4', 'Fgf4', 'Esrrb', 'Tcf3', 'Tbx3')

FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_naivPluripotent.pdf'), 
       height = 12, width = 18)


## pluripotent markers
genes = c('Zfp42', 'Tcf15', 'Skil', 'Lef1','Sox2', 'Pou5f1', 'Sall4', 'Tdgf1')

FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_pluripotentMarkers.pdf'), 
       height = 12, width = 18)

## RA response
genes = c('Dhrs3', 'Rarg', 'Cyp26a1', 'Cyp26b1', 'Zfp703', 'Otx2', 'Stra8', 'Tead1', 'Tead2', 'Tead3', 'Tead4',
          'Yap1', 'Tbx1', 'Tbx2', 'T', 'Runx1', 'Rbp4', 'Rxra', 'Rxrb', 'Myb', 'Lefty1', 'Krt13', 
          'Igf2r', 'Gsk3b', 'Gata4', 'Gata6', 'Creb1', 'Asxl1', 'Ascl1', 'Bmp6')

FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_RAresponse.pdf'), 
       height = 36, width = 27)

## Hox genes
genes = rownames(aa)[grep('Hox', rownames(aa))]
genes = c('Hoxa1', 'Hoxa3', 'Hoxb4')

FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_HoxGenes.pdf'), 
       height = 12, width = 18)

## WNT genes
genes = c('Sfrp2', 'Zeb2', 'Zfp703', 'Sox17', 'Snai2', 'Sall1', 'Rspo2', 'Rspo4', 'Mdk', 'Lrrk1', 
          'Lmx1a', 'Lgr4', 'Lgr5', 'Dkk1', 'Dkk4', 'Tcf15', 'Wnt9a', 'Wnt4', 'Lef1')

FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_Wnt.pdf'), 
       height = 21, width = 27)

FeaturePlot(aa, features = c('Lef1', 'Tcf15', 'Zeb2', 'Zfp703', 'Axin2', 'Tcf7', 'Tcf7l1', 'Tcf7l2', 
                             'Dvl1', 'Dvl2'), 
            reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_Wnt_2.pdf'), 
       height = 16, width = 27)

genes = c('Gdf3',
          'Etv5', 'Fgf4', 'Otx2', 'Zscan10', 'Apoe', 'Peg10', 'Klf9', 'Tshz1')
FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_HannahsExamples.pdf'), 
       height = 16, width = 27)

## FGF genes
genes = c('Fgf4', 'Frs2', 'Akt1', 'Spry1', 'Spry2',
          'Etv5', 'Etv4', 'Dusp6', 'Egr1')
FeaturePlot(aa, features = genes, reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_FGFgenes.pdf'), 
       height = 16, width = 27)

## BMP and TGFb genes
genes = c('Bmp2', 'Bmp4', 'Bmp7', 'Gdf5', 'Gdf6', 'Id1', 'Id2', 'Id3', 'Msx1', 'Msx2', 'Fst', 'Smad1', 
          'Smad5', 'Smad4', 'Smad2', 'Smad3', 'Lefty1', 'Lefty2', 'Nog', 'Grem1')

p1 = FeaturePlot(aa, features = genes, reduction = 'wnn.umap')
plot(p1)
ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_BMP.TGFbeta.pdf'), 
       height = 21, width = 27)



########################################################
########################################################
# Section II: VarID method to identify the high expression noise
# 
########################################################
########################################################
outDir = paste0(resDir, '/mNTs_multiome_RA_searching_earlyBias/',
                'gene_Noise_day2_2.5_3_scRNAseq_rmCellCycleGenes')

system(paste0('mkdir -p ', outDir))


cc = c('day2_beforeRA', 'day2.5_RA', 'day3_RA.rep1')
aa = readRDS(file = paste0('../results/Rdata/', 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', '_timePoint_',
                           'mNT_scRNAseq', 
                           '_R13547_10x_mNT_20220813', '.rds'))

p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'clusters', label = TRUE, repel = TRUE)

p1 + p2

DimPlot(aa, group.by = 'Phase', label = TRUE, repel = TRUE)

#ggsave(filename = paste0(outDir, '/UMAP_conditions_clusters.pdf'), 
#       width = 16, height = 6)

Idents(aa) = aa$condition

cat(which(rownames(aa) == 'Dhrs3'), '\n')

set.seed(2023)

# remove the cluster 9 and 10
aa = subset(aa, cells = colnames(aa)[which(!is.na(match(aa$condition, cc))
                                           & as.character(aa$clusters) != '9' 
                                           & as.character(aa$clusters) != '10')])
#aa = subset(aa, cells = colnames(aa)[which(!is.na(match(aa$condition, cc)))])

aa$condition = droplevels(aa$condition)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)
aa <- RunPCA(aa, verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

p2 = DimPlot(aa, group.by = 'Phase', label = TRUE, repel = TRUE) 

p1 + p2

saveRDS(aa, file = paste0(outDir, 'scRNAseq_d2_d2.5_d3_cleaned.rds'))


##########################################
# test clustering and find markers 
##########################################
subclustering_eachTimepoint = FALSE
if(subclustering_eachTimepoint){
  
  aa = readRDS(file = paste0(outDir, 'scRNAseq_d2_d2.5_d3_cleaned.rds'))
  #aa = readRDS(file = paste0(outDir, 'scRNAseq_d2_d2.5_d3_cleaned_dropcellCycleGenes.rds'))
  
  c = "day3_RA.rep1"
  aa = subset(aa, cells = colnames(aa)[which(aa$condition == c)])
  aa$condition = droplevels(aa$condition)
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 1000)
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 30, min.dist = 0.1)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, group.by = 'Phase', label = TRUE, repel = TRUE) 
  
  p1 + p2
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.4)
  
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p1 + p2
  
  #p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'RNA_snn_res.0.5', raster=FALSE)
  #p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'RNA_snn_res.0.7', raster=FALSE)
  #p1 + p2
  
  aa$clusters = aa$seurat_clusters
  #aa$clusters[which(aa$clusters == '7')] = '2'
  #aa$clusters = droplevels(aa$clusters)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  
  p1 + p2
  
  ggsave(filename = paste0(outDir, '/scRNAseq_dropcellCycleGenes_clusters_condition_', c, '.pdf'), 
         height = 12, width = 18)
  
  
  FeaturePlot(aa, features = c('Pax6', 'Foxa2'))
  
  #markers = FindMarkers(aa, ident.1 = c('2'), ident.2 = c('5'))
  
  all.markers <- FindAllMarkers(aa, only.pos = TRUE)
  all.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top10
  
  DoHeatmap(aa, features = top10$gene) + NoLegend()
  
  ggsave(filename = paste0(outDir, '/mNTs_scRNA_clusters_heatmap_markerGenes_condition_', c, '.pdf'), 
         width = 14, height = 24)
  
  FeaturePlot(aa, features = c('Cdh1', 'Sox17'))
  
  ggsave(filename = paste0(outDir, '/scRNAseq_dropcellCycleGenes_featureExamples_', c, '.pdf'), 
         height = 12, width = 18)
  
  
  #FeaturePlot(aa, features = c('Atoh1'))
  
}


##########################################
# gene noise with VarID2
##########################################
aa = readRDS(file = paste0(outDir, 'scRNAseq_d2_d2.5_d3_cleaned.rds'))

#outDir = paste0(resDir, '/mNTs_multiome_RA_searching_earlyBias/gene_Noise_day2_2.5_3_scRNAseq')
#system(paste0('mkdir -p ', outDir))

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, group.by = 'Phase', label = TRUE, repel = TRUE) 

p1 + p2

ggsave(filename = paste0(outDir, '/UMAP_conditions_cellcyclePhase.pdf'), 
       width = 16, height = 6)

FeaturePlot(aa, features = c('Nr5a2', 'Klf4', 'Fgf4', 'Esrrb', 'Tbx3'))

ggsave(filename = paste0(outDir, '/scRNAseq_featurePlot_naivPluripotency.pdf'), 
       height = 12, width = 18)


FeaturePlot(aa, features = c('Zfp42', 'Tcf15', 'Skil', 'Lef1', 'Sox2', 'Pou5f1', 'Sall4', 'Tdgf1'))

ggsave(filename = paste0(outDir, '/scRNAseq_featurePlot_PluripotencyMarkers.pdf'), 
       height = 12, width = 18)


### drop the cell cycle related genes or not
Drop_cellCycleGenes = FALSE
if(Drop_cellCycleGenes){
  
  ## remove cell cycle related genes
  scaledMatrix = GetAssayData(aa, slot = c("scale.data"))
  
  diff <- scater::getVarianceExplained(scaledMatrix, data.frame(phase = aa$Phase))
  diff = data.frame(diff, gene = rownames(diff))
  diff = diff[order(-diff$phase), ]
  
  hist(diff$phase, breaks = 100); abline(v = c(1:5), col = 'red')
  
  genes_discard = diff$gene[which(diff$phase > 5)]
  cat(length(genes_discard), 'genes to discard \n')
  
  aa = subset(aa, features = setdiff(rownames(aa), genes_discard))
  
}

Idents(aa) = aa$condition
aa = subset(aa, downsample = 500)

cat(which(rownames(aa) == 'Dhrs3'), '\n')

x <- aa@assays$RNA@counts
rownames(x) <- rownames(aa)
colnames(x) <- colnames(aa)
sc <- SCseq(x)

## filter the genes 
sc <- filterdata(sc, 
                 mintotal = 1000,
                 minexpr = 1, # 
                 minnumber = 5,
                 FGenes=NULL,
                 CGenes=NULL,
                 verbose = FALSE
)

expData <- getExpData(sc)
which(rownames(expData) == 'Dhrs3')

parallel::detectCores()

tic()
res   <- pruneKnn(expData, no_cores = 32)
toc()

saveRDS(res, file = paste0(outDir, '/RA_d2.beforeRA_no.day6_noNeurons_varID2_pruneKnn_all.rds'))


res = readRDS(file = paste0(outDir, '/RA_d2.beforeRA_no.day6_noNeurons_varID2_pruneKnn_all.rds'))

plotRegNB(expData, res, "(Intercept)")

plotRegNB(expData,res, "beta")

plotRegNB(expData,res,"theta")

plotPearsonRes(res,log=TRUE,xlim=c(-.1,.2))


tic()
cl    <- graphCluster(res, pvalue= 0.01, leiden.resolution = 0.3)
toc()

test = cl$partition
mm = match(names(test), colnames(aa))
newtest = test
newtest[which(aa$condition[mm] == 'day2_beforeRA')] = 1
newtest[which(aa$condition[mm] == 'day2.5_RA')] = 2
newtest[which(aa$condition[mm] == 'day3_RA.rep1')] = 3
cl$partition = newtest

probs <- transitionProbs(res,cl)

plotPC(res)

tic()
sc <- updateSC(sc,res=res,cl=cl)
toc()

#plotmap(sc,fr=TRUE)

#plotmap(sc,fr=TRUE)
# Alternatively, a umap representation can be computed for visualization:
sc <- compumap(sc, min_dist=0.2, n_neighbors =20)
plotmap(sc,um=TRUE)

pdfname = paste0(outDir, '/plot_transition_probability.pdf')
pdf(pdfname, width=10, height = 8)

plotTrProbs(sc,probs,um=TRUE)

dev.off()

save(sc, cl, res, 
     file = paste0(outDir, '/out_varID2_scSeq_pruneKnn_all.Rdata'))

##########################################
# ## compute noise from corrected variance
# https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html#varid2
# The most important argument of the compTBNoise is the gamma parameter 
# which determines the strength of the prior for the suppression of spurious inflation of noise levels
# at low expression
# If a positive correlation is observed, gamma should be increased in order to weaken the prior. 
# If the correlation is negative, gamma should be decreased in order to increase the strength of the prior.
##########################################
for(gamma in c(2, 3, 4, 5, 6)) # search for the optimal value of gamma
{
  # gamma = 5
  cat('gamma -- ', gamma, '\n')
  
  load(file = paste0(outDir, '/out_varID2_scSeq_pruneKnn_all.Rdata'))
  
  #nn <- inspectKNN(20,expData,res,cl,object=sc,pvalue=0.01,plotSymbol=TRUE,um=TRUE,cex=1)
  #head(nn$pv.neighbours)
  #head(nn$expr.neighbours)
  x <- getFilteredCounts(sc, minexpr=1, minnumber=10)
  
  tic()
  noise <- compTBNoise(res, x, gamma = gamma, 
                       pvalue=0.01, no_cores=64) # take 40 minutes, required ~80G memory  
  toc()
  
  
  pdfname = paste0(outDir, '/plot_varID2_variance_mean_gamma_', gamma, '.pdf')
  pdf(pdfname, width=8, height = 6)
  
  plotUMINoise(sc, noise,log.scale=TRUE)
  
  dev.off()
  
  #noise <- compTBNoise(res,expData,pvalue=0.01,no_cores=5) 
  sc <- updateSC(sc,res=res,cl=cl, noise=noise)
  
  #saveRDS(sc, file = paste0(outDir, '/RA_d2.beforeRA_no.day6_noNeurons_varID2_noise_gamma.', gamma,'.rds'))
  save(sc, noise, cl, res, 
       file = paste0(outDir, '/out_varID2_noise_gamma_', gamma, '_all.Rdata'))
  
}


########################################################
# ## reload the results and make analysis
########################################################
# The prior parameter gamma should be chosen such that 
# the correlation between mean noise across cells and total UMI count per cell is minimal. 
# This dependence can be analysed with the function plotUMINoise
gamma = 3
load(file = paste0(outDir, '/out_varID2_noise_gamma_', gamma, '_all.Rdata'))
plotUMINoise(sc,noise,log.scale=TRUE)

plotExpNoise("Foxa2", sc, noise,norm=TRUE, log="xy")

pdfname = paste0(outDir, '/plot_umap_varID_conditions.pdf')
pdf(pdfname, width=8, height = 6)

#p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
#plot(p1)
plotmap(sc, um=TRUE, tp = 1, cex = 0.3)

dev.off()

#tic()
#cl    <- graphCluster(res, pvalue= 0.01, leiden.resolution = 1.5)
#toc()

#sc <- updateSC(sc,res=res,cl=cl)
#plotmap(sc, um=TRUE, tp = 1, cex = 0.3)
#sc <- updateSC(sc,res=res,cl=cl,noise=noise)

plotExpNoise("Foxa2", sc, noise,norm=TRUE,log="xy")


plotexpmap(sc, 'Hdac1', logsc=TRUE, um=TRUE, noise=FALSE, cex=1)


# Gene expression (on logarithmic scale):
plotexpmap(sc, 'Pax6', logsc=TRUE, um=TRUE, noise=FALSE, cex=1)



plotexpmap(sc, 'Pax6', logsc=TRUE, um=TRUE, cex=1, noise = TRUE)

plotexpmap(sc, 'Foxa2', logsc=TRUE, um=TRUE, noise=FALSE,cex=1)
plotexpmap(sc, 'Foxa2', logsc=TRUE, um=TRUE, noise=TRUE,cex=1)

# Biological noise (on logarithmic scale):
plotexpmap(sc, 'Zfp42', logsc=TRUE, um=TRUE, noise=TRUE,cex=1)
plotexpmap(sc, 'Zfp42', logsc=TRUE, um=TRUE, noise=TRUE,cex=1)

plotexpmap(sc, 'Mki67', logsc=TRUE, um=TRUE, noise=TRUE,cex=0.5)

p1 = plotexpmap(sc, 'Pax6', logsc=TRUE, um=TRUE, noise=FALSE,cex=1)
p2 = plotexpmap(sc, 'Foxa2', logsc=TRUE, um=TRUE, noise=FALSE,cex=1)

p1 / p2


pdfname = paste0(outDir, '/plot_noise_geneExamples.pdf')
pdf(pdfname, width=8, height = 6)

for(n in 1:length(gene_examples))
{
  if(!is.na(match(gene_examples[n], sc@genes))){
    cat(n, '--', gene_examples[n], '\n')
    plotexpmap(sc, gene_examples[n], logsc=TRUE, um=TRUE, noise=FALSE, cex=1)
    plotexpmap(sc, gene_examples[n], logsc=TRUE, um=TRUE, noise=TRUE, cex=1)
  }else{
    cat(n, '--', gene_examples[n], ' missing \n')
  }
  
}

dev.off()

genes <- gene_examples
mm = match(genes, rownames(sc@noise))
genes = genes[!is.na(mm)]

pdfname = paste0(outDir, '/plot_expression_noise_geneExamples_heatmap.pdf')
pdf(pdfname, width=8, height = 12)

ph <- plotmarkergenes(sc,genes=genes,noise=FALSE, )

plotmarkergenes(sc,genes=genes[ph$tree_row$order], noise=TRUE,cluster_rows=FALSE)

dev.off()

##########################################
# Select noisy genes for each time point 
##########################################
mgenes1 <- maxNoisyGenesTB(noise,cl=cl, set=1)
head(mgenes1)

mgenes2 <- maxNoisyGenesTB(noise,cl=cl, set=2)

mgenes3 <- maxNoisyGenesTB(noise,cl=cl, set=3)

save(mgenes1, mgenes2, mgenes3, 
     file = paste0(outDir, '/out_varID2_noise_maxNoisyGenesTB.Rdata'))


plotmarkergenes(sc,genes=head(names(mgenes1),50), noise=TRUE)

# test = noise$epsilon
# ss = apply(test, 1, function(x) sum(x>1, na.rm = TRUE))
# genes = rownames(test)[which(ss>50)]
# mm = match(genes, c(tfs, sps, gene_examples))
# genes = genes[which(!is.na(mm))]

#genes <- c("Lyz1","Agr2","Clca3","Apoa1","Aldob","Clca4","Mki67","Pcna")
#ph <- plotmarkergenes(sc,genes=genes,noise=FALSE)

#genes <- c("Lyz1","Agr2","Clca3","Apoa1","Aldob","Clca4","Mki67","Pcna")
#ph <- plotmarkergenes(sc,genes=genes,noise=FALSE)
#plotmarkergenes(sc,genes=genes[ph$tree_row$order],noise=TRUE,cluster_rows=FALSE)
#fractDotPlot(sc, genes, zsc=TRUE)

# ngenes_1 <- diffNoisyGenesTB(noise, cl, set=c(2), bgr = c(1), no_cores=32)
# ngenes_2 <- diffNoisyGenesTB(noise, cl, set=c(3), bgr = c(1), no_cores=32)
# ngenes_3 <- diffNoisyGenesTB(noise, cl, set=c(3), bgr = c(2), no_cores=32)
# ngenes_4 <- diffNoisyGenesTB(noise, cl, set=c(1), bgr = c(2), no_cores=32)
# ngenes_5 <- diffNoisyGenesTB(noise, cl, set=c(1), bgr = c(3), no_cores=32)
# ngenes_6 <- diffNoisyGenesTB(noise, cl, set=c(2), bgr = c(3), no_cores=32)

#xx = ngenes_1[order(ngenes_1$pvalue), ] 
#maxNoisyGenesTB(noise,cl=cl,set=3)

#plotDiffNoise(ngenes_1, pthr = 10^-40, lthr = 0.5)
#plotexpmap(sc, 'Lef1', logsc=TRUE, um=TRUE, noise=FALSE, cex=1)
#ngenes = ngenes[order(ngenes$pvalue), ]
#head(ngenes, 50)

# cells = colnames(sc@ndata)
# cl_new = cl$partition
# #newOrder = c(1, 2, 4, 5, 6, 3, 7, 8,9, 10, 11)
# cl_new[which(cl$partition == 4)] = 3
# cl_new[which(cl$partition == 5)] = 4
# cl_new[which(cl$partition == 6)] = 5
# cl_new[which(cl$partition == 3)] = 6
# 
# cl_new = cl_new[order(cl_new)]
#load(file = paste0(outDir, '/out_varID2_noise_diffNoisyGenes_allpairComps.Rdata'))



# select_highNoisey_genes = function(ngenes, logfc_cutoff = 1, pval_cutoff = 10^-20)
# {
#   ngenes = ngenes[which(abs(ngenes$log2FC) > logfc_cutoff & ngenes$pvalue < pval_cutoff), ]
#   mm = match(rownames(ngenes), c(tfs, sps, gene_examples))
#   ngenes = ngenes[which(!is.na(mm)), ]
#   return(ngenes)
# }
# 
# ngenes_1 = select_highNoisey_genes(ngenes_1)
# ngenes_2 = select_highNoisey_genes(ngenes_2)
# ngenes_3 = select_highNoisey_genes(ngenes_3)
# ngenes_4 = select_highNoisey_genes(ngenes_4)
# ngenes_5 = select_highNoisey_genes(ngenes_5)
# ngenes_6 = select_highNoisey_genes(ngenes_6)

# genes = unique(c(rownames(ngenes_1), rownames(ngenes_2), rownames(ngenes_3), 
#                  c(rownames(ngenes_4), rownames(ngenes_5), rownames(ngenes_6),
#                  'Pax6', 'Dhrs3', 'Lef1')))

genes = unique(c(names(mgenes1)[which(mgenes1 > 0.7)], 
                 names(mgenes2)[which(mgenes2 > 0.7)], 
                 names(mgenes3)[which(mgenes3 > 0.7)],          
                 'Pax6', 'Dhrs3', 'Lef1'))
mm = match(genes, c(tfs, sps, gene_examples))
genes = genes[which(!is.na(mm))]        
                 
#genes = head(genes, 400)
plotmarkergenes(sc,genes=head(names(mgenes1),50), noise=TRUE)


pdfname = paste0(outDir, '/plot_highNoise_genes_byTimpoint_tops.pdf')
pdf(pdfname, width=8, height = 16)

ph <- plotmarkergenes(sc, genes=genes, 
                      noise=TRUE,
                      cluster_rows=TRUE, 
                      #cells = names(cl_new), 
                      order.cells = TRUE,
                      cluster_set = FALSE,
                      cluster_cols = FALSE,
                      cap = 5, #flo = -3,
                      zsc=FALSE, 
                      logscale = TRUE)

dev.off()


##########################################
# GO term test for each time point 
##########################################
library(enrichplot)
library(clusterProfiler)
library(stringr)
library(org.Mm.eg.db)


index_cluster = 3
mgenes = mgenes3

genes = unique(c(names(mgenes)[which(mgenes > 0.6)]))
ego <-  enrichGO(gene         = genes,
                 #universe     = bgs0.df$ENSEMBL,
                 #OrgDb         = org.Hs.eg.db,
                 OrgDb         = org.Mm.eg.db,
                 #keyType       = 'ENSEMBL',
                 keyType =  "SYMBOL",
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 #pvalueCutoff  = 0.01,
                 #qvalueCutoff  = 0.2,
                 readable=FALSE)
#head(ego)
barplot(ego, showCategory=20) + ggtitle("Go term enrichment")


pdfname = paste0(outDir, '/Goterm_enrichment_cluster', index_cluster, '.pdf')
pdf(pdfname, width=8, height = 12)

barplot(ego, showCategory=20) + ggtitle("Go term enrichment")

dev.off()

write.csv2(ego, file = paste0(outDir, "/GO_term_enrichmenet_highNoisyGenes_cluster", index_cluster, ".csv"), 
          row.names = TRUE)


##########################################
# Visualize gene examples from noisy gene analysis by VarID2
##########################################
# Nodal/Msx1/Neurod1/Ets2/Ttn/Gpc3
# /Id1/Msx1/Smad7/Egr1/
genes = c('Nodal', "Msx1", "Neurod1", 'Ets2', 'Ttn', 'Gpc3', 'Id1', 'Egr1', 'Smad7')

p1 = FeaturePlot(aa, features = genes, reduction = 'wnn.umap')
plot(p1)

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_cluster1.pdf'), 
       height = 21, width = 27)

#Dppa3/Tdrd12/Fkbp6/Mov10l1
genes = c('Dppa3', "Tdrd12", "Fkbp6", 'Mov10l1', 'Piwil2', 'Exd1', 'Dnmt1', 'Dnmt3a', 'Dnmt3b',
          "Tet1", "Tet2", "Tet3")

p1 = FeaturePlot(aa, features = genes, reduction = 'wnn.umap')
plot(p1)

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_cluster1_DNAmethylation.pdf'), 
       height = 21, width = 27)


# Aldh1a3/Cyp26b1/Crabp2/Cyp26a1 ALDH1A1, ALDH1A2
genes = c("Aldh1a3", "Aldh1a1", "Aldh1a2",
          "Cyp26b1", "Cyp26a1", "Crabp2", 
          "Rara", "Rarb", "Rarg")

p1 = FeaturePlot(aa, features = genes, reduction = 'wnn.umap')
plot(p1)

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_cluster_RAresponse.pdf'), 
       height = 16, width = 27)

# Pbld2/Cflar/Zeb2/Cdkn1c/Fgfbp1/Gata3/Fgfbp3/Fgf4/Folr1/Hes1
genes = c("Pbld2", "Cflar", "Zeb2",
          "Cdkn1c", "Fgfbp1", "Gata3", 
          "Fgfbp3", "Fgf4", "Folr1", "Hes1", 'Foxa2')

p1 = FeaturePlot(aa, features = genes, reduction = 'wnn.umap')
plot(p1)

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_cluster_GrowthFactor.pdf'), 
       height = 16, width = 27)

# Tcf15/Zeb2/Cyp26b1/Gbx2/Sox17/Hoxb1/Cdx1/Hoxb2/Pax6/Cyp26a1/Tshz1/Ooep/Folr1/Hes1/Lef1

genes = c("Tcf15", "Zeb2", "Gbx2",
          "Sox17", "Hoxb1", "Cdx1", 
          "Hoxb2", "Ooep", "Lef1")

p1 = FeaturePlot(aa, features = genes, reduction = 'wnn.umap')
plot(p1)

ggsave(filename = paste0(outDir, 
                         '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_cluster_patternSpecification.pdf'), 
       height = 16, width = 27)

# Msx1/Mafb/Sox17/Foxa2/Lhx1/Fst/Nr2f2/Cyp26a1/Fgf3/Hoxb1/Pax6/Nr2f1/Ifitm1/Hoxb2/
# Gdf3/Hes1/Folr1/Zic3/Robo2/Gbx2/Sema3a/Sp8/Hes3/Sfrp1/Dll1
genes = c("Mafb", "Lhx1", "Fst",
          "Nr2f2", "Fgf3", "Nr2f1", 
          "Ifitm1", "Gdf3", "Hes1", 
          "Zic3", "Robo2", "Gbx2", "Sema3a", "Sp8", "Sfrp1", "Dll1")

p1 = FeaturePlot(aa, features = genes, reduction = 'wnn.umap')
plot(p1)

ggsave(filename = paste0(outDir, 
                         '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_cluster_cluster3_patterning.pdf'), 
       height = 21, width = 27)



########################################################
########################################################
# Section III : prepare the scRNA-seq data for OT analysis
# 
########################################################
########################################################
outDir = paste0(resDir, '/mNTs_multiome_RA_searching_earlyBias/',
                'Test_OT/')

system(paste0('mkdir -p ', outDir))


##########################################
# import the RA samples and cleaning
# subseting the cells for OT test
##########################################
aa = readRDS(file = paste0('../results/Rdata/', 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', '_timePoint_',
                           'mNT_scRNAseq', 
                           '_R13547_10x_mNT_20220813', '.rds'))

p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'clusters', label = TRUE, repel = TRUE)

p1 + p2

DimPlot(aa, group.by = 'Phase', label = TRUE, repel = TRUE)

Idents(aa) = aa$condition

set.seed(2023)

# remove the cluster 9 and 10
aa = subset(aa, cells = colnames(aa)[which(as.character(aa$clusters) != '9' 
                                           & as.character(aa$clusters) != '10')])
#aa = subset(aa, cells = colnames(aa)[which(!is.na(match(aa$condition, cc)))])

aa$condition = droplevels(aa$condition)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)

aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

p2 = DimPlot(aa, group.by = 'Phase', label = TRUE, repel = TRUE) 

p1 + p2

saveRDS(aa, file = paste0(outDir, 'scRNAseq_d2_d2.5_d3_d4_d5_RA_cleaned.rds'))

##########################################
#  test metacell to have small clusters for OT
##########################################
## test metacell approach 
## original code from 
## https://gfellerlab.github.io/MetacellAnalysisTutorial/Metacell-
## construction-chapter.html#SuperCell-construction
library(SuperCell)

## subsampling cell at each time point
aa = readRDS(file = paste0(outDir, 'scRNAseq_d2_d2.5_d3_d4_d5_RA_cleaned.rds'))

set.seed(2026)
Idents(aa) = aa$condition
aa = subset(aa, downsample = 2000)


sc_data = aa
annotation_label = 'condition'

#sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize")
sc_data <- FindVariableFeatures(sc_data, nfeatures = 3000)
#sc_data <- ScaleData(sc_data)
#> Centering and scaling data matrix
sc_data <- RunPCA(sc_data, npcs = 100, verbose = F)

sc_data <- RunUMAP(sc_data, reduction = "pca", dims = c(1:30), n.neighbors = 30, verbose = F, 
                   min.dist = 0.1)
UMAPPlot(sc_data, group.by = "condition")

ggsave(filename = paste0(outDir, '/plot_UMAP_for_metacell.pdf'), 
       width = 10, height = 6)

FeaturePlot(sc_data, features = c('Foxa2', 'Pax6'))

ggsave(filename = paste0(outDir, '/plot_UMAP_for_metacell_FoxA2_Pax6.pdf'), 
       width = 14, height = 6)


gamma = 200 # the requested graining level.
k_knn = 20 # the number of neighbors considered to build the knn network.
nb_var_genes = 3000 # number of the top variable genes to use for dimensionality reduction 
nb_pc = 30 # the number of principal components to use.   
cell_condition = sc_data$condition

MC <- SuperCell::SCimplify(Seurat::GetAssayData(sc_data, slot = "data"),  
                           k.knn = k_knn,
                           gamma = gamma,
                           #n.var.genes = nb_var_genes,  
                           cell.split.condition = cell_condition,
                           n.pc = nb_pc,
                           genes.use = Seurat::VariableFeatures(sc_data)
)

MC.GE <- supercell_GE(Seurat::GetAssayData(sc_data, slot = "counts"),
                      MC$membership,
                      mode =  "sum"
)
dim(MC.GE) 


print(annotation_label)

MC$annotation <- supercell_assign(clusters = sc_data@meta.data[, annotation_label], # single-cell annotation
                                  supercell_membership = MC$membership, # single-cell assignment to metacells
                                  method = "absolute"
)

head(MC$annotation)

table(MC$annotation)

pdfname = paste0(outDir, '/plot_metacell_graph_v2.pdf')
pdf(pdfname, width=10, height = 6)

supercell_plot(
  MC$graph.supercells, 
  group = MC$annotation,
  lay.method = 'fr',
  seed = 1, 
  alpha = -pi/2,
  min.cell.size = 20,
  main  = "Metacells colored by condition"
)

dev.off()


## save the metacell into Seurat object
mcs  = MC$membership
aa$metacell = mcs[match(colnames(aa), names(mcs))]

cc = as.character(unique(aa$condition))
for(n in 1:length(cc))
{
  cat(cc[n], '-- nb of metacells ', 
      length(unique(aa$metacell[which(aa$condition == cc[n])])), '\n')
}

mc_size = table(aa$metacell)
mc_small = names(mc_size)[which(mc_size < 50)]

aa$metacell[!is.na(match(aa$metacell, mc_small))] = NA

cc = as.character(unique(aa$condition))
for(n in 1:length(cc))
{
  cat(cc[n], '-- nb of metacells ', 
      length(unique(aa$metacell[which(aa$condition == cc[n] & !is.na(aa$metacell))])), '\n')
}

DimPlot(aa, group.by = "metacell", label = TRUE, repel = TRUE) + NoLegend()

## processing the metacell labels
aa = subset(aa, cells = colnames(aa)[which(!is.na(aa$metacell))])

cc = as.character(unique(aa$condition))
aa$metacell_label = NA
for(c in cc)
{
  # c = "day2_beforeRA"
  kk = which(aa$condition == c)
  index = unique(aa$metacell[kk])
  
  for(n in 1:length(index)) 
  {
    aa$metacell_label[which(aa$condition == c & aa$metacell == index[n])] = paste0(c, '_', n)
  }
  
}

aa$metacell_label = gsub('day2_beforeRA_', 'd2_', aa$metacell_label)
aa$metacell_label = gsub('day2.5_RA_', 'd2.5_', aa$metacell_label)
aa$metacell_label = gsub('day3_RA.rep1_', 'd3_', aa$metacell_label)
aa$metacell_label = gsub('day3.5_RA_', 'd3.5_', aa$metacell_label)
aa$metacell_label = gsub('day4_RA_', 'd4_', aa$metacell_label)
aa$metacell_label = gsub('day5_RA_', 'd5_', aa$metacell_label)

DimPlot(aa, group.by = "metacell_label", label = TRUE, repel = TRUE) + NoLegend()

ggsave(filename = paste0(outDir, '/plot_UMAP_for_metacell_index_forOT.pdf'), 
       width = 12, height = 8)


saveRDS(aa, file = paste0(outDir, 'scRNAseq_RA_d2_d2.5_d3_d4_d5_metacellAnnot_forOTtest.rds'))


##########################################
# save file for OT moscot  
##########################################
library(SeuratDisk)

aa = readRDS(file = paste0(outDir, 'scRNAseq_RA_d2_d2.5_d3_d4_d5_metacellAnnot_forOTtest.rds'))

aa$condition = as.character(aa$condition)

FeaturePlot(aa, features = c("Foxa2", "Pax6"))

DimPlot(aa, group.by = "metacell_label", label = TRUE, repel = TRUE) + NoLegend()

ggsave(filename = paste0(outDir, '/plot_UMAP_for_metacell_index_forOT_beforeMerge.pdf'), 
       width = 12, height = 8)

mm = match(aa$metacell_label, c("d5_1", "d5_2", "d5_3", "d5_5", "d5_6", "d5_7", "d5_9", 
                                      "d5_10", "d5_15", "d5_18"))
aa$metacell_label[!is.na(mm)] = 'd5_NP'

DimPlot(aa, group.by = "metacell_label", label = TRUE, repel = TRUE) + NoLegend()

mm = match(aa$metacell_label, c("d5_4", "d5_8", "d5_11", "d5_12", "d5_13", "d5_14", "d5_16", 
                                "d5_17"))
aa$metacell_label[!is.na(mm)] = 'd5_FP'

DimPlot(aa, group.by = "metacell_label", label = TRUE, repel = TRUE) + NoLegend()


mm = match(aa$metacell_label, c("d4_4", "d4_5", "d4_7", "d4_9"))
aa$metacell_label[!is.na(mm)] = 'd4_NP'

mm = match(aa$metacell_label, c("d4_1", "d4_2", "d4_3", "d4_6", "d4_8", "d4_10"))
aa$metacell_label[!is.na(mm)] = 'd4_FP'

DimPlot(aa, group.by = "metacell_label", label = TRUE, repel = TRUE) + NoLegend()

saveRDS(aa, file = paste0(outDir, 'scRNAseq_RA_d2_d2.5_d3_d4_d5_metacellAnnot_manualPPFP_forOTtest.rds'))


library(SeuratDisk)
aa = readRDS(file = paste0(outDir, 'scRNAseq_RA_d2_d2.5_d3_d4_d5_metacellAnnot_manualPPFP_forOTtest.rds'))

aa$condition = as.character(aa$condition)

mnt = aa
VariableFeatures(mnt) = NULL

DefaultAssay(mnt) = 'RNA'

mnt = DietSeurat(mnt, 
                 counts = TRUE, 
                 data = TRUE,
                 scale.data = FALSE,
                 features = rownames(mnt), 
                 assays = c('RNA'), 
                 dimreducs = c('umap', 'pca'), graphs = NULL, 
                 misc = TRUE
)

DefaultAssay(mnt) = 'RNA'
VariableFeatures(mnt)

Idents(mnt) = mnt$condition

#mnt = subset(mnt, downsample = 1000)

saveFile =  'scRNAseq_RA_d2_d2.5_d3_d4_d5_metacellAnnot_forOTtest_v2.h5Seurat'

SaveH5Seurat(mnt, filename = paste0(outDir, saveFile), 
             overwrite = TRUE)

Convert(paste0(outDir, saveFile), dest = "h5ad", overwrite = TRUE)
