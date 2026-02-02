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

FeaturePlot(aa, features = c('Pax6', 'Foxa2'), reduction = 'wnn.umap')

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_wnnUMAP_featurePlot_Pax6_Foxa2.pdf'), 
       height = 6, width = 17)

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
library(RColorBrewer)
require(Matrix)
library(RaceID)
require(tictoc)

cc = c('day2_beforeRA', 'day2.5_RA', 'day3_RA.rep1')

outDir = paste0(resDir, '/mNTs_multiome_RA_searching_earlyBias/gene_Noise_day2_2.5_3_scRNAseq')
system(paste0('mkdir -p ', outDir))

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
aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'RNA_snn_res.0.6', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'RNA_snn_res.0.7', raster=FALSE)

p1 + p2

aa$clusters = aa$RNA_snn_res.0.7
aa$clusters[which(aa$clusters == '7')] = '2'
aa$clusters = droplevels(aa$clusters)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)

p1 + p2


markers = FindMarkers(aa, ident.1 = c('2'), ident.2 = c('5'))

all.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.4)
all.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10

DoHeatmap(aa, features = top10$gene) + NoLegend()

ggsave(filename = paste0(outDir, '/mNTs_scRNA_clusters_heatmap_markerGenes.pdf'), 
       width = 14, height = 30)

all.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20



##########################################
# gene noise with VarID2
##########################################
aa = readRDS(file = paste0(outDir, 'scRNAseq_d2_d2.5_d3_cleaned.rds'))
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
for(gamma in c(2.0, 2.5, 3, 3.5, 4)) # search for the optimal value of gamma
{
  cat('gamma -- ', gamma, '\n')
  
  load(file = paste0(outDir, '/out_varID2_scSeq_pruneKnn_all.Rdata'))
  
  #nn <- inspectKNN(20,expData,res,cl,object=sc,pvalue=0.01,plotSymbol=TRUE,um=TRUE,cex=1)
  #head(nn$pv.neighbours)
  #head(nn$expr.neighbours)
  x <- getFilteredCounts(sc, minexpr=1, minnumber=10)
  
  tic()
  noise <- compTBNoise(res, x, gamma = gamma, 
                       pvalue=0.01, no_cores=64) # take 40 minutes  
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
gamma = 4
load(file = paste0(outDir, '/out_varID2_noise_gamma_', gamma, '_all.Rdata'))

pdfname = paste0(outDir, '/plot_umap_varID_conditions.pdf')
pdf(pdfname, width=8, height = 6)

plotmap(sc, um=TRUE, tp = 1, cex = 0.3)

dev.off()

#tic()
#cl    <- graphCluster(res, pvalue= 0.01, leiden.resolution = 1.5)
#toc()

#sc <- updateSC(sc,res=res,cl=cl)
#plotmap(sc, um=TRUE, tp = 1, cex = 0.3)
#sc <- updateSC(sc,res=res,cl=cl,noise=noise)

# Gene expression (on logarithmic scale):
plotexpmap(sc, 'Pax6', logsc=TRUE, um=TRUE, cex=1, noise = TRUE)

plotexpmap(sc, 'Foxa2', logsc=TRUE, um=TRUE, noise=TRUE,cex=1)

# Biological noise (on logarithmic scale):
plotexpmap(sc, 'Zfp42', logsc=TRUE, um=TRUE, noise=TRUE,cex=0.5)

plotexpmap(sc, 'Mki67', logsc=TRUE, um=TRUE, noise=TRUE,cex=0.5)


p1 = plotexpmap(sc, 'Foxa2', logsc=TRUE, um=TRUE, noise=TRUE,cex=0.5)
p2 = plotexpmap(sc, 'Foxa2', logsc=TRUE, um=TRUE, noise=FALSE,cex=0.5)

p1 + p2

pdfname = paste0(outDir, '/plot_noise_geneExamples.pdf')
pdf(pdfname, width=8, height = 6)

for(n in 1:length(gene_examples))
{
  if(!is.na(match(gene_examples[n], sc@genes))){
    cat(n, '--', gene_examples[n], '\n')
    plotexpmap(sc, gene_examples[n], logsc=TRUE, um=TRUE, noise=TRUE, cex=0.5)
  }else{
    cat(n, '--', gene_examples[n], ' missing \n')
  }
  
}

dev.off()

##########################################
# how to select noisy genes
##########################################
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

ngenes_1 <- diffNoisyGenesTB(noise, cl, set=c(2), bgr = c(1), no_cores=32)
ngenes_2 <- diffNoisyGenesTB(noise, cl, set=c(3), bgr = c(1), no_cores=32)
ngenes_3 <- diffNoisyGenesTB(noise, cl, set=c(3), bgr = c(2), no_cores=32)
ngenes_4 <- diffNoisyGenesTB(noise, cl, set=c(1), bgr = c(2), no_cores=32)
ngenes_5 <- diffNoisyGenesTB(noise, cl, set=c(1), bgr = c(3), no_cores=32)
ngenes_6 <- diffNoisyGenesTB(noise, cl, set=c(2), bgr = c(3), no_cores=32)

save(ngenes_1, ngenes_2, ngenes_3, ngenes_4, ngenes_5, ngenes_6,
     file = paste0(outDir, '/out_varID2_noise_diffNoisyGenes_allpairComps.Rdata'))

#xx = ngenes_1[order(ngenes_1$pvalue), ] 
#maxNoisyGenesTB(noise,cl=cl,set=3)

plotDiffNoise(ngenes_1, pthr = 10^-40, lthr = 0.5)

plotexpmap(sc, 'Lef1', logsc=TRUE, um=TRUE, noise=FALSE, cex=1)
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
load(file = paste0(outDir, '/out_varID2_noise_diffNoisyGenes_allpairComps.Rdata'))

select_highNoisey_genes = function(ngenes, logfc_cutoff = 1, pval_cutoff = 10^-20)
{
  ngenes = ngenes[which(abs(ngenes$log2FC) > logfc_cutoff & ngenes$pvalue < pval_cutoff), ]
  mm = match(rownames(ngenes), c(tfs, sps, gene_examples))
  ngenes = ngenes[which(!is.na(mm)), ]
  return(ngenes)
}

ngenes_1 = select_highNoisey_genes(ngenes_1)
ngenes_2 = select_highNoisey_genes(ngenes_2)
ngenes_3 = select_highNoisey_genes(ngenes_3)
ngenes_4 = select_highNoisey_genes(ngenes_4)
ngenes_5 = select_highNoisey_genes(ngenes_5)
ngenes_6 = select_highNoisey_genes(ngenes_6)


genes = unique(c(rownames(ngenes_1), rownames(ngenes_2), rownames(ngenes_3), 
                 c(rownames(ngenes_4), rownames(ngenes_5), rownames(ngenes_6),
                 'Pax6', 'Dhrs3', 'Lef1')))

#genes = head(genes, 400)

pdfname = paste0(outDir, '/plot_highNoise_genes_conditionsComparisions_tops.pdf')
pdf(pdfname, width=8, height = 20)

ph <- plotmarkergenes(sc, genes=genes, noise=TRUE,
                      cluster_rows=TRUE, 
                      #cells = names(cl_new), 
                      order.cells = TRUE,
                      cluster_set = FALSE,
                      cluster_cols = FALSE,
                      cap = 5, #flo = -3,
                      zsc=TRUE, 
                      logscale = TRUE)

dev.off()


# genes = rownames(ngenes)[which(abs(ngenes$log2FC)>1 | ngenes$pvalue<10^-10)]
# genes <- genes[which(!is.na(match(genes, unique(c(tfs, sps)))))]
# 
# pdfname = paste0(outDir, '/plot_noise_markGenes_all_v3.pdf')
# pdf(pdfname, width=8, height = 24)
# 
# ph <- plotmarkergenes(sc, genes=genes, noise=TRUE,
#                       cluster_rows=TRUE, 
#                       cells = names(cl_new), order.cells = TRUE,
#                       cluster_set = FALSE,
#                       cluster_cols = FALSE,
#                       cap = 5, #flo = -3,
#                       zsc=TRUE, 
#                       logscale = TRUE)
# dev.off()


# ##########################################
# # To further investigate transcriptome variability and related quantities, 
# # the function quantKnn allows computation of average noise levels across cells, 
# # cell-to-cell transcriptome correlation, and total UMI counts.
# ##########################################
# library(parallel)
# load(file =  paste0(outDir, '/out_varID2_noise_diffNoisyGenes_v2.Rdata'))
# 
# parallel::detectCores()
# 
# tic()
# qn <- quantKnn(res, noise, sc, pvalue = 0.01, minN = 5, no_cores = 64)
# toc()
# 
# 
# #sc <- compumap(sc, min_dist=0.1, n_neighbors =20)
# 
# StemCluster <- 1
# #plotQuantMap(qn,"noise.av",sc,um=TRUE,ceil=.6,cex=1)
# 
# plotQuantMap(qn,"noise.av",sc,box=TRUE,cluster=StemCluster,  set = c(1, 2, 4, 6, 3, 7, 8, 9, 10, 11))
# 
# #plotQuantMap(qn,"local.corr",sc,um=TRUE,logsc=TRUE,cex=1)
# 
# plotQuantMap(qn,"umi",sc,box=TRUE,logsc=TRUE,cluster=StemCluster)
# 
# pdfname = paste0(outDir, '/boxplot_cell_cell_correlation_neighborhood.pdf')
# pdf(pdfname, width=7, height = 5)
# plotQuantMap(qn,"local.corr",sc,box=TRUE,logsc=TRUE,cluster=StemCluster, set = c(1, 2, 4,5, 6, 3, 7, 8, 9, 10, 11))
# 
# dev.off()
# 
