##########################################################################
##########################################################################
# Project: RA competence  
# Script purpose: search for genes with high noise at the symmetry breaking point
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed May  3 10:29:39 2023
##########################################################################
##########################################################################
names(cols) = levels

levels_sels = c("day2_beforeRA",
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/RA_symetryBreaking/gene_Noise/')
system(paste0('mkdir -p ', outDir))

version.analysis = paste0('_R13547_10x_mNT_20220813', '_ES.beforeRA.and.RA')

#library(slingshot, quietly = FALSE)
#library(destiny, quietly = TRUE)
#library(mclust, quietly = TRUE)
#library(scater)
#library(SingleCellExperiment)
#library(scran)
library(RColorBrewer)
require(Matrix)
library(RaceID)
require(tictoc)

########################################################
########################################################
# Section I:
# 
########################################################
########################################################
aa = readRDS(file = 
               paste0('../results/scRNAseq_R13547_10x_mNT_20220813/RA_symetryBreaking/branching_genes_BGP',
                           '/RA_d2.beforeRA_no.day6_noNeurons.rds'))
Idents(aa) = aa$condition

set.seed(2023)
aa = subset(aa, downsample = 1000)

x <- aa@assays$RNA@counts
rownames(x) <- rownames(aa)
colnames(x) <- colnames(aa)

sc <- SCseq(x)

## filter the genes 
sc <- filterdata(sc, mintotal=1000, 
                 FGenes=grep("^Gm\\d",rownames(intestinalData),value=TRUE),
                 CGenes=rownames(x)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(x))])

expData <- getExpData(sc)

parallel::detectCores()
tic()
res   <- pruneKnn(expData, no_cores = 64)
toc()


saveRDS(res, file = paste0(outDir, '/RA_d2.beforeRA_no.day6_noNeurons_varID2_pruneKnn.rds'))

tic()
cl    <- graphCluster(res, pvalue= 0.01)
toc()

probs <- transitionProbs(res,cl)

plotPC(res)

tic()
sc <- updateSC(sc,res=res,cl=cl)
toc()

#plotmap(sc,fr=TRUE)

#plotmap(sc,fr=TRUE)
# Alternatively, a umap representation can be computed for visualization:
sc <- compumap(sc, min_dist=0.1, n_neighbors =20)
plotmap(sc,um=TRUE)

pdfname = paste0(outDir, '/plot_transition_probability_v2.pdf')
pdf(pdfname, width=10, height = 8)

plotTrProbs(sc,probs,um=TRUE)
#ggsave(filename = paste0(outDir, '/plot_transition_probability.pdf'), width = 10, height = 8)
dev.off()

## compute noise from corrected variance
#nn <- inspectKNN(20,expData,res,cl,object=sc,pvalue=0.01,plotSymbol=TRUE,um=TRUE,cex=1)
#head(nn$pv.neighbours)
#head(nn$expr.neighbours)

x <- getFilteredCounts(sc, minexpr=5, minnumber=20)
tic()
noise <- compTBNoise(res,x,pvalue=0.01, no_cores=64) # take 40 minutes  
toc()

#noise <- compTBNoise(res,expData,pvalue=0.01,no_cores=5) 
sc <- updateSC(sc,res=res,cl=cl,noise=noise)

saveRDS(sc, file = paste0(outDir, '/RA_d2.beforeRA_no.day6_noNeurons_varID2_noise_v2.rds'))
save(sc, noise, cl, res, file = paste0(outDir, '/out_varID2_noise_v2.Rdata'))

########################################################
########################################################
# Section II : 
# ## reload the results and make analysis
########################################################
########################################################
sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

load(file = paste0(outDir, '/out_varID2_noise_v2.Rdata'))

pdfname = paste0(outDir, '/plot_umap_varID.clusters.pdf')
pdf(pdfname, width=8, height = 6)
plotmap(sc, um=TRUE, tp = 1, cex = 0.3)

dev.off()

plotUMINoise(sc, noise,log.scale=TRUE)
#sc <- updateSC(sc,res=res,cl=cl,noise=noise)

# Gene expression (on logarithmic scale):
gg = 'Pax6'
plotexpmap(sc, gg, logsc=TRUE,um=TRUE,cex=1)

# Biological noise (on logarithmic scale):
plotexpmap(sc, 'Zfp42', logsc=TRUE, um=TRUE, noise=TRUE,cex=0.5)

plotexpmap(sc, 'Foxa2', logsc=TRUE, um=TRUE, noise=TRUE,cex=0.5)
plotexpmap(sc, 'Foxa2', logsc=TRUE, um=TRUE, noise=FALSE,cex=0.5)

p1 + p2

pdfname = paste0(outDir, '/plot_noise_example_Pax6_expression.pdf')
pdf(pdfname, width=8, height = 6)
plotexpmap(sc, 'Pax6', logsc=TRUE, um=TRUE, noise=FALSE,cex=0.5)
dev.off()

pdfname = paste0(outDir, '/plot_noise_example_Pax6.pdf')
pdf(pdfname, width=8, height = 6)
plotexpmap(sc, 'Pax6', logsc=TRUE, um=TRUE, noise=TRUE,cex=0.5)
dev.off()

#genes <- c("Lyz1","Agr2","Clca3","Apoa1","Aldob","Clca4","Mki67","Pcna")
#ph <- plotmarkergenes(sc,genes=genes,noise=FALSE)
#plotmarkergenes(sc,genes=genes[ph$tree_row$order],noise=TRUE,cluster_rows=FALSE)
#fractDotPlot(sc, genes, zsc=TRUE)

ngenes <- diffNoisyGenesTB(noise, cl,set=1, no_cores=32)

ngenes = ngenes[order(-abs(ngenes$log2FC)), ]

save(sc, noise, cl, res, ngenes, file = paste0(outDir, '/out_varID2_noise_diffNoisyGenes_v2.Rdata'))

head(ngenes, 50)

#ngenes = ngenes[order(ngenes$pvalue), ]
#head(ngenes, 50)

cells = colnames(sc@ndata)
cl_new = cl$partition
#newOrder = c(1, 2, 4, 5, 6, 3, 7, 8,9, 10, 11)
cl_new[which(cl$partition == 4)] = 3
cl_new[which(cl$partition == 5)] = 4
cl_new[which(cl$partition == 6)] = 5
cl_new[which(cl$partition == 3)] = 6

cl_new = cl_new[order(cl_new)]

#sc <- updateSC(sc,cl=cl_new)
#plotmap(sc,um=TRUE)

genes <- rownames(ngenes)[which(!is.na(match(rownames(ngenes), unique(c(tfs, sps)))))]
genes = head(genes, 70)

pdfname = paste0(outDir, '/plot_noise_markGenes_ntop.100_v3.pdf')
pdf(pdfname, width=8, height = 12)

ph <- plotmarkergenes(sc, genes=genes, noise=TRUE,
                      cluster_rows=TRUE, 
                      cells = names(cl_new), order.cells = TRUE,
                      cluster_set = FALSE,
                      cluster_cols = FALSE,
                      cap = 5, #flo = -3,
                      zsc=TRUE, 
                      logscale = TRUE)

dev.off()


genes = rownames(ngenes)[which(abs(ngenes$log2FC)>1 | ngenes$pvalue<10^-10)]
genes <- genes[which(!is.na(match(genes, unique(c(tfs, sps)))))]


pdfname = paste0(outDir, '/plot_noise_markGenes_all_v3.pdf')
pdf(pdfname, width=8, height = 24)

ph <- plotmarkergenes(sc, genes=genes, noise=TRUE,
                      cluster_rows=TRUE, 
                      cells = names(cl_new), order.cells = TRUE,
                      cluster_set = FALSE,
                      cluster_cols = FALSE,
                      cap = 5, #flo = -3,
                      zsc=TRUE, 
                      logscale = TRUE)
dev.off()


##########################################
# To further investigate transcriptome variability and related quantities, 
# the function quantKnn allows computation of average noise levels across cells, 
# cell-to-cell transcriptome correlation, and total UMI counts.
##########################################
library(parallel)
load(file =  paste0(outDir, '/out_varID2_noise_diffNoisyGenes_v2.Rdata'))

parallel::detectCores()

tic()
qn <- quantKnn(res, noise, sc, pvalue = 0.01, minN = 5, no_cores = 64)
toc()


#sc <- compumap(sc, min_dist=0.1, n_neighbors =20)

StemCluster <- 1
#plotQuantMap(qn,"noise.av",sc,um=TRUE,ceil=.6,cex=1)

plotQuantMap(qn,"noise.av",sc,box=TRUE,cluster=StemCluster,  set = c(1, 2, 4, 6, 3, 7, 8, 9, 10, 11))

#plotQuantMap(qn,"local.corr",sc,um=TRUE,logsc=TRUE,cex=1)

plotQuantMap(qn,"umi",sc,box=TRUE,logsc=TRUE,cluster=StemCluster)

pdfname = paste0(outDir, '/boxplot_cell_cell_correlation_neighborhood.pdf')
pdf(pdfname, width=7, height = 5)
plotQuantMap(qn,"local.corr",sc,box=TRUE,logsc=TRUE,cluster=StemCluster, set = c(1, 2, 4,5, 6, 3, 7, 8, 9, 10, 11))

dev.off()

##########################################
# test VarID2 (not used anymore) 
# original code from https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html
##########################################
## library(devtools)
## install_github("dgrun/RaceID3_StemID2_package")
# vignette("RaceID")
sc <- SCseq(intestinalData)
sc <- filterdata(sc, mintotal=1000, 
                 FGenes=grep("^Gm\\d",rownames(intestinalData),value=TRUE),
                 CGenes=grep("^(mt|Rp(l|s))",rownames(intestinalData),value=TRUE)
                 )

expData  <- getExpData(sc)
res      <- pruneKnn(expData,no_cores=1)

plotRegNB(expData,res,"(Intercept)")
plotRegNB(expData,res,"beta")
plotRegNB(expData,res,"theta")

plotPearsonRes(res,log=TRUE,xlim=c(-.1,.2))

plotPC(res)

cl <- graphCluster(res, pvalue=0.01)

#install.packages("reticulate")
# library(reticulate)
# conda_list()
# reticulate::py_config()
# conda_list()[[1]][6] %>% 
#   use_condaenv(required = TRUE)

#use_condaenv("cellrank") 
#reticulate::use_python("/groups/tanaka/People/current/jiwang/local/anaconda3/envs/cellrank.py3.9/bin/python", 
#                       required=TRUE)

#confirm that leiden, igraph, and python are available (should return TRUE).
#reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
#reticulate::py_available()

cl <- graphCluster(res,pvalue=0.01,use.leiden=TRUE,leiden.resolution=1.5)
sc <- updateSC(sc,res=res,cl=cl)
plotmap(sc,fr=TRUE)

# After computing a t-SNE map, the clustering can be also be visualized in t-SNE space:
sc <- comptsne(sc,perplexity=50)
plotmap(sc)

# Alternatively, a umap representation can be computed for visualization:
sc <- compumap(sc,min_dist=0.5)
plotmap(sc,um=TRUE)

probs <- transitionProbs(res,cl,pvalue=0.01)
plotTrProbs(sc,probs,um=TRUE)

nn <- inspectKNN(20,expData,res,cl,object=sc,pvalue=0.01,plotSymbol=TRUE,um=TRUE,cex=1)

head(nn$pv.neighbours)
head(nn$expr.neighbours)

nn <- inspectKNN(20,expData,res,cl,object=sc,pvalue=0.01,plotSymbol=FALSE)

# Alternatively, the coefficient of variation (CV) can be plotted with the same model fits:
nn <- inspectKNN(20,expData,res,cl,object=sc,pvalue=0.01,plotSymbol=FALSE,cv=TRUE)

x <- getFilteredCounts(sc,minexpr=5,minnumber=20)
noise <- compTBNoise(res,x,pvalue=0.01,gamma = 0.5,no_cores=1) 

plotUMINoise(sc,noise,log.scale=TRUE)

sc <- updateSC(sc,res=res,cl=cl,noise=noise)

# Gene expression (on logarithmic scale):
plotexpmap(sc,"Clca4",logsc=TRUE,um=TRUE,cex=1)

# Biological noise (on logarithmic scale):
plotexpmap(sc,"Clca4",logsc=TRUE,um=TRUE,noise=TRUE,cex=1)

plotExpNoise("Clca4",sc,noise,norm=TRUE,log="xy")

genes <- c("Lyz1","Agr2","Clca3","Apoa1","Aldob","Clca4","Mki67","Pcna")
ph <- plotmarkergenes(sc,genes=genes,noise=FALSE)

plotmarkergenes(sc,genes=genes[ph$tree_row$order],noise=TRUE,cluster_rows=FALSE)

fractDotPlot(sc, genes, zsc=TRUE)
ngenes <- diffNoisyGenesTB(noise,cl,set=1,no_cores=1)
head(ngenes)

genes <- head(rownames(ngenes),50)
ph <- plotmarkergenes(sc,genes=genes,noise=TRUE,cluster_rows=FALSE,zsc=TRUE)


qn <- quantKnn(res, noise, sc, pvalue = 0.01, minN = 5, no_cores = 1)
StemCluster <- 2
plotQuantMap(qn,"noise.av",sc,um=TRUE,ceil=.6,cex=1)

plotQuantMap(qn,"noise.av",sc,box=TRUE,cluster=StemCluster)
plotQuantMap(qn,"local.corr",sc,um=TRUE,logsc=TRUE,cex=1)

plotQuantMap(qn,"local.corr",sc,box=TRUE,logsc=TRUE,cluster=StemCluster)