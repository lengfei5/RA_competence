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

levels_sels = c("day0_beforeRA", "day1_beforeRA", "day2_beforeRA",
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA",
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA")

# levels_sels = c("day0_beforeRA", "day1_beforeRA", "day2_beforeRA",
#                 "day2.5_RA", "day3_RA.rep1", "day3.5_RA",
#                 "day2.5_noRA", "day3_noRA",  "day3.5_noRA",
#                 'day4_noRA', 'day4_RA', 'day5_noRA', 'day5_RA')

# levels_sels = c("day2_beforeRA",  
#                 "day2.5_RA", "day3_RA.rep1", "day3.5_RA", 
#                 "day2.5_noRA", "day3_noRA",  "day3.5_noRA")

cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/RA.vs.noRA_firstBifurcation/cell_heterogeneity_start.day0')
system(paste0('mkdir -p ', outDir))


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
require(parallel)

parallel::detectCores()
nb_cores = 64

##########################################
# prepare and process the samples
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                           'cellCycleScoring_annot.v1_', 
                           species, version.analysis, '.rds'))

Idents(aa) = factor(aa$condition, levels = levels)

levels_sels = c("day0_beforeRA", "day1_beforeRA", "day2_beforeRA",
                 "day2.5_RA", "day3_RA.rep1", "day3.5_RA",
                 "day2.5_noRA", "day3_noRA",  "day3.5_noRA",
                 'day4_noRA', 'day4_RA', 'day5_noRA', 'day5_RA')

#levels_sels = c("day2_beforeRA",  
#                 "day2.5_RA", "day3_RA.rep1", "day3.5_RA", 
#                 "day2.5_noRA", "day3_noRA",  "day3.5_noRA")

cols_sel = cols[match(levels_sels, names(cols))]


#Idents(aa) = aa$condition
table(aa$condition)

set.seed(2023)
aa = subset(aa, downsample = 500)

aa = subset(aa, idents = levels_sels)

table(aa$condition)

########################################################
########################################################
# Section : test the heteregeneity quantification in 
# Mohammed et al., 2017
########################################################
########################################################
Test_Heteregeneity_pairwiseDist = FALSE
if(Test_Heteregeneity_pairwiseDist){
  library(scater)
  library(SingleCellExperiment)
  library(scran)
  
  Idents(aa) = factor(aa$condition, levels = levels_sels)
  #cc = unique(aa$condition)
  
  Run_HVGs_perTimePoint = FALSE
  if(Run_HVGs_perTimePoint){
    hete = c()
    nb_features = 1000
    for(n in 1:length(levels_sels))
    {
      # n = 1
      cat(n, ' -- ', levels_sels[n], '\n')
      sce = as.SingleCellExperiment(subset(aa, idents = levels_sels[n]))
      
      dec <- modelGeneVar(sce)
      top.hvgs <- getTopHVGs(dec, n=nb_features)
      sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 30)
      # reducedDimNames(sce)
      #ll.pca = reducedDim(sce, 'PCA')[, c(1:30)]
      ll.pca = logcounts(sce)
      ll.pca = as.matrix(ll.pca[match(top.hvgs, rownames(ll.pca)), ])
      dists = sqrt((1-cor((ll.pca), method = 'pearson'))/2)
        
      #dists = dists[which(dists>0)]
      
      hete = rbind(hete, data.frame(condition = rep(levels_sels[n], length(dists)), dists))
      
    }
  }else{
    
    nb_features = 1000
    sce = as.SingleCellExperiment(aa)
    
    dec <- modelGeneVar(sce)
    top.hvgs <- getTopHVGs(dec, n=nb_features)
    sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 100)
    
    # reducedDimNames(sce)
    
    ll.pca = reducedDim(sce, 'PCA')[, c(1:30)]
    
    hete = c()
    #ll.pca = logcounts(sce)
    #ll.pca = as.matrix(ll.pca[match(top.hvgs, rownames(ll.pca)), ])
    for(n in 1:length(levels_sels))
    {
      # n = 1
      cat(n, ' -- ', levels_sels[n], '\n')
      kk = match(colnames(sce)[which(sce$condition == levels_sels[n])], rownames(ll.pca))
      
      dists =  dist(ll.pca[kk, ], method = "euclidean", diag = FALSE, upper = TRUE) #sqrt((1-cor((ll
      #dists = dists[which(dists>0)]
      
      hete = rbind(hete, data.frame(condition = rep(levels_sels[n], length(dists)), as.numeric(dists)))
      
    }
    
  }
  
  #hete$dists = log10(hete$dists)
  colnames(hete)[2] = 'dists'
  hete = as.data.frame(hete)
  
  ggplot(hete, aes(x=condition, y=dists, fill=condition)) + 
    geom_violin() + 
    scale_fill_manual(values = cols_sel) + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, size = 10, vjust = 0.6),
          axis.text.y = element_text(angle = 0, size = 10)) +
    labs( x = '', y = 'pairwise distance' )
 
  ggsave(filename = paste0(outDir, '/pairwise_distance_', nb_features,
                           'HVGs_allConditions.pdf'), width = 12, height = 8) 
  
  
  
}

########################################################
########################################################
# Section II: VarID2 
# 
########################################################
########################################################
Run_VarID2 = FALSE
if(Run_VarID2){
  x <- aa@assays$RNA@counts
  rownames(x) <- rownames(aa)
  colnames(x) <- colnames(aa)
  
  sc <- SCseq(x)
  
  ## filter the genes 
  sc <- filterdata(sc, mintotal=1000, 
                   FGenes=grep("^Gm\\d",rownames(intestinalData),value=TRUE),
                   CGenes=rownames(x)[grep("^(mt|Rp(l|s)|Gm\\d)",rownames(x))])
  
  expData <- getExpData(sc)
  
  
  tic()
  res   <- pruneKnn(expData, no_cores = 64)
  toc()
  
  saveRDS(res, file = paste0(outDir, '/varID2_pruneKnn.rds'))
  
  tic()
  cl    <- graphCluster(res, pvalue= 0.01, use.leiden = TRUE)
  toc()
  
  probs <- transitionProbs(res,cl)
  
  plotPC(res)
  
  tic()
  sc <- updateSC(sc,res=res,cl=cl)
  toc()
  
  #plotmap(sc,fr=TRUE)
  
  #plotmap(sc,fr=TRUE)
  # Alternatively, a umap representation can be computed for visualization:
  sc <- compumap(sc, min_dist=0.1, n_neighbors = 20)
  
  pdfname = paste0(outDir, '/plot_umap_varID.clusters.pdf')
  pdf(pdfname, width=8, height = 6)
  
  plotmap(sc,um=TRUE)
  
  dev.off()
  
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
  noise <- compTBNoise(res,x,pvalue=0.01, no_cores=nb_cores) # take 40 minutes  
  toc()
  
  #noise <- compTBNoise(res,expData,pvalue=0.01,no_cores=5) 
  sc <- updateSC(sc,res=res,cl=cl,noise=noise)
  
  saveRDS(sc, file = paste0(outDir, '/varID2_noise_v2.rds'))
  save(sc, noise, cl, res, file = paste0(outDir, '/out_varID2_noise_v2.Rdata'))
  
  ##########################################
  # To further investigate transcriptome variability and related quantities, 
  # the function quantKnn allows computation of average noise levels across cells, 
  # cell-to-cell transcriptome correlation, and total UMI counts.
  ##########################################
  load(file = paste0(outDir, '/out_varID2_noise_v2.Rdata'))
  #load(file =  paste0(outDir, '/out_varID2_noise_diffNoisyGenes_v2.Rdata'))
  
  parallel::detectCores()
  
  tic()
  qn <- quantKnn(res, noise, sc, pvalue = 0.01, minN = 5, no_cores = 64)
  toc()
  
  save(sc, noise, cl, res, qn, file = paste0(outDir, '/out_varID2_noise_quntKnn_v2.Rdata'))
  
  
  load(file =  paste0(outDir, '/out_varID2_noise_quntKnn_v2.Rdata'))
  #qn = readRDS(paste0(outDir, '/varID2_pruneKnn.rds'))
  
  cells = colnames(sc@ndata)
  cl_new = cl$partition
  cl_new = data.frame(cluster = cl_new[match(cells, names(cl_new))], 
                      condition = aa$condition[match(cells, colnames(aa))])
  
  #cl_new = data.frome(local.corr = qn$local.corr)
  #sc <- compumap(sc, min_dist=0.1, n_neighbors =20)
  xx = qn$local.corr
  cl_new$local.corr = xx[match(rownames(cl_new), names(xx))]
  
  xx = qn$noise.av
  cl_new$noise.av = xx[match(rownames(cl_new), names(xx))]
  
  cl_new$dist = sqrt((1 - cl_new$local.corr)/2.0)
  ggplot(cl_new, aes(x=condition, y=dist, fill=condition)) + 
    geom_boxplot()
  
  ggplot(cl_new, aes(x=condition, y=local.corr, fill=condition)) + 
    geom_boxplot()
  
  ggplot(cl_new, aes(x=condition, y=noise.av, fill=condition)) + 
    geom_boxplot()
  
  
  StemCluster <- 4
  #plotQuantMap(qn,"noise.av",sc,um=TRUE,ceil=.6,cex=1)
  
  plotQuantMap(qn,"noise.av",sc,box=TRUE,cluster=StemCluster,  set = c(1, 2, 4, 6, 3, 7, 8, 9, 10, 11))
  
  #plotQuantMap(qn,"local.corr",sc,um=TRUE,logsc=TRUE,cex=1)
  
  plotQuantMap(qn,"umi",sc,box=TRUE,logsc=TRUE,cluster=StemCluster)
  
  pdfname = paste0(outDir, '/boxplot_cell_cell_correlation_neighborhood.pdf')
  pdf(pdfname, width=7, height = 5)
  
  plotQuantMap(qn,"local.corr",sc,box=TRUE,logsc=TRUE,cluster=StemCluster, set = c(4, 2, 1, 3, 7, 10, 6, 5, 9))
  
  dev.off()
  
  
}


########################################################
########################################################
# Section II : 
# ## reload the results and analyze gene level noise
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