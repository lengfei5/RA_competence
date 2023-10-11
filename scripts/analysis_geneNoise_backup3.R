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

cc = c('day2_beforeRA', 'day2.5_RA', 'day3_RA.rep1', 
       'day1_beforeRA', 'day2.5_noRA')

outDir = paste0(resDir, '/RA_symetryBreaking/gene_Noise/RAday2_2.5_3_bgs_noRA.day2.5_day1')
system(paste0('mkdir -p ', outDir))

########################################################
########################################################
# Section I:
# 
########################################################
########################################################
#aa = readRDS(file = 
#               paste0('../results/scRNAseq_R13547_10x_mNT_20220813/RA_symetryBreaking/branching_genes_BGP',
#                           '/RA_d2.beforeRA_no.day6_noNeurons.rds'))
#aa = readRDS(file = paste0(RdataDir, 
#                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
#                           'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', '_timePoint_',
#                           species, version.analysis, '.rds'))
aa =  readRDS(file = paste0('../results/Rdata/',  
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', 'mNT_scRNAseq',
                            '_R13547_10x_mNT_20220813', '.rds'))

# bb  = aa
p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'clusters', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(outDir, '/UMAP_conditions_clusters.pdf'), 
       width = 16, height = 6)

Idents(aa) = aa$condition

cat(which(rownames(aa) == 'Dhrs3'), '\n')

set.seed(2023)

# remove the cluster 9 and 10
#aa = subset(aa, cells = colnames(aa)[which(!is.na(match(aa$condition, cc))
#                                           & as.character(aa$clusters) != '9' 
#                                           & as.character(aa$clusters) != '10')])

aa = subset(aa, cells = colnames(aa)[which(!is.na(match(aa$condition, cc)))])
Idents(aa) = factor(aa$condition)
#aa = subset(aa, idents = levels_sels)


## remove cell cycle related genes
scaledMatrix = GetAssayData(aa, slot = c("scale.data"))

diff <- scater::getVarianceExplained(scaledMatrix, data.frame(phase = aa$Phase))
diff = data.frame(diff, gene = rownames(diff))
diff = diff[order(-diff$phase), ]

hist(diff$phase, breaks = 100); abline(v = c(1:5), col = 'red')
rm(scaledMatrix)

genes_discard = diff$gene[which(diff$phase>5)]
cat(length(genes_discard), 'genes to discard \n')

aa = subset(aa, features = setdiff(rownames(aa), genes_discard))

Idents(aa) = aa$condition
aa = subset(aa, downsample = 200)

cat(which(rownames(aa) == 'Dhrs3'), '\n')

x <- aa@assays$RNA@counts
rownames(x) <- rownames(aa)
colnames(x) <- colnames(aa)
sc <- SCseq(x)
rm(x)
## filter the genes 
sc <- filterdata(sc, 
                 mintotal = 3000,
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

saveRDS(res, file = paste0(outDir, '/RA_d2.beforeRA_no.day6_noNeurons_varID2_pruneKnn.rds'))

res = readRDS(file = paste0(outDir, '/RA_d2.beforeRA_no.day6_noNeurons_varID2_pruneKnn.rds'))

#plotRegNB(expData, res,"(Intercept)")
#plotRegNB(expData,res,"beta")
#plotRegNB(expData,res,"theta")
#plotPearsonRes(res,log=TRUE,xlim=c(-.1,.2))

tic()
cl    <- graphCluster(res, pvalue= 0.01, leiden.resolution = 0.3)
toc()

test = cl$partition
mm = match(names(test), colnames(aa))

newtest = test
newtest[which(aa$condition[mm] == 'day2_beforeRA')] = 1
newtest[which(aa$condition[mm] == 'day2.5_RA')] = 2
newtest[which(aa$condition[mm] == 'day3_RA.rep1')] = 3
newtest[which(aa$condition[mm] == 'day1_beforeRA')] = 4
newtest[which(aa$condition[mm] == 'day2.5_noRA')] = 5

cl$partition = newtest

probs <- transitionProbs(res,cl)

plotPC(res)

tic()
sc <- updateSC(sc,res=res,cl=cl)
toc()

#plotmap(sc,fr=TRUE)

#plotmap(sc,fr=TRUE)
# Alternatively, a umap representation can be computed for visualization:
sc <- compumap(sc, min_dist= 0.1, n_neighbors =30)

pdfname = paste0(outDir, '/plot_umap.pdf')
pdf(pdfname, width=10, height = 8)

plotmap(sc,um=TRUE)

dev.off()

pdfname = paste0(outDir, '/plot_transition_probability.pdf')
pdf(pdfname, width=10, height = 8)

plotTrProbs(sc,probs,um=TRUE)

dev.off()

save(sc, cl, res, 
     file = paste0(outDir, '/out_varID2_scSeq_pruneKnn.Rdata'))

##########################################
# ## compute noise from corrected variance
# https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html#varid2
# The most important argument of the compTBNoise is the gamma parameter 
# which determines the strength of the prior for the suppression of spurious inflation of noise levels
# at low expression
# If a positive correlation is observed, gamma should be increased in order to weaken the prior. 
# If the correlation is negative, gamma should be decreased in order to increase the strength of the prior.
##########################################
for(gamma in c(1.5, 2.0, 2.5, 3, 3.5,4)) # search for the optimal value of gamma
{
  cat('gamma -- ', gamma, '\n')
  
  load(file = paste0(outDir, '/out_varID2_scSeq_pruneKnn.Rdata'))
  
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
       file = paste0(outDir, '/out_varID2_noise_gamma_', gamma, '.Rdata'))
  
}

########################################################
########################################################
# Section II : 
# ## reload the results and make analysis
########################################################
########################################################
gamma = 4
load( file = paste0(outDir, '/out_varID2_noise_gamma_', gamma, '.Rdata'))


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
gg = 'Pax6'
plotexpmap(sc, gg, logsc=TRUE, um=TRUE, cex=1)

# Biological noise (on logarithmic scale):
plotexpmap(sc, 'Zfp42', logsc=TRUE, um=TRUE, noise=TRUE,cex=0.5)

plotexpmap(sc, 'Foxa2', logsc=TRUE, um=TRUE, noise=TRUE,cex=0.5)
plotexpmap(sc, 'Foxa2', logsc=TRUE, um=TRUE, noise=FALSE,cex=0.5)


pdfname = paste0(outDir, '/plot_noise_geneExamples.pdf')
pdf(pdfname, width=8, height = 6)

genes_noise = rownames(noise$epsilon)
for(n in 1:length(gene_examples))
{
  if(!is.na(match(gene_examples[n], genes_noise))){
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
plotmap(sc, um=TRUE, tp = 1, cex = 0.3)

ngenes_1 <- diffNoisyGenesTB(noise, cl, set=c(1), bgr = c(4), no_cores=32)
ngenes_2 <- diffNoisyGenesTB(noise, cl, set=c(2), bgr = c(5), no_cores=32)
ngenes_3 <- diffNoisyGenesTB(noise, cl, set=c(3), bgr = c(5), no_cores=32)
ngenes_4 <- diffNoisyGenesTB(noise, cl, set=c(2), bgr = c(1), no_cores=32)
ngenes_5 <- diffNoisyGenesTB(noise, cl, set=c(3), bgr = c(2), no_cores=32)

save(ngenes_1, ngenes_2, ngenes_3,
     ngenes_4, ngenes_5,
     file = paste0(outDir, '/out_varID2_noise_diffNoisyGenes_bgs.Rdata'))

#xx = ngenes_1[order(ngenes_1$pvalue), ] 
#maxNoisyGenesTB(noise,cl=cl,set=3)

pdfname = paste0(outDir, '/plot_diffNoise_comparisons.pdf')
pdf(pdfname, width=12, height = 8)

plotDiffNoise(ngenes_1, pthr = 10^-10, lthr = 1, set.name = 'd2', bgr.name = 'd1')
plotDiffNoise(ngenes_2, pthr = 10^-10, lthr = 1, set.name = 'RA_d2.5', bgr.name = 'noRA.d2.5')
plotDiffNoise(ngenes_3, pthr = 10^-10, lthr = 1, set.name = 'RA_d3', bgr.name = 'noRA.d2.5')
plotDiffNoise(ngenes_4, pthr = 10^-10, lthr = 1, set.name = 'RA_d2.5', bgr.name = 'd2')
plotDiffNoise(ngenes_5, pthr = 10^-10, lthr = 1, set.name = 'RA_d3', bgr.name = 'd2')

dev.off()

# ngenes = ngenes[order(-abs(ngenes$log2FC)), ]
# ngenes = ngenes[order(ngenes$pvalue), ]
# 
# save(sc, noise, cl, res, ngenes, 
#      file = paste0(outDir, '/out_varID2_noise_diffNoisyGenes_v3.Rdata'))
# 
# load(file = paste0(outDir, '/out_varID2_noise_diffNoisyGenes_v3.Rdata'))
# head(ngenes, 50)

plotexpmap(sc, 'Zfp42', logsc=TRUE, um=TRUE, noise=TRUE,cex=0.5)
#ngenes = ngenes[order(ngenes$pvalue), ]
#head(ngenes, 50)

# cells = colnames(sc@ndata)
# cl_new = cl$partition
# #newOrder = c(1, 2, 4, 5, 6, 3, 7, 8,9, 10, 11)
# cl_new[which(cl$partition == 4)] = 3
# cl_new[which(cl$partition == 5)] = 4
# cl_new[which(cl$partition == 6)] = 5
# cl_new[which(cl$partition == 3)] = 6

# cl_new = cl_new[order(cl_new)]
gamma = 4
load( file = paste0(outDir, '/out_varID2_noise_gamma_', gamma, '.Rdata'))
load(file = paste0(outDir, '/out_varID2_noise_diffNoisyGenes_bgs.Rdata'))

select_highNoisey_genes = function(ngenes, logfc_cutoff = 1, pval_cutoff = 10^-10)
{
  ngenes = ngenes[which(ngenes$log2FC > logfc_cutoff & ngenes$pvalue < pval_cutoff), ]
  mm = match(rownames(ngenes), c(tfs, sps, gene_examples))
  ngenes = ngenes[which(!is.na(mm)), ]
  return(ngenes)
}

ngenes_1 = select_highNoisey_genes(ngenes_1)
ngenes_2 = select_highNoisey_genes(ngenes_2)
ngenes_3 = select_highNoisey_genes(ngenes_3)
ngenes_4 = select_highNoisey_genes(ngenes_4)
ngenes_5 = select_highNoisey_genes(ngenes_5)

genes = unique(c(rownames(ngenes_1), rownames(ngenes_2), rownames(ngenes_3), 
                 rownames(ngenes_4), rownames(ngenes_5))
                 #'Pax6', 'Dhrs3')
                 )

saveRDS(genes, file = paste0(RdataDir, 'topGenes_localVaribility.gene.expression_VarID2.rds'))
#test = cl$partition
#newtest = test
#newtest[which(test == 4)] = 0
#newtest[which(test == 5)] = 4
#cl$partition = newtest
#xx <- updateSC(sc,res = res, cl=cl)
#plotmap(xx, um=TRUE, tp = 1, cex = 0.3)

#genes = head(genes, 400)

pdfname = paste0(outDir, '/plot_highNoise_genes_conditionsComparisions_bgs_137topgenes.pdf')
pdf(pdfname, width=8, height = 18)

ph <- plotmarkergenes(sc, genes=genes, noise=TRUE,
                      cluster_rows=TRUE, 
                      cl = c(1, 2, 3),
                      #cells = names(cl_new), 
                      order.cells = FALSE,
                      cluster_set = TRUE,
                      cluster_cols = FALSE,
                      cap = 5, 
                      zsc=TRUE, 
                      logscale = TRUE
                      )

dev.off()


pdfname = paste0(outDir, '/plot_noise_topGenes.pdf')
pdf(pdfname, width=8, height = 6)

for(n in 1:length(genes))
{
  plotexpmap(sc, genes[n], logsc=TRUE, um=TRUE, noise=TRUE, cex=0.5)
  
}

dev.off()
