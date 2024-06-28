##########################################################################
##########################################################################
# Project: RA competence project
# Script purpose: try the data integration across time points to estimate the pseudotime in palantir or scanpy
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Oct 11 10:58:41 2023
##########################################################################
##########################################################################
suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome)
  library(eisaR)
  library(GenomicFeatures)
  library(SummarizedExperiment)
  library(tximeta)
  library(rjson)
  library(reticulate)
  library(SingleCellExperiment)
  library(scater)
})

names(cols) = levels

levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA",  "day5_RA",
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA"
                )

cols_sel = cols[match(levels_sels, names(cols))]

#levels_sels = c("day2_beforeRA", 
#                "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA")
#cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/RA_symetryBreaking/dataIntegration_timePoints_4pseudotime/')
system(paste0('mkdir -p ', outDir))


##########################################
# import data and select the samples  
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObj_clustersFiltered_umapOverview.rds'))
Idents(aa) = factor(aa$condition, levels = levels)

aa = subset(aa, idents = levels_sels)

## filter the weird clusters identified in make_plots_v0.R
cells_2filter = readRDS(file = paste0(RdataDir, 'subObj_clusters_to_filter.rds'))
mm = match(colnames(aa), colnames(cells_2filter))
aa = subset(aa, cells = colnames(aa)[which(is.na(mm))])

# aa = readRDS(file = paste0(RdataDir, 
#                            'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
#                            'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d6_',
#                            species, version.analysis, '.rds'))

Idents(aa) = aa$condition
p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)

p1 + p2

## remove mature neurons 
aa = subset(aa, cells = colnames(aa)[which(is.na(aa$celltypes) | aa$celltypes != "Neurons")])


p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)

p1 + p2

ggsave(filename = paste0(outDir, 'UMAP_RAtreatment_d2.to.d5_noNeurons.pdf'), 
       width = 18, height = 8)


# rerun the umap 
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
Idents(aa) = aa$condition

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)

ggsave(filename = paste0(outDir, 'UMAP_RAtreatment_d2.to.d5.noMatureNeurons.pdf'), 
       width = 10, height = 6)

saveRDS(aa, file = paste0(RdataDir, 'seuratObj_clustersFiltered_RA.noRA_d2_d5_forHeterogeenity.rds'))


Discard_cellCycle.corrrelatedGenes = TRUE
if(Discard_cellCycle.corrrelatedGenes){
  library(scater)
  Idents(aa) = aa$condition
  
  # Identifying the likely cell cycle genes between phases,
  # using an arbitrary threshold of 5%.
  scaledMatrix = GetAssayData(aa, slot = c("scale.data"))
  
  diff <- getVarianceExplained(scaledMatrix, data.frame(phase = aa$Phase))
  diff = data.frame(diff, gene = rownames(diff))
  diff = diff[order(-diff$phase), ]
  
  hist(diff$phase, breaks = 100); 
  abline(v = c(1:5), col = 'red')
  
  rm(scaledMatrix)
  
  genes_discard = diff$gene[which(diff$phase>5)]
  cat(length(genes_discard), 'genes to discard \n')
  
  print(intersect(genes_discard, gene_examples))
  
  aa = subset(aa, features = setdiff(rownames(aa), genes_discard))
  
}

aa = FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)
Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)

saveRDS(aa, file = paste0(outDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d5.noNeurons.rds'))

aa = readRDS(file = paste0(outDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d5.noNeurons.rds'))


########################################################
########################################################
# Section : test data integration with seurat_rpca
# 
########################################################
########################################################
Test_CSS_integration = FALSE
if(Test_CSS_integration){
  aa = readRDS(file = paste0(outDir, 
                             'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                             'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d5.noNeurons.rds'))
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)
  
  aa$dataset = 'afterRA'
  aa$dataset[which(aa$condition == 'day2.5_RA')] = 'RA'
  aa$dataset[which(aa$condition == 'day2_beforeRA')] = 'beforeRA'
  
  aa = FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  
  ggsave(filename = paste0(outDir, 'UMAP_RAtreatment_d2.to.d5.noMatureNeurons_rmCellCycle.correlatedGenes.pdf'), 
         width = 10, height = 6)
  
  ##########################################
  # test CSS
  ##########################################
  library(simspec)
  aa = FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000)
  aa <- RunPCA(aa, npcs = 50, weight.by.var = FALSE, 
               features = VariableFeatures(object = aa), verbose = FALSE)
  
  aa <- cluster_sim_spectrum(object = aa, 
                             label_tag = "condition",
                             use_dr = 'pca', 
                             #dims_use = 1:30, 
                             cluster_resolution = 0.5
  )
  
  aa <- RunUMAP(aa, reduction = "css", dims = 1:ncol(Embeddings(aa, "css")),
                n.neighbors = 30, min.dist = 0.1)
  DimPlot(aa, reduction = "umap", label = TRUE, repel = TRUE,  cols = cols_sel)
  
  
  FeaturePlot(aa, features = c('Zfp42', 'Foxa2', 'Pax6', 'Rarg', 'Dhrs3'))
  
  ggsave(filename = paste0(outDir, 
                           'CSS_integration_beforeRA_RA_afterRA.pdf'), 
         width = 10, height = 6)
  
  
  #aa <- FindNeighbors(aa, reduction = "css", dims = 1:ncol(Embeddings(aa, "css")))
  #aa <- FindClusters(aa, resolution = 1)
  #UMAPPlot(aa, group.by = "batch") + UMAPPlot(aa, group.by = "celltype_RSS") + UMAPPlot(aa)
  
  aa <- cluster_sim_spectrum(object = aa, 
                             label_tag = "condition",
                             #cluster_col = 'condition',
                             use_dr = 'pca', 
                             #dims_use = 1:10, 
                             spectrum_type = "corr_kernel",
                             corr_method = "spearman",
                             cluster_resolution = 0.5
                             
  )
  
  nb_css = ncol(Embeddings(aa, "css"))
  aa <- RunUMAP(aa, reduction = "css", dims = 1:nb_css, n.neighbors = 30, min.dist = 0.1)
  DimPlot(aa, reduction = "umap", label = TRUE, repel = TRUE, cols = cols_sel)
  
  ggsave(filename = paste0(outDir, 
                           'CSS_integration_kernelProb_condition_test_v2.pdf'), 
         width = 10, height = 6)
  
  FeaturePlot(aa, features = c('Zfp42', 'Foxa2', 'Pax6', 'Rarg', 'Dhrs3'))
  
  
}

##########################################
# save file for palantir
##########################################
library(SeuratDisk)

mnt = aa

VariableFeatures(mnt) = NULL
#mnt@assays$RNA@scale.data = NULL
#mnt@assays$RNA@data = NULL

DefaultAssay(mnt) = 'RNA'
mnt = DietSeurat(mnt, 
                 counts = TRUE, 
                 data = TRUE,
                 scale.data = FALSE,
                 features = rownames(mnt), 
                 assays = c('RNA'), 
                 dimreducs = c('pca', 'umap'), graphs = NULL, 
                 misc = TRUE
)


DefaultAssay(mnt) = 'RNA'
VariableFeatures(mnt)


#saveRDS(mnt, file = paste0(outDir, '/SeuratObj_splice.unspliced_RA_day2_to_day6_all.rds'))
Idents(mnt) = mnt$condition

#mnt = subset(mnt, downsample = 1000)
#saveDir = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/",
#                "results/scRNAseq_R13547_10x_mNT_20220813/RA_symetryBreaking/")

saveFile = '/RNAmatrix_RA_noRA_d2_d5.h5Seurat'

SaveH5Seurat(mnt, filename = paste0(outDir, saveFile), 
             overwrite = TRUE)
Convert(paste0(outDir, saveFile), 
        dest = "h5ad", overwrite = TRUE)


## manually select starting cell and terminal cells
DimPlot(aa, reduction = "umap", label = TRUE, repel = TRUE, cols = cols_sel)

FeaturePlot(aa, features = c('Nanog', 'Pou5f1', 'Foxa2', 'Pax6'))
xx = Embeddings(aa, reduction = 'umap')
kk = which(xx[,1] < (0) & xx[,2] < -6)


DimPlot(aa, cells.highlight = c('ACTCCCACAGCTCTGG-1_1_1_1_1_1_1_1_1', 'TTGTGTTAGATTGACA-1_2_1_1',
                                'GCTCAAACATCGGAAG-1_1_1_1'
                                ), 
        cols.highlight = "#DE2D26", sizes.highlight = 2)


##########################################
# test Seurat integration
##########################################
Test_Seurat_Integration = FALSE
if(Test_Seurat_Integration){
  ref.list <- SplitObject(aa, split.by = "dataset")
  #rm(list = c('refs.merged', 'aa', 'srat')) # remove big seurat objects to clear memory
  
  # normalize and identify variable features for each dataset independently
  ref.list <- lapply(X = ref.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize")
    x <- FindVariableFeatures(x, selection.method = "vst")
  })
  
  # select features that are repeatedly variable across datasets for integration run PCA on each
  # dataset using these features
  #features <- SelectIntegrationFeatures(object.list = ref.list, nfeatures = 5000)
  #hvgs1 = FindVariableFeatures(ref.list[[1]], selection.method = "vst", nfeatures = 500)
  #hvgs2 = FindVariableFeatures(ref.list[[2]], selection.method = "vst", nfeatures = 500)
  #hvgs3 = FindVariableFeatures(ref.list[[3]], selection.method = "vst", nfeatures = 2000)
  #features = unique(c(VariableFeatures(aa), hvgs1, hvgs2, hvgs3))
  features = VariableFeatures(aa)
  
  ref.list <- lapply(X = ref.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = TRUE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                        anchor.features = features, 
                                        #reference = c(2),
                                        reduction = "cca", 
                                        #reduction = 'rpca',
                                        k.anchor = 5,
                                        dims = 1:30)
  
  #rm(ref.list)
  
  # this command creates an 'integrated' data assay
  ref.combined <- IntegrateData(anchorset = ref.anchors,
                                #sample.tree = matrix(c(-5, 1, 2, 3, 4, -6, -4, -3, -2,-1), ncol = 2),
                                sample.tree = matrix(c(-2, 1, -3, -1), ncol = 2),
                                dims = 1:30,
                                verbose = TRUE
  ) ## take ~100G memory
  
  #rm(ref.anchors)
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(ref.combined) <- "integrated"
  
  ref.combined <- ScaleData(ref.combined, verbose = FALSE)
  ref.combined <- RunPCA(ref.combined, npcs = 50, verbose = FALSE)
  
  ElbowPlot(ref.combined, ndims = 50)
  
  #ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
  #ref.combined <- FindClusters(ref.combined, resolution = 0.2)
  ref.combined <- RunUMAP(ref.combined, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
  DimPlot(ref.combined, reduction = "umap", cols = cols_sel)
  
  ggsave(filename = paste0(outDir, 'UMAP_integrationTest_CCA.pdf'), 
         width = 10, height = 6)
  
  #saveRDS(ref.combined, file = paste0(outDir, '/integrated_mNT_mouseGastrulation.rds'))
  rm(ref.combined)
  
}

##########################################
# test package psupertime 
##########################################
Test_psupertime = FALSE
if(Test_psupertime){
  
  # load psupertime package
  suppressPackageStartupMessages({
    library('psupertime')
    library('SingleCellExperiment')
    library(scuttle)
    library(scran)
    library(tictoc)
  })
  
  Idents(aa) = aa$condition
  sce = as.SingleCellExperiment(aa)
  
  dec <- modelGeneVar(sce)
  plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
  curve(metadata(dec)$trend(x), col="blue", add=TRUE)
  top.hvgs <- getTopHVGs(dec, prop=0.05, var.threshold = 0.0)
  
  # top.hvgs = getTopHVGs(dec, fdr.threshold=0.1)
  cat(length(top.hvgs), ' HVGs \n')
  
  # load the data
  #data(acinar_hvg_sce)
  
  # run psupertime
  tic()
  y           = as.numeric(sce$time)
  #psuper_obj  = psupertime(sce, y, sel_genes='all')
  #psuper_hvg      = psupertime(sce, y, sel_genes='hvg')
  psuper_hvg_custom  = psupertime(sce, y, 
                                  sel_genes='list', 
                                  gene_list=top.hvgs)
  #psuper_hvg_custom2  = psupertime(acinar_hvg_sce, y, sel_genes=list(hvg_cutoff=0.1, bio_cutoff=0.5, span=0.1))
  #psuper_hvg
  
  toc()
  
  saveRDS(psuper_hvg_custom, file = paste0(outDir, 'pseudotime_psupertime_saved_all.cells_custome.hvgs.rds'))
  rm(sce)
  
  xx = readRDS(file = paste0(outDir, 'pseudotime_psupertime_saved_all.rds'))
  
  psuper_obj = readRDS(file = paste0(outDir, 'pseudotime_psupertime_saved_all.cells_custome.hvgs.rds'))
  
  g       = plot_train_results(psuper_obj)
  (g)
  
  g       = plot_labels_over_psupertime(psuper_obj, label_name='time')
  (g)
  
  g       = plot_identified_genes_over_psupertime(psuper_obj, label_name='time')
  (g)
  
}


########################################################
########################################################
# Section IV: process the Palantir pseudotime and double check the imputation by MAGIC
# 
########################################################
########################################################
require(data.table)

if(USE_CSS){
  pt = read.csv(file = paste0(outDir, 'palantir_pseudotime_all.csv'), header = TRUE, row.names = 1)
  impute = fread(paste0(outDir, 'MAGIC_imputed_data.csv'))
  impute = impute[, -1]
  impute = impute[-1, ]
  impute = t(impute)
  
  cells = read.csv(paste0(outDir, 'palantir_cells.csv'), header = TRUE, row.names = 1)
  genes = read.csv(paste0(outDir, 'palantir_genes.csv'), header = TRUE, row.names = 1)
  
  #fates = read.csv(paste0(outDir, 'palantir_pseudotime_cellfates.csv'), header = TRUE, row.names = 1)
  #pt = data.frame(pt, fates[match(rownames(pt), rownames(fates)), ])
  #pt = pt[, c(1, 3, 5:6)]
  #rm(fates)
  colnames(impute) = rownames(cells)
  rownames(impute) = rownames(genes)
  
  rm(genes)
  
  mm = match(rownames(impute), c(gene_examples, tfs, sps))
  mm = which(!is.na(mm))
  impute = impute[mm, ]
  
  pt = pt[order(pt$palantir_pseudotime), ]
  impute = impute[ ,match(rownames(pt), colnames(impute))]
  
  save(impute, pt, file = paste0(outDir, 'MAGIC_imputatedData_palantir_pseudotime_fateProb_CSS.Rdata'))
  
}else{
  pt = read.csv(file = paste0(outDir, 'palantir_pseudotime_all_noCSS.csv'), header = TRUE, row.names = 1)
  
  aa$pdt = NA
  aa$pentropy = NA
  mm = match(colnames(aa), rownames(pt))
  aa$pdt = pt$palantir_pseudotime[mm]
  aa$pentropy = pt$palantir_entropy[mm]
  
  FeaturePlot(aa, features = c('pdt')) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  psuper_obj = readRDS(file = paste0(outDir, 
                                     'pseudotime_psupertime_saved_all.cells_custome.100.hvgs.rds'))
  psuper_obj = psuper_obj$proj_dt
  pt$psuper = psuper_obj$psuper[match(rownames(pt), psuper_obj$cell_id)]
  rm(psuper_obj)
  
  aa$psuper = pt$psuper[match(colnames(aa), rownames(pt))]
  
  FeaturePlot(aa, features = c('psuper')) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  
  impute = fread(paste0(outDir, 'MAGIC_imputed_data_noCSS.csv'))
  impute = impute[, -1]
  impute = impute[-1, ]
  impute = t(impute)
  
  cells = read.csv(paste0(outDir, 'palantir_cells_noCSS.csv'), header = TRUE, row.names = 1)
  genes = read.csv(paste0(outDir, 'palantir_genes_noCSS.csv'), header = TRUE, row.names = 1)
  fates = read.csv(paste0(outDir, 'palantir_pseudotime_cellfates.csv'), header = TRUE, row.names = 1)
  
  pt = data.frame(pt, fates[match(rownames(pt), rownames(fates)), ])
  #pt = pt[, c(1, 3, 5:6)]
  rm(fates)
  
  colnames(impute) = rownames(cells)
  rownames(impute) = rownames(genes)
  
  rm(genes)
  
  mm = match(rownames(impute), c(gene_examples, tfs, sps))
  mm = which(!is.na(mm))
  
  impute = impute[mm, ]
  
  save(impute, pt, file = paste0(outDir, 'MAGIC_imputatedData_palantir_pseudotime_fateProb.Rdata'))
 
  load(file = paste0(outDir, 'MAGIC_imputatedData_palantir_pseudotime_fateProb.Rdata'))
  
  kk = which(rownames(impute) == 'Foxa2')
  plot(pt$palantir_pseudotime, impute[kk, ], cex = 0.1)
  plot(pt$psuper, impute[kk, ], cex = 0.1)
  
}


#plot(cells$palantir_pseudotime, impute[which(rownames(impute) == 'Rarg'), ], cex = 0.1)
##########################################
# plot the distribution of imputated gene expression
##########################################
library(ggplot2)
library(ggridges)
theme_set(theme_minimal())

noisyGenes = readRDS(file = paste0(RdataDir, 'topGenes_localVaribility.gene.expression_VarID2.rds'))

if(!is.null(dev.list())) dev.off()
pdf(paste0(outDir, '/Global_distribution_MAGICimputation_noiseGenes_noCSS.pdf'),
    width =8, height = 6, useDingbats = FALSE)

for(g in noisyGenes)
{
  # g = 'Foxa2'
  cat(g, '--\n')
  xx = data.frame(expression = impute[which(rownames(impute) == g), ], 
                  cells = colnames(impute), 
                  #pseudt = cells$palantir_pseudotime,
                  condition = aa$condition[match(colnames(impute), colnames(aa))])
  p1 = ggplot(xx, aes(x = expression, y = condition)) +
    geom_density_ridges(aes(fill = condition)) +
    scale_fill_manual(values = cols_sel) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, size = 12), 
          axis.text.y = element_text(angle = 0, size = 12), 
          axis.title =  element_text(size = 14),
          legend.text = element_text(size=12),
          legend.title = element_text(size = 12)
    ) + 
    labs(x = "Imputed expression by MAGIC", y = "") +
    ggtitle(g)
  
  plot(p1)
  
}

dev.off()

##########################################
# plot the gene trend for two branches 
##########################################
#pt = pt[order(pt$palantir_pseudotime), ]
#impute = impute[ ,match(rownames(pt), colnames(impute))]

kk = which(rownames(impute) == 'Foxa2')
plot(pt$palantir_pseudotime, impute[kk, ], cex = 0.1)
plot(pt$psuper, impute[kk, ], cex = 0.1)


pt$branch = 'none'
pt$branch[pt$NP == 'False' & pt$FP == 'True'] = 'fp'
pt$branch[pt$NP == 'True' & pt$FP == 'False'] = 'np'

get_smooth_curve_spline = function(x, t, newt, downsample = TRUE)
{
  # x = as.numeric(cds[1, jj]); t = Pseudotime; newt = pseudot_comomon;
  if(downsample){
    nb_t = min(5000, length(t))
    nn = sample(1:length(t), size = nb_t, replace = FALSE)
    t = t[nn]
    x = x[nn]
  }
  
  fit_sel = smooth.spline(t, x, df = 3)
  
  #plot(Pseudotime, cds_sel, cex = 0.5)
  #lines(fit_sel, col = 'red', lwd =2.0)
  newx = predict(fit_sel, newt)
  return(newx$y)
  #VGAM::vglm(~sm.ns(Pseudotime, df=3), family = 'gaussian', data = cds_sel)
  
}

for(branch in c('fp', 'np'))
{
  # branch = 'np'
  jj = which(pt$branch == 'none'|pt$branch == branch)
  tt = pt$palantir_pseudotime[jj];
  bx = seq(from = min(tt), 
           to = max(tt),
           length.out = 100)
  
  trend = matrix(NA, nrow = nrow(impute), ncol = length(bx))
  rownames(trend) = rownames(impute)
  
  for(n in 1:nrow(impute))
  {
    cat(n, '--', rownames(impute)[n], '\n')
    #kk = which(rownames(impute) == g)
    #plot(pt$palantir_pseudotime, impute[kk, ], cex = 0.1)
    #plot(pt$palantir_pseudotime[jj], impute[kk, jj], cex = 0.1)
    #bx <- c(0,50,100,150,200)+0.5
    #yS <- binMeans(y =impute[kk, jj] , x=tt, bx=bx)
    yS = get_smooth_curve_spline(x = impute[n, jj], t = tt, newt = bx, downsample = TRUE)
    #plot(bx, yS, type = 'l')
    trend[n, ] = yS
  }
  
  save(bx, trend, file = paste0(outDir, 'gene_trends_branch_', branch, '.Rdata'))
  
}

load(file = paste0(outDir, 'gene_trends_branch_fp.Rdata'))



