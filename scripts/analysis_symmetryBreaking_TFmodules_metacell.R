##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: search for genes with asymmetric expression at early time points
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Mar 30 16:04:20 2023
##########################################################################
##########################################################################
library(Seurat)
library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(future)

options(future.globals.maxSize = 120000 * 1024^2)

source('functions_utility.R')

outDir = paste0(resDir, '/RA_symetryBreaking/TF_modules/')
system(paste0('mkdir -p ', outDir))

levels_sels = c("day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA")

names(cols) = levels
cols_sel = cols[match(levels_sels, names(cols))]


##########################################
# import the all data for RA treatment
##########################################
#aa = readRDS(file = paste0(RdataDir,
#                           'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
#                           'cellCycleScoring_annot.v1_savedUMAP.subs.v2_', species, version.analysis, '.rds'))
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d6_',
                           species, version.analysis, '.rds'))


DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE, pt.size = 1)

ggsave(filename = paste0(outDir, '/UMAP_conditions.pdf'), 
       width = 10, height = 6)


Idents(aa) = aa$condition
p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)

p1 + p2

ggsave(filename = paste0(outDir, 'UMAP_RAtreatment_d2.to.d6_.pdf'), 
       width = 18, height = 8)

FeaturePlot(aa, features = c('Zfp42', 'Tcf15', 'Skil', 'Lef1'))
ggsave(filename = paste0(outDir, '/asymmetric_feature_expression.pdf'), 
       width = 14, height = 10)

FeaturePlot(aa, features = c('Pax6', 'Foxa2', 'Sox1'))
ggsave(filename = paste0(outDir, '/asymmetric_feature_expression_v2.pdf'), 
       width = 14, height = 10)

p1 = FeaturePlot(aa, features = c("Pax6", "Sox1"), blend = TRUE)
p2 = FeaturePlot(aa, features = c("Foxa2", "Sox1"), blend = TRUE)

p1 /p2

ggsave(filename = paste0(outDir, '/Overlapping_Foxa2_Pax6_Sox1.pdf'), 
       width = 14, height = 8)


FeatureScatter(aa, feature1 = "Foxa2", feature2 = "Pax6", group.by = 'condition')
FeatureScatter(aa, feature1 = "Foxa2", feature2 = "Sox1", group.by = 'condition')
FeatureScatter(aa, feature1 = "Pax6", feature2 = "Sox1", group.by = 'condition')


##########################################
# clean and subset the samples for scFates test 
##########################################
data_version = 'd2.5_d5_TFs_SPs_metacell_v1'
#data_version = 'd2.5_d5_TFs_SPs_metacell_v1.5'

system(paste0('mkdir -p ', outDir, data_version))


Clean_Subset_for_scFates = FALSE
if(Clean_Subset_for_scFates){
  ## remove the cluster 8 and 9 mainly mature neurons and also day6_RA
  # aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '8' & aa$celltypes != '9' & 
  #                                              aa$condition != 'day6_RA')])
  # 
  # aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '8' & aa$celltypes != '9')])
  # 
  # aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '8' & aa$celltypes != '9' & 
  #                                              aa$condition != 'day2_beforeRA')])
  # 
  # aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '8' & aa$celltypes != '9' & 
  #                                              aa$condition != 'day2_beforeRA' & 
  #                                              aa$condition != 'day2.5_RA')])
  
  aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '8' & aa$celltypes != '9' & 
                                               aa$condition != 'day2_beforeRA' & 
                                               aa$condition != 'day6_RA')])
  
  # aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '8' & aa$celltypes != '9' & 
  #                                              aa$condition != 'day2_beforeRA' &
  #                                              aa$condition != 'day2.5_RA' &
  #                                              aa$condition != 'day6_RA')])
  
  # aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != '8' & aa$celltypes != '9' & 
  #                                              aa$condition != 'day2_beforeRA' & 
  #                                              aa$condition != 'day6_RA' & 
  #                                              aa$condition != 'day5_RA')])
  # 
  # first subset the cells 
  table(aa$condition)
  Idents(aa) = aa$condition
  
  ## test metacell approach 
  ## original code from 
  ## https://gfellerlab.github.io/MetacellAnalysisTutorial/Metacell-
  ## construction-chapter.html#SuperCell-construction
  library(SuperCell)
  sc_data = aa
  annotation_label = 'condition'
  
  #sc_data <- NormalizeData(sc_data, normalization.method = "LogNormalize")
  sc_data <- FindVariableFeatures(sc_data, nfeatures = 3000)
  #sc_data <- ScaleData(sc_data)
  #> Centering and scaling data matrix
  sc_data <- RunPCA(sc_data, npcs = 100, verbose = F)
  
  sc_data <- RunUMAP(sc_data, reduction = "pca", dims = c(1:100), n.neighbors = 50, verbose = F, 
                     min.dist = 0.2)
  UMAPPlot(sc_data, group.by = "condition")
  
  ggsave(filename = paste0(outDir, data_version, '/plot_UMAP_for_metacell_v2.pdf'), 
         width = 10, height = 6)
  
  FeaturePlot(sc_data, features = c('Foxa2', 'Pax6'))
  
  
  # source(paste0(functionDir, '/functions_scRNAseq.R'))
  # explore.umap.params.combination(sub.obj = sc_data, 
  #                                 resDir = outDir, 
  #                                 pdfname = 'UMAP_params_for_metacells.pdf',
  #                                 use.parallelization = FALSE,
  #                                 group.by = 'condition',
  #                                 cols = cols_sel, 
  #                                 nfeatures.sampling = c(3000, 5000),
  #                                 nb.pcs.sampling = c(30, 50, 100), 
  #                                 n.neighbors.sampling = c(30, 50, 100),
  #                                 min.dist.sampling = c(0.2, 0.3, 0.5)
  #                                 )
                                  
  
  gamma = 50 # the requested graining level.
  k_knn = 30 # the number of neighbors considered to build the knn network.
  nb_var_genes = 2000 # number of the top variable genes to use for dimensionality reduction 
  nb_pc = 30 # the number of principal components to use.   
  
  MC <- SuperCell::SCimplify(Seurat::GetAssayData(sc_data, slot = "data"),  
                             k.knn = k_knn,
                             gamma = gamma,
                             # n.var.genes = nb_var_genes,  
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
  
  pdfname = paste0(outDir, data_version, '/plot_metacell_graph_v2.pdf')
  pdf(pdfname, width=10, height = 6)
  
  supercell_plot(
    MC$graph.supercells, 
    group = MC$annotation,
    lay.method = 'fr',
    seed = 1, 
    alpha = -pi/2,
    main  = "Metacells colored by condition"
  )
  
  # supercell_plot_UMAP(
  #   MC$graph.supercells, 
  #   group = MC$annotation,
  #   #lay.method = 'fr',
  #   #seed = 1, 
  #   alpha = -pi/2,
  #   title  = "Metacells colored by condition"
  # )
  # 
  dev.off()
  
  
  ## save the metacell result
  colnames(MC.GE) <- as.character(1:ncol(MC.GE))
  
  MC.seurat <- CreateSeuratObject(counts = MC.GE, 
                                  meta.data = data.frame(size = as.vector(table(MC$membership)))
  )
  MC.seurat[[annotation_label]] <- MC$annotation
  
  # save single-cell membership to metacells in the MC.seurat object
  MC.seurat@misc$cell_membership <- data.frame(row.names = names(MC$membership), membership = MC$membership)
  MC.seurat@misc$var_features <- MC$genes.use 
  
  # Save the PCA components and genes used in SCimplify  
  PCA.res <- irlba::irlba(scale(Matrix::t(sc_data@assays$RNA@data[MC$genes.use, ])), nv = nb_pc)
  pca.x <- PCA.res$u %*% diag(PCA.res$d)
  rownames(pca.x) <- colnames(sc_data@assays$RNA@data)
  
  MC.seurat@misc$sc.pca <- CreateDimReducObject(
    embeddings = pca.x,
    loadings = PCA.res$v,
    key = "PC_",
    assay = "RNA"
  )
  
  if(packageVersion("Seurat") >= 5) {
    MC.seurat[["RNA"]] <- as(object = MC.seurat[["RNA"]], Class = "Assay")
  }
  
  saveRDS(MC.seurat, file = paste0(outDir, data_version,
                                   '/seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                                   'cellCycleScoring_annot.v2_newUMAP_clusters_time_metacell_SuperCell',
                                   '_v2.rds'))
  
  
  ## QC of metacell (https://gfellerlab.github.io/MetacellAnalysisTutorial/QCs.html)
  QC_metacells = FALSE 
  if(QC_metacells){
    
    library(reticulate)
    conda_env <-  conda_list()[reticulate::conda_list()$name == "palantir","python"]
    
    use_condaenv(conda_env)
    
    library(MetacellAnalysisToolkit)
    
    mc_data = MC.seurat
    membership_df <- mc_data@misc$cell_membership
    
    mc_data$purity <- mc_purity(membership = membership_df$membership, 
                                annotation = sc_data@meta.data[, annotation_label])
    qc_boxplot(mc.obj = mc_data, qc.metrics = "purity")
    
    ggsave(filename = paste0(outDir, data_version, '/QCs_metacell_purity.pdf'), 
           width = 10, height = 6)
    
    qc_boxplot(mc.obj = mc_data, qc.metrics = "purity", split.by = annotation_label)
    ggsave(filename = paste0(outDir, data_version, '/QCs_metacell_purity_condition.pdf'), 
           width = 10, height = 6)
    
    pdfname = paste0(outDir, data_version, '/plot_metacell_sizeDistribution.pdf')
    pdf(pdfname, width=10, height = 6)
    
    # Seurat::VlnPlot(mc_data, features = "size", pt.size = 2)
    # Seurat::VlnPlot(mc_data, features = "size", pt.size = 2, group.by = annotation_label)
    hist(mc_data$size, main = "Size distribution", xlab = "Size", breaks = 50)
    
    dev.off()
    
    
    # we run diffufion maps on the pca axes saved in the seurat object. 
    # Note that a matrix containing the PCA components can be provided as well
    diffusion_comp <- get_diffusion_comp(sc.obj = sc_data, sc.reduction = "pca", dims = 1:30)
    
    #> Error in python_config_impl(python) : 
    #>   Error running '/users/agabrie4/.virtualenvs/r-reticulate/bin/python': No such file.
    #> The Python installation used to create the virtualenv has been moved or removed:
    #>   '/usr/bin'
    #> Warning: The following arguments are not used: layer
    #> Computing diffusion maps ...
    mc_data$compactness <- mc_compactness(cell.membership = membership_df, 
                                          sc.obj = sc_data,
                                          sc.reduction = diffusion_comp, 
                                          dims = 1:ncol(diffusion_comp))
    qc_boxplot(mc.obj = mc_data, qc.metrics = "compactness")
    
    
    sc_data <- RunUMAP(sc_data, reduction = "pca", dims = c(1:30), n.neighbors = 30, 
                       #min.dist = 0.2, 
                       verbose = F)
    UMAPPlot(sc_data, group.by = "condition", cols = cols_sel)
    
    pdfname = paste0(outDir, data_version, '/plot_metacell_projection_v2.pdf')
    pdf(pdfname, width=10, height = 6)
    
    mc_projection(
      sc.obj = sc_data, 
      mc.obj = mc_data,
      cell.membership = membership_df,
      sc.reduction = "umap", 
      sc.label = unlist(annotation_label), # single cells will be colored according the sc.label  
      metacell.label = unlist(annotation_label) # metacells cell will be colored according the metacell.label
    )
    
    dev.off()
    
  }
  
  
  aa = readRDS(file = paste0(outDir, data_version,
                             '/seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                             'cellCycleScoring_annot.v2_newUMAP_clusters_time_metacell_SuperCell',
                             '_v2.rds'))
  
  #aa = subset(aa, downsample = 2000)
  
  table(aa$condition)
  
  aa <- NormalizeData(aa)
  aa <- FindVariableFeatures(aa, selection.method = "vst", 
                             nfeatures = 3000)
  
  aa <- ScaleData(aa, vars.to.regress = 'nCount_RNA')
  
  ## double check the Foxa2 and Pax6 levels 
  xx = data.frame(aa@assays$RNA@data)
  jj = which(rownames(xx) == 'Foxa2'|rownames(xx) == 'Pax6')
  plot(as.numeric(xx[jj[1],]), as.numeric(xx[jj[2], ]))
  
  
  ## calculate the cell cycle scoreing
  s.genes <- firstup(cc.genes$s.genes)
  s.genes[which(s.genes == 'Mlf1ip')] = 'Cenpu'
  g2m.genes <- firstup(cc.genes$g2m.genes)
  
  aa <- CellCycleScoring(aa, s.features = s.genes, 
                         g2m.features = g2m.genes, 
                         set.ident = FALSE)
  
  # Visualize the distribution of cell cycle markers across
  RidgePlot(aa, features = c("Pcna", "Top2a", "Mcm6", "Mki67"),
            ncol = 2)
  
  # Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
  # phase
  aa <- RunPCA(aa, features = c(s.genes, g2m.genes))
  DimPlot(aa, reduction = 'pca', group.by = 'Phase')
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), 
               verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 50, min.dist = 0.3)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', 
          cols = cols_sel, raster=FALSE)
  
  ggsave(filename = paste0(outDir, data_version, '/UMAP_RAtreatment', 
                           '_metacell.pdf'), 
         width = 10, height = 6)
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', 
          raster=FALSE)
  
  Discard_cellCycle.corrrelatedGenes = FALSE
  if(Discard_cellCycle.corrrelatedGenes){
    library(scater)
    Idents(aa) = aa$condition
    
    DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
    
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
    
    tfs_sels = intersect(genes_discard, gene_examples)
    print(tfs_sels)
    
    if(length(tfs_sels)>0) genes_discard = setdiff(genes_discard, tfs_sels)
    
    tfs_sels = intersect(genes_discard, gene_examples)
    print(tfs_sels)
    
    aa = subset(aa, features = setdiff(rownames(aa), genes_discard))
    
  }else{
    
    DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase',  
            raster=FALSE)
    
    # subset first maybe
    Idents(aa) = aa$condition
    #aa = subset(aa, downsample = 2000)
    
    aa <- ScaleData(aa, vars.to.regress = c("nCount_RNA", 
                                            "S.Score", 
                                            "G2M.Score"), 
                    features = rownames(aa))
    
    saveRDS(aa, file = paste0(outDir, data_version,
                          '/seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                              'cellCycleScoring_annot.v2_newUMAP_clusters_time_cellCycole.regression_',
                              'metacells_gamma50', '.rds'))
    
    aa = readRDS(file = paste0(outDir, data_version,
                               '/seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                               'cellCycleScoring_annot.v2_newUMAP_clusters_time_cellCycole.regression_',
                               'metacells_gamma50', '.rds'))
    
    
  }
  
  
  aa = readRDS(file = paste0(outDir, data_version,
                             '/seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                             'cellCycleScoring_annot.v2_newUMAP_clusters_time_cellCycole.regression_',
                             data_version, '.rds'))
  
  aa = FindVariableFeatures(aa, selection.method = "vst", 
                            nfeatures = 2000)
  
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), 
               verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.2)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', 
          cols = cols_sel, raster=FALSE)
  
  
  saveRDS(aa, file = paste0(outDir, data_version,
                            '/seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                            'cellCycleScoring_annot.v2_newUMAP_clusters_time_cellCycole.regression_',
                            'metacells_gamma50', '.rds'))
  
  
  ## save the object for scFates
  aa$dataset = 'afterRA'
  aa$dataset[which(aa$condition == 'day2.5_RA')] = 'RA'
  aa$dataset[which(aa$condition == 'day2_beforeRA')] = 'beforeRA'
  
  mm = match(rownames(aa), c(tfs, sps, gene_examples))
  mnt = subset(aa, features = rownames(aa)[which(!is.na(mm))])
  
  library(SeuratDisk)
  
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
  
  saveFile = paste0('/RNAmatrix_RA_TFs_SPs_regressedCellCycle', 
                    data_version, '.h5Seurat')
  
  SaveH5Seurat(mnt, filename = paste0(outDir, data_version,  saveFile), 
               overwrite = TRUE)
  Convert(paste0(outDir, data_version,  saveFile), 
          dest = "h5ad", overwrite = TRUE)
  
  saveRDS(aa, file = paste0(outDir, data_version,
                            '/seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                            'cellCycleScoring_annot.v2_newUMAP_clusters_time_metacell_',
                            'cellcycleRegressed.rds'))
  
  
}


########################################################
########################################################
# Section I: plot the gene-gene correlation of scFates  
# 
########################################################
########################################################
library(data.table)

dataDir = paste0(outDir, 'd2.5_d5_TFs_SPs_metacell_v1/')
aa = readRDS(file = paste0(outDir, data_version,
                           '/seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_time_cellCycole.regression_',
                           'metacells_gamma50', '.rds'))

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', 
        cols = cols_sel, raster=FALSE)

ggsave(filename = paste0(outDir, data_version, '/UMAP_metacell.pdf'), width = 10, height = 6)

#fit = read.csv(file = paste0(dataDir, 'annData_layer_fitted.csv'), header = TRUE, row.names = c(1))
fit = read.csv(file = paste0(dataDir, 'annData_geneExpr_metacells.csv'), 
                    header = TRUE, row.names = c(1))

USE_MAGIC_allGenes = FALSE
if(USE_MAGIC_allGenes){
  fit.all = fread(file = paste0(dataDir, 'annData_magic_impuated_allGenes.csv'), header = TRUE)
  fit.all = data.frame(fit.all)
  rownames(fit.all) = fit.all$V1
  fit.all = fit.all[, -1]
  fit.sel = fit.all[, match(colnames(fit), colnames(fit.all))]
  #fit.sel = data.frame(fit.sel)
  fit = fit.sel
  
  rm(fit.all)
  rm(fit.sel)
  #rownames(impuated) = rownames(fit)
  #colnames(impuated) = colnames(fit)
  #fit = impuated
}


pst = read.csv(file = paste0(dataDir, 'annData_pseudotime_segments_milestones.csv'), header = TRUE,
               row.names = c(1))

assign = read.csv(file = paste0(dataDir, 'annData_cellAssignment_to_nonIntersectingWindows_8.csv'),
                  header = TRUE, row.names = c(1))

modules = read.csv(file = paste0(dataDir, 'scFates_early_late_modules_NP_FP_tl.test.fork.rescaled.csv'),
                   header = TRUE, row.names = c(1))

rownames(pst) = paste0('X', rownames(pst))
rownames(fit) = paste0('X', rownames(fit))

#fit = t(aa@assays$RNA@data)
mm = match(paste0('X', colnames(aa)), colnames(assign))
assign = assign[,mm]
assignment = rownames(assign)[apply(assign, 2, which.max)] 
aa$window = paste0('w', assignment)

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=8)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'window', pt.size = 1.5, label.size = 6.0,
        cols = color_list, raster=FALSE)

ggsave(filename = paste0(outDir, data_version, '/UMAP_metacell_scFates_slidingWindow.pdf'), 
       width = 10, height = 6)



ggs = rownames(modules)
ggs[which(ggs == "Nkx6-2")] = "Nkx6.2"
jj = match(ggs, colnames(fit))
fit = fit[, jj]

source('functions_utility.R')
p = make_scatterplot_genePair_scFates(geneA = 'Pax6', geneB = 'Foxa2', fit, assign, pst)
ggsave(filename = paste0(outDir, 'gene_gene_correlation_scFates_Foxa2_Pax6_v2.pdf'),
       width = 10, height = 6)

pdfname = paste0(outDir, data_version, '/genePairs_scatterplot_scFates_vs_Foxa2.pdf')
pdf(pdfname, width=12, height = 8)

source('functions_utility.R')
make_scatterplot_genePair_scFates_all(geneToCompare = c('Foxa2'), fit, assign, pst)

dev.off()

pdfname = paste0(outDir, data_version, '/genePairs_scatterplot_scFates_vs_Pax6.pdf')
pdf(pdfname, width=12, height = 6)
source('functions_utility.R')
make_scatterplot_genePair_scFates_all(geneToCompare = c('Pax6'), fit, assign, pst)

dev.off()


########################################################
########################################################
# Section II : narrow down the TFs and SPs
# 
########################################################
########################################################
require(scater)
require(scran)

ggs = rownames(aa)
ggs = ggs[which(!is.na(match(ggs, c(tfs, sps))))]

sce <- as.SingleCellExperiment(subset(aa, features = ggs))

#rownames(sce) = toupper(rownames(sce))
ave.counts <- calculateAverage(sce, assay.type = "counts")

hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))

num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
              xlab=expression(Log[10]~"average count"))

# detected in >= 5 cells, ave.counts >=5 but not too high
genes.to.keep <- num.cells > 100 & ave.counts >= 10^-2
summary(genes.to.keep)

sce <- sce[genes.to.keep, ]

ggs = rownames(sce)

length(ggs[!is.na(match(ggs, tfs))])



########################################################
########################################################
# Section : subclustering each time points
# 
########################################################
########################################################
p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
#p2 = DimPlot(aa, group.by = 'clusters', label = TRUE, repel = TRUE)

conditions = unique(aa$condition)
print(conditions)
bb = aa;
rm(aa)

conditions = c('day2_beforeRA', 'day2.5_RA', 'day3_RA.rep1')
noisyGenes = readRDS(file = paste0(RdataDir, 'topGenes_localVaribility.gene.expression_VarID2.rds'))
gene_examples = unique(c(gene_examples, noisyGenes))

for(cc in conditions)
{
  #cc = c('day2_beforeRA')
  # cc = c('day2.5_RA')
  # cc = c('day3_RA.rep1')
  # cc = c('day3.5_RA')
  # cc = c('day4_RA')
  # cc = c('day5_RA')
  #cc = c('day3_RA.rep1', 'day3.5_RA')
  cc = c('day2_beforeRA', 'day2.5_RA', 'day3_RA.rep1')
  
  outDir_cc = paste0(outDir, '/', paste0(cc, collapse = "_"), '/')
  system(paste0('mkdir -p ', outDir_cc))
  
  aa = subset(bb, cells = colnames(bb)[which(!is.na(match(bb$condition, cc)))])
  
  Idents(aa) = aa$condition
  table(aa$condition)
 
  DimPlot(aa, group.by = 'condition', label = TRUE)
  FeaturePlot(aa, features = c('Pax6', 'Foxa2'))
  
  ##########################################
  # UMAP and clustering without sparse feature selection
  ##########################################
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000) # find subset-specific HVGs
  variableGenes = VariableFeatures(object = aa)
  select.method = 'HVGs2000_pca.weighted_'
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  #variableGenes = setdiff(VariableFeatures(object = aa), 'Lypd2')
  # variableGenes = intersect(variableGenes, c(tfs, sps))
  
  #variableGenes = VariableFeatures(object = aa)
  cat(length(variableGenes), ' variable genes using \n ')
  
  aa <- RunPCA(aa, features = variableGenes, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) + NoLegend()
  p11 = FeaturePlot(aa, features = 'nCount_RNA')
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(outDir_cc, '/UMAP_RA_condition_clustering_cellcyclePhase_', select.method, '.pdf'), 
         width = 12, height = 6)
  
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(outDir_cc, '/featureExamples_originalUMAP.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = aa, 
                           features = gene_examples[!is.na(match(gene_examples, rownames(aa)))])
  
  dev.off()
  
  ##########################################
  # regress out the cell cycle  
  ##########################################
  xx = aa
  xx$CC.Difference <- xx$S.Score - xx$G2M.Score
  #xx <- ScaleData(xx, vars.to.regress = c('nCount_RNA', "CC.Difference"), features = variableGenes)
  
  xx <- ScaleData(xx, vars.to.regress = c('nCount_RNA', "S.Score", "G2M.Score"), 
                  features = variableGenes)
  
  # cell cycle effects strongly mitigated in PCA
  xx <- RunPCA(xx, features = variableGenes, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(xx, ndims = 50)
  
  xx <- RunUMAP(xx, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
  DimPlot(xx, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  
  xx <- FindNeighbors(xx, dims = 1:10)
  xx <- FindClusters(xx, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  p1 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) + NoLegend()
  p11 = FeaturePlot(xx, features = 'nCount_RNA')
  p2 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(outDir_cc, '/UMAP_RA_condition_clustering_cellcyclePhaseRegression_', 
                           select.method, '.pdf'), 
         width = 12, height = 6)
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(outDir_cc, '/featureExamples_UMAP_cellcycleRegression.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = xx, 
                           features = gene_examples[!is.na(match(gene_examples, rownames(aa)))])
  
  dev.off()
  
  rm(xx)
  
  ##########################################
  # Removing cell cycle-related genes in the HVGs
  # original code :
  # https://bioconductor.org/books/3.12/OSCA/cell-cycle-assignment.html#removing-cell-cycle-effects
  ##########################################
  library(scater)
  #require(batchelor)
  
  Idents(aa) = aa$condition
  
  # Identifying the likely cell cycle genes between phases,
  # using an arbitrary threshold of 5%.
  scaledMatrix = GetAssayData(aa, slot = c("scale.data"))
  
  diff <- getVarianceExplained(scaledMatrix, data.frame(phase = aa$Phase))
  diff = data.frame(diff, gene = rownames(diff))
  diff = diff[order(-diff$phase), ]
  
  hist(diff$phase, breaks = 100); abline(v = c(1:5), col = 'red')
  
  genes_discard = diff$gene[which(diff$phase>5)]
  cat(length(genes_discard), 'genes to discard \n')
  
  hvgs = VariableFeatures(FindVariableFeatures(subset(aa, features = setdiff(rownames(aa), genes_discard)), 
                                               selection.method = "vst", nfeatures = 2000))
  # 
  # discard <- diff > diff_cutoff
  # summary(discard)
  # cat('cut off -- ', diff_cutoff, '--', sum(discard), ' genes related to cell cycle \n')
  # genes.discard = rownames(diff)[which(discard)]
  aa <- RunPCA(aa, 
               features = hvgs, 
               verbose = FALSE, 
               weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) + NoLegend()
  p11 = FeaturePlot(aa, features = 'nCount_RNA')
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(outDir_cc, '/UMAP_RA_condition_clustering_cellcycle.rmCellCycleGenes_', 
                           select.method, '.pdf'),  width = 12, height = 6)
  
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(outDir_cc, '/featureExamples_UMAP_cellcycle.rmCellCycleGenes.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = aa, 
                           features = gene_examples[!is.na(match(gene_examples, rownames(aa)))])
  
  dev.off()
  
  # Idents(aa) = aa$seurat_clusters
  # all.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
  # all.markers %>%
  #   group_by(cluster) %>%
  #   top_n(n = 10, wt = avg_log2FC) -> top10
  # 
  # DoHeatmap(aa, features = top10$gene) + NoLegend()
  # ggsave(filename = paste0(outDir_cc, '/Heatmap_clusterMarkers_', select.method, '.pdf'), 
  #        width = 14, height = 20)
  
}

########################################################
########################################################
# Section : test the local gene-gene correlation
# original code from https://github.com/hms-dbmi/crestree/blob/master/R/bifurcation.functions.R
#
########################################################
########################################################
# library(igraph)
# library(mgcv)
# library(quadprog) 
# library(pcaMethods) 
# library(Rcpp) 
# library(inline) 
# library(RcppArmadillo) 
# library(pbapply)
# library(crestree)
# library(ggplot2); library(gridExtra); library(grid);
# library(parallel)
# 
# data(crest)
# emb <- crest$emb
# clcol <- crest$clcol
# nc.cells <- crest$nc.cells
# wgm <- crest$wgm
# wgwm <- crest$wgwm # matrix of expression weights
# fpm <- read.table("http://pklab.med.harvard.edu/ruslan/neural_crest/fpm.txt",header=TRUE)
# fpm <- as.matrix(fpm)
# genes.tree <- crest$genes.tree
# 
# ppt <- readRDS(url("http://pklab.med.harvard.edu/ruslan/neural_crest/tree_structure_full.rds"))
# 
# plotppt(ppt,emb,tips=TRUE,forks=FALSE,cex.tree = 0.2,lwd.tree = 2)
# 
# root <- 355
# leaves <- c(165,91)
# 
# subtree <- extract.subtree(ppt,c(root,leaves))
# plotppt(ppt,emb,tips=TRUE,forks=FALSE,cex.tree = 0.3,lwd.tree = 3,subtree=subtree)
# 
# reRun_forkDE = FALSE
# if(reRun_forkDE){
#   fork.de <- test.fork.genes(ppt,fpm,root=root,leaves=leaves,n.mapping = 5)
#   saveRDS(fork.de, file = paste0(outDir, 'saved_forkDE_crestree.rds' ))
# }else{
#   fork.de = readRDS(file = paste0(outDir, 'saved_forkDE_crestree.rds' ))
# }
# 
# head(fork.de[order(fork.de$p),], )
# 
# fork.de <- branch.specific.genes(fork.de,effect.b1 = 0.1,effect.b2 = 0.3)
# 
# gene <- "Neurog2"
# 
# ppt2 <- fit.associated.genes(ppt,fpm,n.map=1)
# ## [1] "fit gene expression for mapping 1"
# ## 
# ##     branch-monotonous      complex patterns transiently expressed 
# ##                   673                   112                   263
# 
# visualise.trajectory(ppt, gene, fpm[gene,], subtree = subtree, cex.main = 3,lwd.t2=0.5)
# 
# genes.sensory <- rownames(fork.de)[fork.de$state==1]
# genes.autonomic  <- rownames(fork.de)[fork.de$state==2]
# 
# genes.sensory <- intersect(genes.sensory,genes.tree)
# str(genes.sensory)
# ##  chr [1:98] "Rdh10" "Hes6" "Cxcr4" "Nfasc" "5730559C18Rik" "Rgs16" "Nhlh1" "Zbtb18" "Utrn" "Gamt" "Neurod4" "Upp1" ...
# 
# genes.autonomic <- intersect(genes.autonomic,genes.tree)
# str(genes.autonomic)
# ##  chr [1:122] "Serpine2" "Lrrfip1" "Cdh19" "Ralb" "Angptl1" "Pbx1" "Mpz" "Lama4" "Cnn2" "Timp3" "Ascl1" "Elk3" ...
# 
# cells <- rownames(ppt$cell.summary)[ppt$cell.summary$seg %in% extract.subtree(ppt,c(root,leaves))$segs]
# par(mfrow=c(1,2))
# plot(t(programs[c(1,3),cells]),col=ppt$cell.summary[cells,]$color,pch=19,cex=0.5)
# plot(t(programs[c(2,4),cells]),col=ppt$cell.summary[cells,]$color,pch=19,cex=0.5)
# 
# 
# fork.de.act <- activation.fork(ppt,fork.de,fpm,root,leaves,deriv.cutoff = 0.015,n.mapping=10)
# 
# fork.pt(ppt,root,leaves)
# 
# cutoff <- 16.0
# 
# genes.sensory.late <- genes.sensory[fork.de.act[genes.sensory,]$activation > cutoff]
# genes.sensory.early <- setdiff(genes.sensory,genes.sensory.late)
# 
# genes.autonomic.late <- genes.autonomic[fork.de.act[genes.autonomic,]$activation > cutoff]
# genes.autonomic.early <- setdiff(genes.autonomic,genes.autonomic.late)
# 
# cells <- rownames(ppt$cell.summary)[ppt$cell.summary$seg %in% extract.subtree(ppt,c(root,leaves))$segs]
# par(mfrow=c(1,2))
# plot(t(programs[c(1,3), cells]),col=ppt$cell.summary[cells,]$color,pch=19,cex=0.5)
# plot(t(programs[c(2,4), cells]),col=ppt$cell.summary[cells,]$color,pch=19,cex=0.5)
# 
# 
# freq <- slide.cells(ppt,root,leaves,wind=50)
# 
# 
# regions = list(list(7,151,200,1),
#                list(7,101,151,1),
#                list(7,51,100,1),
#                list(7,1,50,1),
#                list(list(6,5,1,2),1,50, -1),
#                list(list(6,5,1,2),51,100, -1),
#                list(5,1,50,1),
#                list(1,1,50,1))
# 
# 
# freq <- slide.cells(ppt, root, leaves, wind=50, regions=regions)
# fig_cells <- fig.cells(emb,freq)
# marrangeGrob(c(fig_cells), ncol=length(fig_cells),nrow=1,top=NA)
# 
# cors <- slide.cors(freq,fpm,genes.sensory.early,genes.autonomic.early)
# 
# fig_cor <- fig.cors(cors,genes.sensory.early,genes.autonomic.early)
# marrangeGrob( c(fig_cells,fig_cor),ncol=length(fig_cells),nrow=2,
#               layout_matrix = matrix(seq_len(2*length(fig_cells)), nrow = 2, 
#                                      ncol = length(fig_cells),byrow=TRUE),
#               top=NA)
# 
# 
# ## Re-estimation average window-specific correlations for cleaned up sets of genes genesetA and genesetB:
# corA <- cors[[5]][,1]
# genesetA <- names(which(corA[genes.sensory.early] > 0.07))
# 
# corB <- cors[[5]][,2]
# genesetB <- names(which(corB[genes.autonomic.early] > 0.07))
# 
# cors <- slide.cors(freq, fpm, genesetA, genesetB)
# fig_cor <- fig.cors(cors, genesetA, genesetB)
# marrangeGrob( c(fig_cells,fig_cor),ncol=length(fig_cells),nrow=2,
#               layout_matrix = matrix(seq_len(2*length(fig_cells)), nrow = 2, 
#                                      ncol = length(fig_cells),byrow=TRUE),top=NA)
# 
# 
