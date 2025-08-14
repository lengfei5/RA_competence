##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: test sparse feature selection methods from scRNA-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Feb 21 11:36:15 2023
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

# levels_sels = c("day3_RA.rep1", "day3.5_RA", "day4_RA")
# data_version = "_d3_d3.5_d4"

#levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1")
#data_version = "_d2_d2.5_d3"

#levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA")
#data_version = "_d2_d2.5_d3_d3.5_d4"

#levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA")
#data_version = "_d2_d2.5_d3_d3.5_d4"

# levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1",  "day3_RA.rep2",
#                 "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")

levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA")
data_version = "_d2_d2.5_d3_d3.5_d4_d5"

names(cols) = levels
cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/RA_symetryBreaking/sparse_featureSelection', data_version)
system(paste0('mkdir -p ', outDir))

##########################################
# tfs and sps annotations 
##########################################
sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

# all 16 samples
# aa = readRDS(file = paste0(RdataDir, 
#                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
#                            'cellCycleScoring_annot.v1_savedUMAP.v1_', species, version.analysis, '.rds'))

##########################################
# import the RA data without day6 and mature neurons 
##########################################
# only RA samples incl. dya2_beforeRA
Import_prepare_RAsamples = FALSE
if(Import_prepare_RAsamples){
  aa = readRDS(file = paste0(RdataDir,
                             'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                             'cellCycleScoring_annot.v1_savedUMAP.subs.v2_', species, version.analysis, '.rds'))
  
  Idents(aa) = aa$condition
  DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)
  
  FeaturePlot(aa, features = c('Zfp42', 'Tcf15', 'Skil', 'Lef1'))
  
  ggsave(filename = paste0(outDir, '/asymmetric_feature_expression.pdf'), 
         width = 14, height = 10)
  
  
  aa = subset(aa, idents = levels_sels)
  cell_sels = colnames(aa)[which(aa$celltypes != 'Neurons'|is.na(aa$celltypes))]
  aa = subset(aa, cells = cell_sels)
  
  aa$condition = droplevels(as.factor(aa$condition))
  Idents(aa) = aa$condition
  DimPlot(aa, cols = cols_sel, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)
  
  # downsample for each time point 
  Idents(aa) = aa$condition
  aa = subset(aa, downsample = 2000)
  
  DimPlot(aa, group.by = 'clusters', label = TRUE)
  
  
  # rerun the umap 
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.05)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  # quickly run clustering
  #ElbowPlot(aa, ndims = 50)
  aa <- FindNeighbors(aa, dims = 1:20)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p1 + p2
  
  ggsave(filename = paste0(outDir, '/UMAP_RA_symmetryBreaking.pdf'), 
         width = 14, height = 6)
  
  #DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)
  aa$clusters = aa$seurat_clusters
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                            'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', data_version, '_',
                            species, version.analysis, '.rds'))
  
  
}

########################################################
########################################################
# Section I: test sparse feature selection methods
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 
              'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
              'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', data_version, '_',
              species, version.analysis, '.rds'))

DimPlot(aa, group.by = 'clusters', label = TRUE, repel = TRUE)
FeaturePlot(aa, features = c('Wnt6'))

FeaturePlot(aa, features = c('Nr2f1', 'Peg3', 'Zeb1', 'Irx3', 'Pax6', 'Dusp5'))

FeaturePlot(aa, features = c('Zfp42', 'Tcf15', 'Skil', 'Lef1'))

ggsave(filename = paste0(outDir, '/asymmetric_feature_expression.pdf'), 
       width = 14, height = 10)

Test_geneBasis = FALSE
if(Test_geneBasis){
  # install scater https://bioconductor.org/packages/release/bioc/html/scater.html
  library(scater)
  # install loomR from GitHub using the remotes package remotes::install_github(repo =
  # 'mojaveazure/loomR', ref = 'develop')
  #library(loomR)
  library(Seurat)
  library(patchwork)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(geneBasisR)
  
  sce = as.SingleCellExperiment(aa)
  sce$cell = colnames(sce)
  sce$celltype = sce$clusters
  
  # sce - SingleCellExperiment object, wehre normalized counts stored in 'logcounts' assay
  # discard definetely uninteresting genes
  sce = retain_informative_genes(sce)
  
  mm = match(rownames(sce), c(tfs, sps))
  tf_keep = rownames(sce)[which(!is.na(mm))]
  #sce = sce[which(!is.na(mm)),]
  
  # run gene selection
  genes = gene_search(sce, 
                      n_genes_total = 50, 
                      genes.discard = setdiff(rownames(sce), tf_keep),
                      nPC.all = 100,
                      n.neigh = 10)
  
  saveRDS(genes, file = paste0(outDir, 'sparse_features_test_geneBasis_50TFs.SPs_v2.rds')) 
  save(genes, sce, file = paste0(RdataDir, 'sparse_features_test_geneBasis_50TFs.SPs_v2.Rdata')) 
  
  genes = gene_search(sce, 
                      n_genes_total = 100, 
                      genes.discard = setdiff(rownames(sce), tf_keep),
                      nPC.all = 100,
                      n.neigh = 10)
  
  saveRDS(genes, file = paste0(outDir, 'sparse_features_test_geneBasis_100TFs.SPs_v2.rds')) 
  save(genes, sce, file = paste0(RdataDir, 'sparse_features_test_geneBasis_100TFs.SPs_v2.Rdata')) 
  
  genes = gene_search(sce, 
                      n_genes_total = 200, 
                      genes.discard = setdiff(rownames(sce), tf_keep),
                      nPC.all = 100,
                      n.neigh = 10)
  saveRDS(genes, file = paste0(outDir, 'sparse_features_test_geneBasis_200TFs.SPs_v2.rds')) 
  save(genes, sce, file = paste0(RdataDir, 'sparse_features_test_geneBasis_200TFs.SPs_v2.Rdata')) 
  
  #saveRDS(genes, file = paste0(outDir, 'sparse_features_test_geneBasis_100TFs.SPs_v2.rds'))
  #saveRDS(genes, file = paste0(outDir, 'sparse_features_test_geneBasis_100TFs_v1.rds')) 
  #load(paste0(RdataDir, 'sparse_features_test_geneBasis_v0.Rdata'))
  
  ##########################################
  # evaluation of gene panels
  ##########################################
  ## reload the v0 with 50TFs
  genes = readRDS(file = paste0(outDir, 'sparse_features_test_geneBasis_50TFs.SPs_v2.rds')) 
  
  # cell scores
  stat = evaluate_library(sce, genes$gene, 
                          genes.all = rownames(sce), 
                          #batch = "sample", 
                          n.neigh = 10,
                          library.size_type = "single", celltype.id = "celltype",
                          return.cell_score_stat = TRUE, 
                          return.gene_score_stat = TRUE, 
                          return.celltype_stat = F, verbose = FALSE)
  
  saveRDS(stat, file = paste0(outDir, 'sparse_features_test_geneBasis_50TFs.SPs_evaluationScores_v2.rds')) 
  
  ##########################################
  # reload the output of geneBasis 
  ##########################################
  load(file = paste0(RdataDir, 'sparse_features_test_geneBasis_100TFs.SPs_v2.Rdata')) 
  write.csv2(genes, file = paste0(outDir, '/geneBasis_top100_tfs_sps_final.csv'), quote = FALSE, row.names = FALSE)
  
  load(file = paste0(RdataDir, 'sparse_features_test_geneBasis_50TFs.SPs_v2.Rdata')) 
  stat = readRDS(file = paste0(outDir, 'sparse_features_test_geneBasis_50TFs.SPs_evaluationScores_v2.rds')) 
  n_genes_total = 50
  
  write.csv2(genes, file = paste0(outDir, '/geneBasis_top50_tfs_sps_final.csv'), quote = FALSE, row.names = FALSE)
  
  
  cell_score_stat = stat$cell_score_stat[stat$cell_score_stat$n_genes == n_genes_total , ] 
  metadata = aa@meta.data
  cell_score_stat = data.frame(cell_score_stat, 
                               metadata[match(rownames(cell_score_stat), rownames(metadata)), ]) 
  
  p = ggplot(cell_score_stat , aes(x = condition , y = cell_score, fill = condition)) + 
    geom_boxplot() + 
    scale_fill_manual(values = cols_sel) + 
    theme_classic() + 
    labs(y = "Cell neighborhood preservation score" , x = "# genes") + 
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  p 
  
  ggsave(filename = paste0(outDir, '/geneBasis_cell.neighborhood.presevartion.scores_ntop.50.pdf'), 
         width = 10, height = 8)
  
  FeaturePlot(aa, features = c('Sox11', 'Sox17', 'Sox4', 'Sox1', 'Sox2'))
  ggsave(filename = paste0(outDir, '/sparse_features_sox.pdf'), 
         width = 10, height = 8)
  
  FeaturePlot(aa, features = c('Pou3f3', 'Pou3f2', 'Pou5f1', 'Meis2'))
  ggsave(filename = paste0(outDir, '/sparse_features_Oct_Meis.pdf'), 
         width = 10, height = 8)
  
  FeaturePlot(aa, features =genes$gene[1:20])
  ggsave(filename = paste0(outDir, '/sparse_features_geneBasis_top20.pdf'), 
         width = 24, height = 24)
  
  FeaturePlot(aa, features =genes$gene)
  ggsave(filename = paste0(outDir, '/sparse_features_geneBasis_top50.pdf'), 
         width = 24, height = 40)
  
  FeaturePlot(aa, features = c('Pou3f3', 'Pou3f2', 'Pou5f1', 'Meis2',  "Foxp1", 
                               'Zfp703', 'Hoxa1', 'Tet2'))
  
  FeaturePlot(aa, features = c('Meis2', 'Zfp703', 'Hoxa1'))
  #genes = genes_stat$gene
  
  p = plot_expression_heatmap(sce, genes = genes$gene, value.type = "mean", celltype.id = 'condition') 
  p 
  ggsave(filename = paste0(outDir, '/sparse_features_geneBasis_heatmap_ntop50.pdf'), 
         width = 10, height = 6)
  
  p = plot_coexpression(sce, genes = genes$gene) 
  p 
  ggsave(filename = paste0(outDir, '/sparse_features_geneBasis_coexpression_ntop50.pdf'), 
         width = 10, height = 6)
  
  ##########################################
  # Let’s compare UMAP plots when using whole transcriptome and 
  # using the selected panel (UMAP-coordinates can be calculated using get_umap_coordinates).
  # For consistency, e will re-calculate UMAP coordinates for the whole transcriptome as well.
  ##########################################
  load(file = paste0(RdataDir, 'sparse_features_test_geneBasis_200TFs.SPs_v2.Rdata'))
  
  write.csv(genes, file = paste0(outDir, '/geneBasis_top200_tfs_sps.csv'), quote = FALSE, row.names = FALSE)
  
  # test seurat umap using selected sparse features
  sub_obj = subset(aa, features = genes$gene[1:50])
  
  sub_obj <- RunPCA(sub_obj, features = rownames(sub_obj), verbose = FALSE, weight.by.var = FALSE,
                    npcs = 50)
  ElbowPlot(sub_obj, ndims = 50)
  
  sub_obj <- RunUMAP(sub_obj, dims = 1:10, n.neighbors = 30, min.dist = 0.1)
  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) 
  
  p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) 
  p1 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  #FeaturePlot(sub_obj, features = c('Foxa2', 'Pax6', 'Meis2'))
  p0 + p1
  
  ggsave(filename = paste0(outDir, '/UMAP_origine_vs.sparseFeatures_geneBasis_ntop.50tfs.sps.pdf'), 
         width = 14, height = 6)
  
  FeaturePlot(sub_obj, features = c('Pax6', 'Foxa2', 'Zfp42', 'Tcf15', 'Skil', 'Lef1', 
                                     'Tfap2c', 'Cebpb', 'Sox17', 'Prdm1', 'Cdh1', 'Irx3', 'Sox1'))
  
  ggsave(filename = paste0(outDir, '/asymmetric_feature_expression_umap_geneBasis_ntop200.tfs.sps.pdf'), 
         width = 14, height = 10)
  
  pdf(paste0(outDir, '/sparse_features_geneBasis_top200.tfs.sps.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  
  for(n in 1:20)
  {
    kk = c(((n-1)*10+1):(n*10))
    cat(n, ' -- ', kk, '\n')
    p = FeaturePlot(sub_obj, features =genes$gene[kk])
    plot(p)                                         
  }
  
  dev.off()
  
}

##########################################
# test DUBStepR
# exmaple code from 
# https://cran.r-project.org/web/packages/DUBStepR/vignettes/dub-vignette.html
##########################################
test_DUBStepR = FALSE
if(test_DUBStepR){
  aa = readRDS(file = paste0(RdataDir, 
                'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', data_version, '_',
                species, version.analysis, '.rds'))
  
  library(DUBStepR)
  
  require(tictoc)
  tic()
  dubstepR.out <- DUBStepR(input.data = aa@assays$RNA@data, 
                           min.cells = 0.05*ncol(aa), 
                           optimise.features = TRUE, 
                           k = 20, num.pcs = 30, 
                           error = 0)
  toc()
  
  saveRDS(dubstepR.out, file = paste0(RdataDir, 'sparse_features_test_dubstepR_v0.rds'))
  
  dubstepR.out = readRDS(file = paste0(RdataDir, 'sparse_features_test_dubstepR_v0.rds'))
  
  seuratObj = aa;
  seuratObj@assays$RNA@var.features <- dubstepR.out$optimal.feature.genes
  seuratObj
  
  genes.dub = dubstepR.out$optimal.feature.genes
  
  genes.basis =  readRDS(file = paste0(outDir, 'sparse_features_test_geneBasis_100TFs.SPs_v2.rds'))
  
  intersect(genes.dub, unique(c(tfs, sps)))
  intersect(genes.dub, genes.basis$gene)
  
  # test seurat umap using selected sparse features
  sub_obj = subset(aa, features = genes.dub)
  sub_obj <- RunPCA(sub_obj, features = rownames(sub_obj), verbose = FALSE, weight.by.var = FALSE, 
                    npcs = 50)
  ElbowPlot(sub_obj, ndims = 50)
  
  sub_obj <- RunUMAP(sub_obj, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
  
  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  ggsave(filename = paste0(outDir, '/sparse_features_dubstepR_448genes.pdf'), 
         width = 10, height = 6) 
  
  FeaturePlot(sub_obj, features = c('Pax6', 'Foxa2', 'Zfp42', 'Tcf15', 'Skil', 'Lef1'))
  
  ggsave(filename = paste0(outDir, '/asymmetric_feature_expression_umap_dubstep_448genes.pdf'), 
         width = 14, height = 10)
  
  
  ## test the selected sparse tfs and sps
  sub_obj = subset(aa, features = intersect(genes.dub, unique(c(tfs, sps))))
  sub_obj <- RunPCA(sub_obj, features = rownames(sub_obj), verbose = FALSE, weight.by.var = FALSE,
                    npcs = 50)
  ElbowPlot(sub_obj, ndims = 50)
  
  
  sub_obj <- RunUMAP(sub_obj, dims = 1:10, n.neighbors = 30, min.dist = 0.1
                     #metric = "euclidean"
                     )
  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  ggsave(filename = paste0(outDir, '/sparse_features_dubstepR_74.tfs.spss.pdf'), 
         width = 10, height = 6)
  
  #FeaturePlot(sub_obj, features = rownames(sub_obj))
  FeaturePlot(sub_obj, features = c('Pax6', 'Foxa2', 'Zfp42', 'Tcf15', 'Skil', 'Lef1'))
  ggsave(filename = paste0(outDir, '/asymmetric_feature_expression_umap_dubstep.74.tfs.sps.pdf'), 
         width = 14, height = 10)
  
  write.csv(genes.dub, file = paste0(outDir, '/dubstep_sparse_448genes.csv'), row.names = FALSE,
             quote = FALSE)
  
  write.csv2(intersect(genes.dub, unique(c(tfs, sps))), 
            file = paste0(outDir, '/dubstep_sparse_74.tfs.sps_final.csv'), row.names = FALSE, col.names = FALSE,
            quote = FALSE)
  
  ggs = intersect(genes.dub, unique(c(tfs, sps)))
  pdf(paste0(outDir, '/sparse_features_dubStep_top74.tfs.sps.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  
  for(n in 1:8)
  {
    kk = c(((n-1)*10+1):(n*10))
    cat(n, ' -- ', kk, '\n')
    p = FeaturePlot(sub_obj, features =ggs[kk])
    plot(p)                                         
  }
  
  dev.off()
  
}

##########################################
# test SMD
##########################################
Test_SMD = FALSE
if(Test_SMD){
  reDownsample_cells = FALSE
  
  if(reDownsample_cells){
    # only RA samples incl. dya2_beforeRA
    aa = readRDS(file = 
                   paste0(RdataDir, 
                          'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', data_version, '_',
                            species, version.analysis, '.rds'))
    # aa = subset(aa, idents = levels_sels)
    # cell_sels = colnames(aa)[which(aa$celltypes != 'Neurons'|is.na(aa$celltypes))]
    # aa = subset(aa, cells = cell_sels)
    
    aa$condition = droplevels(as.factor(aa$condition))
    Idents(aa) = aa$condition
    DimPlot(aa, cols = cols_sel, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)
    
    # downsample for each time point 
    Idents(aa) = aa$condition
    aa = subset(aa, downsample = 3000)
    
    DimPlot(aa, group.by = 'condition', label = TRUE)
    
    saveRDS(aa, file = paste0(RdataDir, 
                              'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                              'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures_SMD', data_version, '_',
                              species, version.analysis, '.rds'))
    
  }
  
  aa = readRDS(file = 
                 paste0(RdataDir, 
                        'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                        'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', data_version, '_',
                        species, version.analysis, '.rds'))
  Idents(aa) = aa$condition
  #aa = subset(aa, downsample = 2000)
  
  DimPlot(aa, group.by = 'condition', repel = TRUE)
  
  exp.mat = aa@assays$RNA@data
  mm = match(rownames(exp.mat), unique(c(tfs, sps)))
  
  exp.mat = exp.mat[which(!is.na(mm)), ]
  
  cat(nrow(exp.mat), 'tfs detected \n')
  
  means = apply(exp.mat, 1, mean)
  sds = apply(exp.mat, 1, sd)
  
  # plot(means, sds)
  
  exp.mat = exp.mat[means > 0.05 & sds > 0.05,  ]
   
  cat(nrow(exp.mat), 'tfs after filtering \n')
  
  exp.mat = t(apply(exp.mat, 1, scale, center = TRUE, scale = TRUE))
  exp.mat = as.matrix(t(exp.mat))
  
  write.csv(exp.mat, file = paste0(outDir, '/exp_matrix_TFs_SPs_4SMD_12k.cells.csv'), 
            quote = FALSE, row.names = FALSE)
  
  
  ##########################################
  # process the SMD output 
  ##########################################
  smd = read.csv(file = paste0(outDir, '/output_SMD_12k.cells_v3.csv'), sep = '\t')
  #smd = read.csv(file = paste0(outDir, '/output_SMD_12k.cells_tfs.sps_v4.csv'), sep = '\t')
  smd = smd[, -1]
  smd = smd[order(-smd$SMD_z), ]
  
  saveRDS(smd, file = paste0(outDir, '/output_SMD_12k.cells_tfs_v3.rds'))
  saveRDS(smd, file = paste0(outDir, '/output_SMD_12k.cells_tfs.sps_v4.rds'))
  
  ## save the final list of smd
  Save_final.list_SMD = FALSE
  if(Save_final.list_SMD){
    
    smd = readRDS(file = paste0(outDir, '/output_SMD_12k.cells_tfs_v3.rds'))
    genes.smd = smd$gene[which(smd$SMD_z>0.5)]
    
    write.csv(genes.dub, file = paste0(outDir, '/SMD_sparse_21tfs_final.csv'), row.names = FALSE,
              quote = FALSE)
    
    smd = readRDS(file = paste0(outDir, '/output_SMD_12k.cells_tfs.sps_v4.rds'))
    genes.smd = smd$gene[which(smd$SMD_z>0.5)]
    
    write.csv(genes.dub, file = paste0(outDir, '/SMD_sparse_31tfs.sps_final.csv'), row.names = FALSE,
              quote = FALSE)
    
  }
  
  
  
  
  smd = readRDS(file = paste0(outDir, '/output_SMD_12k.cells_tfs.sps_v4.rds'))
  smd = smd[order(-smd$SMD_z), ]
  FeaturePlot(aa, features = smd$gene[which(smd$SMD_z>0.5)])
  
  ggsave(filename = paste0(outDir, '/sparse_features_SMD_31TFs_SPs.pdf'), 
         width = 24, height = 18) 
  
  ## test umap with SMD features
  genes.smd = smd$gene[which(smd$SMD_z>0.5)]
  cat(length(genes.smd), 'feature selected by smd \n')
  
  
  sub_obj = subset(aa, features = genes.smd)
  sub_obj <- RunPCA(sub_obj, features = rownames(sub_obj), verbose = FALSE, weight.by.var = FALSE, 
                    npcs = 30)
  ElbowPlot(sub_obj, ndims = 30)
  
  sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) +
    NoLegend()
  p1 =  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)+
    NoLegend()
  
  p0 + p1
  
  ggsave(filename = paste0(outDir, '/Umap_comparison_sparseFeatures_smd_31.tfs.sps.pdf'), 
         width = 12, height = 6) 
  
  FeaturePlot(sub_obj, features = smd$gene[which(smd$SMD_z>0.6)], max.cutoff = 'q1')
  
  ggsave(filename = paste0(outDir, '/sparse_features_umap.smd.features_SMD_31TFs_SPs.pdf'), 
         width = 24, height = 18) 
  
  # compare the features from geneBasis and dupstep
  Compare.selected.features.with.dubstep = FALSE
  if(Compare.selected.features.with.dubstep){
    dubstepR.out = readRDS(file = paste0(RdataDir, 'sparse_features_test_dubstepR_v0.rds'))
    genes.dub = dubstepR.out$optimal.feature.genes
    
    genes.smd = smd$gene[which(smd$SMD_z>1.0)]
    mm = match(genes.smd, genes.dub)
    genes.smd[which(is.na(mm))]
    
    genes.basis =  readRDS(file = paste0(outDir, 'sparse_features_test_geneBasis_100TFs.SPs_v2.rds'))
    
    mm = match(genes.smd, genes.basis$gene)
    genes.smd[which(is.na(mm))]
    
  }
  
}

########################################################
########################################################
# Section : Compare different methods: geneBasis, DUBStepR and SMD  
# merge the gene list from DUBStepR, geneBasis and SMD
########################################################
########################################################
Combine_sparse_featureSelection = FALSE
if(Combine_sparse_featureSelection){
  
  dubstepR.out = readRDS(file = paste0(RdataDir, 'sparse_features_test_dubstepR_v0.rds'))
  genes.dub = intersect(dubstepR.out$optimal.feature.genes, unique(c(tfs, sps)))
  
  genes.basis =  readRDS(file = paste0(outDir, 'sparse_features_test_geneBasis_50TFs.SPs_v2.rds'))
  genes.basis = genes.basis$gene
  
  xx.basis = readRDS(file = paste0(outDir, 'sparse_features_test_geneBasis_50TFs_v0.rds'))
  genes.basis = unique(c(genes.basis, xx.basis$gene))
    
  genes.smd = readRDS(file = paste0(outDir, '/output_SMD_12k.cells_tfs.sps_v4.rds'))
  genes.smd = genes.smd$gene[which(genes.smd$SMD_z>0.6)]
  
  #ggs = unique(c(genes.dub, genes.basis, genes.smd))
  ggs.addtions =  c('Pax6', 'Foxa2', 'Zfp42', 'Tcf15', 'Skil', 'Lef1', 
                 'Tfap2c', 'Cebpb', 'Sox17', 'Prdm1', 'Cdh1', 'Irx3', 'Sox1')
  
  source('functions_utility.R')
  ggs.merged = merge_sparseFeatures(list(genes.dub, genes.basis, genes.smd, ggs.addtions), manual_addtion = TRUE)
  FeaturePlot(aa, features = rownames(ggs.merged)[c(1:15)])
  
  ggsave(filename = paste0(outDir, '/originalUMAP_commonSparseFeatures.pdf'), 
         width = 14, height = 10)
  
  FeaturePlot(aa, features = c('Zfp42', 'Fgf4', 'Tcf15'))
  ##########################################
  # save the results for Hannah and Elnea 
  ##########################################
  saveDir = paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/scRNAseq_mNT/',
                   'RA_symetryBreaking/sparse_featureSelection_d2_d2.5_d3_d3.5_d4_d5/merged_table_plots/')
  
  source('functions_utility.R')
  ggs.merged = merge_sparseFeatures(list(genes.dub, genes.basis, genes.smd, ggs.addtions), manual_addtion = TRUE)
  
  saveRDS(ggs.merged, file = paste0(RdataDir, 'sparseFeatures_merged_dubstepR.geneBasis.SMD_v2.rds'))
  
  # test seurat umap using selected sparse features
  ggs = rownames(ggs.merged)
  # ggs = unique(ggs, res$genes)
  sub_obj = subset(aa, features = ggs)
  sub_obj <- RunPCA(sub_obj, features = rownames(sub_obj), verbose = FALSE, weight.by.var = FALSE, 
                    npcs = 50)
  ElbowPlot(sub_obj, ndims = 50)
  
  sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 20, min.dist = 0.1)
  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  sub_obj <- FindNeighbors(sub_obj, dims = 1:20)
  sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  p1 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p11 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(saveDir, '/updatedUMAP_condition_clustering_cellcyclePhase_all.pdf'), 
         width = 20, height = 12)
  
  write.csv2(gg.merged, file = paste0(saveDir, 'merged_sparse_features_withRanks.csv'), 
             row.names = TRUE, quote = FALSE)
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(saveDir, '/updatedUMAP_merged_sparseFeatures_all.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = sub_obj, features = ggs)
  
  dev.off()
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p11 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(saveDir, '/originalUMAP_condition_clustering_cellcyclePhase_all.pdf'), 
         width = 20, height = 12)
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(saveDir, '/originalUMAP_merged_sparseFeatures_all.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = aa, features = ggs)
  
  dev.off()
  
  
  p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) +
    NoLegend()
  p1 =  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)+
    NoLegend()
  
  p0 + p1
  ggsave(filename = paste0(outDir, '/Umap_comparison_origine_sparseFeatures_dubstepR.geneBasis.SMD.merged.pdf'), 
         width = 10, height = 6) 
  
  FeaturePlot(sub_obj, features = c('Pax6', 'Foxa2', 'Zfp42', 'Tcf15', 'Skil', 'Lef1', 
                                    'Tfap2c', 'Cebpb', 'Sox17', 'Prdm1', 'Cdh1', 'Irx3', 'Sox1'))
  ggsave(filename = paste0(outDir, 
                           '/check_representiveFeature_umapEmbedding_mergedFeatures_dubstepR.geneBasis.SMD.pdf'), 
         width = 14, height = 10)
  
  
  
  ##########################################
  # visualize the BGP result 
  ##########################################
  res = readRDS(file = paste0('../results/Rdata/', 
                              'symmetry_breaking_early.tfs.sps_BGP_output.rds'))
  res = res[which(res$bf>5), ]
  
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(saveDir, '/originalUMAP_BranchingGenes_BGP.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = aa, features = res$gene)
  
  dev.off()
  
  write.csv2(res, file = paste0(saveDir, 'BGP_branchingGenes_branchingPseudotime_bayesianFactor.csv'), 
             row.names = FALSE, quote = FALSE)
  
  
}


########################################################
########################################################
# Section : test features of RA signaling pathways and RA targets
# 
########################################################
########################################################
aa = readRDS(file = 
               paste0(RdataDir, 
                      'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                      'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', data_version, '_',
                      species, version.analysis, '.rds'))
# aa = subset(aa, idents = levels_sels)

sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
genes.merged = readRDS(file = paste0(RdataDir, 'sparseFeatures_merged_dubstepR.geneBasis.SMD_v1.rds'))

ggs = sps$gene[which(sps$pathway == 'RA')]
sps = unique(sps$gene)

# manually adding genes from reactome https://reactome.org/PathwayBrowser/#/R-HSA-5362517&DTAB=MT
ra = read.delim('../data/RA_reactome.tsv')
ra = ra$MoleculeName
ra = as.character(sapply(ra, function(x) unlist(strsplit(as.character(x), ' '))[2]))
ra = firstup(ra)
ggs = unique(c(ggs, ra))

ggs = c(ggs, genes.merged)

ggs = intersect(ggs, rownames(aa))

sub_obj = subset(aa, features = ggs)
sub_obj <- RunPCA(sub_obj, features = rownames(sub_obj), verbose = FALSE, weight.by.var = FALSE, 
                  npcs = 50)
ElbowPlot(sub_obj, ndims = 50)

sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 20, min.dist = 0.1)
DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p1 =  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

p0 + p1

ggsave(filename = 
         paste0(outDir, '/Umap_comparison_origine_sparseFeatures_merge.sparseFeature.plus.RA.pathways.pdf'), 
       width = 10, height = 6) 

FeaturePlot(sub_obj, features = c('Pax6', 'Foxa2', 'Zfp42', 'Tcf15', 'Skil', 'Lef1', 
                                  'Tfap2c', 'Cebpb', 'Sox17', 'Prdm1', 'Cdh1', 'Irx3', 'Sox1',
                                  "Cyp26a1", "Cyp26b1"))
ggsave(filename = paste0(outDir, 
                         '/check_representiveFeature_umapEmbedding_mergedFeatures_RA.pathways.pdf'), 
       width = 14, height = 10)


## first found the two split clusters of day2.5_RA with 20 PCs and won't see them wiht 10 PCs
## now test the 30 PCs
sub_obj <- RunUMAP(sub_obj, dims = 1:30, n.neighbors = 20, min.dist = 0.1)
DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

ggsave(filename = 
         paste0(outDir, '/Umap_condition_clusters_merge.sparseFeature.plus.RA.pathways_30PCs_v2.pdf'), 
       width = 16, height = 6) 

sub_obj <- FindNeighbors(sub_obj, dims = 1:30)
sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.7)

p1 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

ggsave(filename = 
         paste0(outDir, '/Umap_merge.sparseFeature.plus.RA.pathways_30PCs_v2.pdf'), 
       width = 10, height = 6) 

Idents(sub_obj) = sub_obj$seurat_clusters

cluster3.markers <- FindMarkers(sub_obj, only.pos = TRUE, ident.1 = 3, ident.2 = 4, 
                                min.pct = 0.1, logfc.threshold = 0.1)
head(cluster3.markers, n = 10)

FeaturePlot(sub_obj, features = c(rownames(cluster3.markers)[1:5],'Foxa2', 'Pax6'))
ggsave(filename = paste0(outDir, '/FeaturePlots_markers_cluster4_30PCs_v2.pdf'), 
       width = 14, height = 8)

cluster4.markers <- FindMarkers(sub_obj, only.pos = TRUE, ident.1 =4 , ident.2 = 3, 
                                min.pct = 0.1, logfc.threshold = 0.1)
head(cluster4.markers, n = 10)

FeaturePlot(sub_obj, features = c(rownames(cluster4.markers)[1:10], 'Foxa2', 'Pax6'), max.cutoff = 'q95')
ggsave(filename = paste0(outDir, '/FeaturePlots_markers_cluster4_30PCs_v2.pdf'), 
       width = 14, height = 8)

FeaturePlot(sub_obj, features = c('Rarg', 'Cyp26a1', 'Dhrs3', 'Foxa2', 'Pax6'), max.cutoff = 'q99')

source('functions_utility.R')

features = c(ra, 'Foxa2', 'Pax6')
pdf(paste0(outDir, '/FeaturePlots_RAgenes_30PCs_v2.pdf'),
    width =16, height = 8, useDingbats = FALSE)

plot_manyFeatures_seurat(seurat_obj = sub_obj, features = features)

dev.off()


##########################################
# import the RA targets 
##########################################
RARtargetDir = '../results/RA_targets_L118404_smartseq3_20221117/Compare.diffRAstimulationTime.Batch1/'
# peaks = readRDS(file = paste0(RARtargetDir, 
#                               'RAR_chipseq.peak.assignment_promoters.genebody.downstream.chiapet.closestTSS.rds'))
# 
# ggs = unlist(peaks[, c(13:17)])
# ggs = ggs[!is.na(ggs)]
# ggs = unique(ggs)
ggs = read.csv2(paste0(RARtargetDir, 
                       'RARtarget_intersection/DESeq2_DEgenes_pairwiseComparison_RA_d2.18h.vs.noRA_d2.18h.csv'),
               row.names = c(1))
ggs = rownames(ggs)[which(abs(ggs$log2FoldChange_RA_d2.18hvs.noRA_d2.18h)>2 & 
                      ggs$padj_RA_d2.18hvs.noRA_d2.18h<0.01)]

saveRDS(ggs, file = paste0(RdataDir, '/RA_targets_smartseq2_chipseq_strigent.rds'))
genes.merged = readRDS(file = paste0(RdataDir, 'sparseFeatures_merged_dubstepR.geneBasis.SMD_v1.rds'))

ggs = unique(c(ggs, genes.merged))

ggs = intersect(ggs, rownames(aa))
ggs = c(ggs, 'Dhrs3')

sub_obj = subset(aa, features = ggs)
sub_obj <- RunPCA(sub_obj, features = rownames(sub_obj), verbose = FALSE, weight.by.var = FALSE, 
                  npcs = 50)
ElbowPlot(sub_obj, ndims = 50)

sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 20, min.dist = 0.1)
DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p1 =  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

p0 + p1

ggsave(filename = 
         paste0(outDir, '/Umap_comparison_origine_selectedFeatures_RAtargets_chipseq_sparseMerged_v2.pdf'), 
       width = 16, height = 6) 

FeaturePlot(sub_obj, features = c('Rarg', 'Cyp26a1', 'Dhrs3', 'Foxa2', 'Pax6'), max.cutoff = 'q99')

ggsave(filename = 
         paste0(outDir, '/RAgenes_selectedFeatures_RAtargets_chipseq_sparseMerged_v2.pdf'), 
       width = 12, height = 8) 


p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
p1 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
p0 + p1

ggsave(filename = 
         paste0(outDir, '/Umap_comparison_origine_selectedFeatures_RAtargets_chipseq_sparseMerged',
                '_clusters_v2.pdf'), 
       width = 16, height = 6) 

##########################################
# check the genes separate clusters in day2.5
##########################################
sub_obj <- FindNeighbors(sub_obj, dims = 1:20)
sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.7)

p1 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

ggsave(filename = 
         paste0(outDir, '/Umap_merge.sparseFeature.plus.RAtargets_20PCs_v2.pdf'), 
       width = 10, height = 6) 

p1 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
p1 + p2

Idents(sub_obj) = sub_obj$seurat_clusters

cluster3.markers <- FindMarkers(sub_obj, only.pos = TRUE, ident.1 = 3, ident.2 = c(7,8), 
                                min.pct = 0.2, logfc.threshold = 0.2)
head(cluster3.markers, n = 10)

FeaturePlot(sub_obj, features = c(rownames(cluster3.markers)[1:10],'Foxa2', 'Pax6'))
ggsave(filename = paste0(outDir, '/RAtargets_chipseq_sparseMerged_FeaturePlots_markers_cluster3_20PCs_v2.pdf'), 
       width = 14, height = 8)

cluster7.markers <- FindMarkers(sub_obj, only.pos = TRUE, ident.1 =7 , ident.2 = c(3, 8), 
                                min.pct = 0.2, logfc.threshold = 0.2)
head(cluster7.markers, n = 10)

FeaturePlot(sub_obj, features = c(rownames(cluster7.markers)[1:10], 'Foxa2', 'Pax6'), max.cutoff = 'q99')
ggsave(filename = paste0(outDir, '/RAtargets_chipseq_sparseMerged_FeaturePlots_markers_cluster7_20PCs_v2.pdf'), 
       width = 14, height = 8)

cluster8.markers <- FindMarkers(sub_obj, only.pos = TRUE, ident.1 =8 , ident.2 = c(3, 7), 
                                min.pct = 0.2, logfc.threshold = 0.2)
head(cluster8.markers, n = 10)

FeaturePlot(sub_obj, features = c(rownames(cluster8.markers)[1:10], 'Foxa2', 'Pax6'), max.cutoff = 'q99')
ggsave(filename = paste0(outDir, '/RAtargets_chipseq_sparseMerged_FeaturePlots_markers_cluster8_20PCs_v2.pdf'), 
       width = 14, height = 8)
