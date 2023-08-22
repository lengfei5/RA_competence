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

levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA")

data_version = "_timePoint"

names(cols) = levels
cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/RA_symetryBreaking/sparse_featureSelection', data_version)
system(paste0('mkdir -p ', outDir))

# tfs and sps annotations 
sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

source('functions_utility.R')

##########################################
# import the RA data without day6 and mature neurons
# wihtout downsampling 
##########################################
Select_samples_processing = FALSE
if(Select_samples_processing){
  aa = readRDS(file = paste0(RdataDir,
                             'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                             'cellCycleScoring_annot.v1_savedUMAP.subs.v2_', species, version.analysis, '.rds'))
  
  Idents(aa) = aa$condition
  DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)
  
  FeaturePlot(aa, features = c('Zfp42', 'Tcf15', 'Skil', 'Lef1'))
  FeaturePlot(aa, features = c('Lhx1', 'Foxa2', 'Pax6', 'Eomes'))
  
  ggsave(filename = paste0(outDir, '/geneExamples_feature_expression.pdf'), 
         width = 14, height = 10)
  
  
  aa = subset(aa, idents = levels_sels)
  cell_sels = colnames(aa)[which(aa$celltypes != 'Neurons'|is.na(aa$celltypes))]
  aa = subset(aa, cells = cell_sels)
  
  aa$condition = droplevels(as.factor(aa$condition))
  Idents(aa) = aa$condition
  DimPlot(aa, cols = cols_sel, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)
  
  # downsample for each time point 
  Idents(aa) = aa$condition
  
  #aa = subset(aa, downsample = 2000)
  #DimPlot(aa, group.by = 'clusters', label = TRUE)
  
  
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
  
  ggsave(filename = paste0(outDir, '/UMAP_RA_symmetryBreaking_condition_clusterss.pdf'), 
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
# Section : test sparse feature selection methods
# 
########################################################
########################################################
bb = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', data_version, '_',
                           species, version.analysis, '.rds'))

# bb  = aa
p1 = DimPlot(bb, group.by = 'condition', label = TRUE, repel = TRUE)
p2 = DimPlot(bb, group.by = 'clusters', label = TRUE, repel = TRUE)

p1 + p2

conditions = unique(bb$condition)
print(conditions)

for(cc in conditions){
  cc = c('day2_beforeRA')
  cc = c('day2.5_RA')
  cc = c('day3_RA.rep1')
  # cc = c('day2_RA', 'day2.5_RA')
  #cc = c('day2.5_RA', 'day3_RA.rep1')
  cc = c('day3.5_RA')
  cc = c('day4_RA')
  cc = c('day5_RA')
  cc = c('day4_RA', 'day5_RA')
  
  # cc = c('day3.5_RA', 'day4_RA')
  cc = c('day3.5_RA', 'day4_RA', 'day5_RA')
  
  #cc = c('day2_beforeRA', 'day2.5_RA', 'day3_RA.rep1')
  #cc = c('day3_RA.rep1', 'day3.5_RA')
  
  # start to remove cluster 9 and cluster 10
  cc = c('day3_RA.rep1', 'day3.5_RA', 'day4_RA', 'day5_RA') # from day3 we start to see Foxa2 and Pax6 clusters
  #cc = c('day3_RA.rep1', 'day3.5_RA', 'day4_RA')  # test end point is day4, found the trajctory became unclear
  
  cc = c('day2.5_RA', 'day3_RA.rep1', 'day3.5_RA', 'day4_RA', 'day5_RA')
  
  cc = c('day2_beforeRA', 'day2.5_RA', 'day3_RA.rep1', 'day3.5_RA', 'day4_RA', 'day5_RA')
  
  Remove_cluster9_cluster10 = FALSE
  if(Remove_cluster9_cluster10){
    
    outDir_cc = paste0(outDir, '/', paste0(cc, collapse = "_"), '_rmCluster9.10', '/')
    system(paste0('mkdir -p ', outDir_cc))
    aa = subset(bb, cells = colnames(bb)[which(!is.na(match(bb$condition, cc))
                                               & as.character(bb$clusters) != '9' 
                                               & as.character(bb$clusters) != '10')])
    
  }else{
    
    outDir_cc = paste0(outDir, '/', paste0(cc, collapse = "_"), '/')
    aa = subset(bb, cells = colnames(bb)[which(!is.na(match(bb$condition, cc)))])
    system(paste0('mkdir -p ', outDir_cc))
    
  }
                                       
  # downsample for each time point
  aa$condition = droplevels(aa$condition)
  Idents(aa) = aa$condition
  table(aa$condition)
  
  if(length(cc)<=3){
    aa = subset(aa, downsample = 3000)
  }else{
    aa = subset(aa, downsample = 2000)
  }
  
  table(aa$condition)
  
  DimPlot(aa, group.by = 'clusters', label = TRUE)
  DimPlot(aa, group.by = 'condition', label = TRUE)
  
  # rm(bb)
  
  FeaturePlot(aa, features = c('Pax6', 'Foxa2'))
  
  ##########################################
  # UMAP and clustering without sparse feature selection
  ##########################################
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000) # find subset-specific HVGs
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  select.method = 'rm.cluster9.cluster10_variableGenes.tfs'
  variableGenes = setdiff(VariableFeatures(object = aa), 'Lypd2')
  # variableGenes = intersect(variableGenes, c(tfs, sps))
  
  #variableGenes = VariableFeatures(object = aa)
  cat(length(variableGenes), ' variable genes using \n ')
  
  aa <- RunPCA(aa, features = variableGenes, verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 20, min.dist = 0.05)
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  FeaturePlot(aa, features = c('Pax6', 'Foxa2', 'Cyp26a1', 'Dhrs3'))
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p11 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p3 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
  (p1 + p11) / (p2 + p3) 
  
  ggsave(filename = paste0(outDir_cc, '/UMAP_RA_condition_clustering_cellcyclePhase_', select.method, '.pdf'), 
         width = 20, height = 12)
  
  
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
                      c('Zfp42', 'Tcf15', 'Skil', 'Lef1', 'Dhrs3', 'Rarg', 'Cyp26a1')
      
  ))
  
  if(!is.null(dev.list())) dev.off()
  pdf(paste0(outDir_cc, '/featureExamples_originalUMAP.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  plot_manyFeatures_seurat(seurat_obj = aa, features = gene_examples)
  
  dev.off()
  
  Idents(aa) = aa$seurat_clusters
  all.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
  all.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  
  DoHeatmap(aa, features = top10$gene) + NoLegend()
  ggsave(filename = paste0(outDir_cc, '/Heatmap_clusterMarkers_', select.method, '.pdf'), 
         width = 14, height = 20)
  
  #FeaturePlot(aa, features = 'Lypd2')
  #ggsave(filename = paste0(outDir_cc, '/featurePlot_cluster_markers.pdf'), width = 10, height = 8)
  
  ##########################################
  # geneBasis to select features
  ##########################################
  Test_geneBasis = TRUE
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
    ntop_sparseFeatures = 50
    
    sce = as.SingleCellExperiment(aa)
    sce$cell = colnames(sce)
    sce$celltype = sce$clusters
    
    # sce - SingleCellExperiment object, wehre normalized counts stored in 'logcounts' assay
    # discard definetely uninteresting genes
    sce = retain_informative_genes(sce)
    
    mm = match(rownames(sce), c(tfs, sps))
    tf_keep = rownames(sce)[which(!is.na(mm))]
    
    # run gene selection
    genes = gene_search(sce, 
                        n_genes_total = ntop_sparseFeatures, 
                        genes.discard = setdiff(rownames(sce), tf_keep),
                        nPC.all = 50,
                        n.neigh = 10)
    
    saveRDS(genes, file = paste0(outDir_cc, 
                                 'sparse_features_test_geneBasis_TFs.SPs_ntop', ntop_sparseFeatures, '.rds')) 
    save(genes, sce, file = paste0(outDir_cc, 'sparse_features_test_geneBasis_TFs.SPs_ntop', 
                                   ntop_sparseFeatures, '.Rdata')) 
    
    genes = readRDS(file = paste0(outDir_cc, 
                                  'sparse_features_test_geneBasis_TFs.SPs_ntop', ntop_sparseFeatures, '.rds')) 
    
    rm(sce)
    
    if(!is.null(dev.list())) dev.off()
    source('functions_utility.R')
    pdf(paste0(outDir_cc, '/sparse_features_geneBasis.tfs.sps_originalUMAP_', 
               ntop_sparseFeatures, '.pdf'),
        width =16, height = 10, useDingbats = FALSE)
    plot_manyFeatures_seurat(seurat_obj = aa, features = genes$gene)
    
    dev.off()
    
    write.csv2(genes, file = paste0(outDir_cc, 'geneBasis_tfs_sps', 
                                   ntop_sparseFeatures, '.csv'), quote = FALSE, row.names = FALSE)
    
    ##########################################
    # Let’s compare UMAP plots when using whole transcriptome and 
    # using the selected panel (UMAP-coordinates can be calculated using get_umap_coordinates).
    # For consistency, e will re-calculate UMAP coordinates for the whole transcriptome as well.
    ##########################################
    # test seurat umap using selected sparse features
    Test_updatedUMAP_sparseFeatures = FALSE
    if(Test_updatedUMAP_sparseFeatures){
      sub_obj <- RunPCA(aa, features = genes$gene, verbose = FALSE, weight.by.var = FALSE,
                        npcs = min(c(25, length(genes$gene)-1)))
      ElbowPlot(sub_obj, ndims = 25)
      
      #Idents(aa) = aa$condition
      sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 20, min.dist = 0.1)
      DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE) 
      
      sub_obj <- FindNeighbors(sub_obj, dims = 1:20)
      sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.5)
      
      p0 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE) 
      p1 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
      
      p0 + p1
      ggsave(filename = paste0(outDir_cc, '/UMAP_origine_vs.sparseFeatures_geneBasis_tfs.sps',
                               ntop_sparseFeatures, '.pdf'), 
             width = 14, height = 6)
      
      p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
      p3 = FeaturePlot(sub_obj, features = c('Pax6', 'Foxa2'))
      
      p2 / p3
      
      ggsave(filename = paste0(outDir_cc, '/UMAP.sparseFeatures_geneBasis_Pax6.vs.Foxa2.sps',
                               ntop_sparseFeatures, '.pdf'), 
             width = 14, height = 6)
      
      source('functions_utility.R')
      
      if(!is.null(dev.list())) dev.off()
      pdf(paste0(outDir_cc, '/sparse_features_geneBasis.tfs.sps_geneBasisUMAP_', ntop_sparseFeatures, '.pdf'),
          width =16, height = 10, useDingbats = FALSE)
      plot_manyFeatures_seurat(seurat_obj = sub_obj, features = genes$gene)
      
      dev.off()
      
    }
    
  }
  
  ##########################################
  # test DUBStepR
  # exmaple code from 
  # https://cran.r-project.org/web/packages/DUBStepR/vignettes/dub-vignette.html
  ##########################################
  test_DUBStepR = TRUE
  if(test_DUBStepR){
        
    library(DUBStepR)
    require(tictoc)
    
    tic()
    dubstepR.out <- DUBStepR(input.data = aa@assays$RNA@data, 
                             min.cells = 0.05*ncol(aa), 
                             optimise.features = TRUE, 
                             k = 20, num.pcs = 30, 
                             error = 0)
    toc()
    
    saveRDS(dubstepR.out, file = paste0(outDir_cc, 'sparse_features_dubstepR_test.v2.rds'))
    
    dubstepR.out = readRDS(file = paste0(outDir_cc, 'sparse_features_dubstepR_test.v2.rds'))
    genes.dub = dubstepR.out$optimal.feature.genes
    
    length(genes.dub)
    length(intersect(genes.dub, unique(c(tfs, sps))))
    
    #intersect(genes.dub, genes.basis$gene)
    tfs.dub = intersect(genes.dub, unique(c(tfs, sps)))
    cat(length(tfs.dub), ' tfs and sps \n')
    
    if(length(tfs.dub)>50) tfs.dub = tfs.dub[c(1:50)]
    
    if(!is.null(dev.list())) dev.off()
    source('functions_utility.R')
    pdf(paste0(outDir_cc, '/sparse_features_dubstepR.tfs.sps_originalUMAP.pdf'),
        width =16, height = 10, useDingbats = FALSE)
    plot_manyFeatures_seurat(seurat_obj = aa, features = tfs.dub)
    
    dev.off()
    
    write.csv2(genes.dub, file = paste0(outDir_cc, '/dubstep_sparse_genes.csv'), row.names = FALSE,
              quote = FALSE)
    write.csv2(tfs.dub, 
              file = paste0(outDir_cc, '/dubstep_sparse_tfs.sps.csv'), row.names = FALSE,
              quote = FALSE)
    
    # test seurat umap using selected sparse features
    Test_updatedUMAP_sparseFeatures = FALSE
    if(Test_updatedUMAP_sparseFeatures){
      # tfs.dub = tfs.dub[1:100]
      sub_obj <- RunPCA(aa, features = tfs.dub, verbose = FALSE, weight.by.var = FALSE, 
                        npcs = 50)
      ElbowPlot(sub_obj, ndims = 50)
      
      sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
      DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
      
      p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
      p3 = FeaturePlot(sub_obj, features = c('Pax6', 'Foxa2'))
      p2 / p3
      
      ggsave(filename = paste0(outDir_cc, '/UMAP.sparseFeatures_geneBasis_Pax6.vs.Foxa2.sps',
                               ntop_sparseFeatures, '.pdf'), 
             width = 14, height = 6)
      
      FeaturePlot(sub_obj, features = c('Pax6', 'Foxa2', 'Zfp42', 'Tcf15', 'Skil', 'Lef1'))
      ggsave(filename = paste0(outDir_cc, '/asymmetric_feature_expression_UMAP.dupstep.features.pdf'), 
             width = 14, height = 10)
      
      
      if(!is.null(dev.list())) dev.off()
      pdf(paste0(outDir_cc, '/sparse_features_dubstepR.tfs.sps_umap.updated.pdf'),
          width =16, height = 10, useDingbats = FALSE)
      plot_manyFeatures_seurat(seurat_obj = sub_obj, features = tfs.dub)
      
      dev.off()
      
    }
    
    
  }
  
  ##########################################
  # test SMD
  ##########################################
  Test_SMD = TRUE
  if(Test_SMD){
    
    DimPlot(aa, group.by = 'condition', repel = TRUE)
    
    exp.mat = aa@assays$RNA@data
    mm = match(rownames(exp.mat), unique(c(tfs)))
    exp.mat = exp.mat[which(!is.na(mm)), ]
    cat(nrow(exp.mat), 'tfs detected \n')
    means = apply(exp.mat, 1, mean)
    sds = apply(exp.mat, 1, sd)
    # plot(means, sds)
    exp.mat = exp.mat[means > 0.05 & sds > 0.05,  ]
    cat(nrow(exp.mat), 'tfs after filtering \n')
    exp.mat = t(apply(exp.mat, 1, scale, center = TRUE, scale = TRUE))
    exp.mat = as.matrix(t(exp.mat))
    
    write.csv(exp.mat, file = paste0(outDir_cc, '/exp_matrix_TFs.SPs_4SMD.csv'), 
              quote = FALSE, row.names = FALSE)
    
    ##########################################
    # process the SMD output 
    ##########################################
    Process_SMD_ouptut = FALSE
    if(Process_SMD_ouptut){
      smd = read.csv(file = paste0(outDir_cc, 'SMD_output_tfs.sps.csv'), sep = '\t')
      #smd = read.csv(file = paste0(outDir, '/output_SMD_12k.cells_tfs.sps_v4.csv'), sep = '\t')
      smd = smd[, -1]
      smd = smd[order(-smd$SMD_z), ]
      
      saveRDS(smd, file = paste0(outDir_cc, 'output_SMD_tfs_sps.rds'))
      #smd = readRDS(file = paste0(outDir, '/output_SMD_12k.cells_tfs.sps_v4.rds'))
      
      ## test umap with SMD features
      genes.smd = smd$gene[which(smd$SMD_z>1)]
      cat(length(genes.smd), 'feature selected by smd \n')
      print(genes.smd)
      
      
      if(!is.null(dev.list())) dev.off()
      source('functions_utility.R')
      pdf(paste0(outDir_cc, '/sparse_features_SMD.tfs.sps_originalUMAP.pdf'),
          width =16, height = 10, useDingbats = FALSE)
      plot_manyFeatures_seurat(seurat_obj = aa, features = genes.smd)
      
      dev.off()
      
      #FeaturePlot(aa, features = smd$gene[which(smd$SMD_z>1)])
      #ggsave(filename = paste0(outDir, '/sparse_features_SMD_TFs_SPs.pdf'), 
      #       width = 24, height = 18) 
      
      
      Test_updatedUMAP_sparseFeatures = FALSE
      if(Test_updatedUMAP_sparseFeatures){
        sub_obj = subset(aa, features = genes.smd)
        sub_obj <- RunPCA(sub_obj, features = rownames(sub_obj), verbose = FALSE, weight.by.var = FALSE, 
                          npcs = 30)
        ElbowPlot(sub_obj, ndims = 30)
        
        sub_obj <- RunUMAP(sub_obj, dims = 1:30, n.neighbors = 50, min.dist = 0.1)
        DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
        
        p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE) +
          NoLegend()
        p1 =  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)+
          NoLegend()
        
        p0 + p1
        
        ggsave(filename = paste0(outDir, '/Umap_comparison_sparseFeatures_smd_45.tfs.sps.pdf'), 
               width = 12, height = 6) 
        
        FeaturePlot(sub_obj, features = smd$gene[which(smd$SMD_z>0.6)], max.cutoff = 'q1')
        
        ggsave(filename = paste0(outDir, '/sparse_features_umap.smd.features_SMD_31TFs_SPs.pdf'), 
               width = 24, height = 18)
        
      }
      
    }
   
  }
  
  ########################################################
  ########################################################
  # Section : Compare different methods: geneBasis, DUBStepR and SMD  
  # merge the gene list from DUBStepR, geneBasis and SMD
  ########################################################
  ########################################################
  merge_sparse_featureSelection = TRUE
  if(merge_sparse_featureSelection){
    
    ntop_sparseFeatures = 50
    
    dubstepR.out = readRDS(file = paste0(outDir_cc, 'sparse_features_dubstepR_test.v2.rds'))
    genes.dub = dubstepR.out$optimal.feature.genes
    length(genes.dub)
    length(intersect(genes.dub, unique(c(tfs, sps))))
    genes.dub = intersect(genes.dub, unique(c(tfs, sps)))
    
    if(length(genes.dub)>ntop_sparseFeatures) genes.dub = genes.dub[1:ntop_sparseFeatures]
    
    genes.basis = readRDS(paste0(outDir_cc,  
                                 'sparse_features_test_geneBasis_TFs.SPs_ntop', ntop_sparseFeatures, '.rds')) 
    genes.basis = genes.basis$gene
    
    cat(length(genes.dub), 'features from DUBStep \n')
    cat(length(genes.basis), 'features from geneBasis \n')
    
    ggs.merged = unique(c('Pax6', 'Foxa2', genes.dub, genes.basis))
    cat(length(ggs.merged), ' merged features from geneBasis and DUBStep \n')
    
    ggs.addition = c('Zfp42', 'Tcf15', 'Skil', 'Lef1', 
    'Tfap2c', 'Cebpb', 'Sox17', 'Prdm1', 'Cdh1', 'Irx3', 'Sox1', 'Dhrs3', 'Lhx1')
    
    setdiff(ggs.addition, ggs.merged)
    
    ggs.merged = unique(c(ggs.merged, ggs.addition))
    
    cat(length(ggs.merged), ' features after merging \n')
    
    genes.smd = readRDS(file = paste0(outDir_cc, 'output_SMD_tfs_sps.rds'))
    genes.smd = genes.smd$gene[which(genes.smd$SMD_z>1)]
    if(length(genes.smd) > ntop_sparseFeatures) genes.smd = genes.smd[1:ntop_sparseFeatures]
    ggs.merged = unique(c(ggs.merged, genes.smd))
    
    cat(length(ggs.merged), ' features after merging \n')
    
    saveRDS(ggs.merged, file = paste0(outDir_cc, 'sparseFeatures_merged_dubstepR.geneBasis.smd.rds'))
    
    Add_sparseFeatures_smd = TRUE
    if(!Add_sparseFeatures_smd){
      cat(length(ggs.merged), ' features after merging \n')
      saveRDS(ggs.merged, file = paste0(outDir_cc, 'sparseFeatures_merged_dubstepR.geneBasis.rds'))
      
      ggs.merged = readRDS(file = paste0(outDir_cc, 'sparseFeatures_merged_dubstepR.geneBasis.rds'))
      
    }
    
    
    Test_updatedUMAP_sparseFeatures = TRUE
    if(Test_updatedUMAP_sparseFeatures){
      # test seurat umap using selected sparse features
      ggs.merged = readRDS(file = paste0(outDir, '/day2.5_RA_day3_RA.rep1_day3.5_RA_day4_RA_day5_RA_rmCluster9.10/',
                                         'sparseFeatures_merged_dubstepR.geneBasis.smd.rds'))
      ggs.merged2 = readRDS(file = paste0(outDir, '/day2_beforeRA/', 
                                           'sparse_features_test_geneBasis_TFs.SPs_ntop30.rds'))
      
      ggs.merged = c(ggs.merged, ggs.merged2$gene)
      
      ggs.merged = readRDS(file = paste0(outDir_cc, 'sparseFeatures_merged_dubstepR.geneBasis.rds'))
      
      ggs.merged = readRDS(file = paste0(outDir_cc, 'sparseFeatures_merged_dubstepR.geneBasis.rds'))
      #ggs.merged = setdiff(ggs.merged, 'Foxq1')
      
      sub_obj <- RunPCA(aa, features = ggs.merged, verbose = FALSE, weight.by.var = FALSE, 
                        npcs = 30)
      ElbowPlot(sub_obj, ndims = 30)
      
      sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 20, min.dist = 0.05)
      DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
      
      sub_obj <- FindNeighbors(sub_obj, dims = 1:20)
      sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.5)
      
      p1 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
      p11 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
      p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
      p3 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
      (p1 + p11) / (p2 + p3) 
      
      ggsave(filename = paste0(outDir_cc, '/UMAP_RA_condition_clustering_cellcyclePhase_mergedSparseFeatures.pdf'), 
             width = 20, height = 12)
      
      
      p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
      p3 = FeaturePlot(sub_obj, features = c('Pax6', 'Foxa2'))
      
      p2 / p3
      
      ggsave(filename = paste0(outDir_cc, '/UMAP.using.merged.sparse.features_Pax6_Foxa2.pdf'), 
             width = 20, height = 12)
      
      if(!is.null(dev.list())) dev.off()
      pdf(paste0(outDir_cc, '/sparse_features_merged.geneBasis.dubstepR.tfs.sps_UMAP.updated.pdf'),
          width =16, height = 10, useDingbats = FALSE)
      plot_manyFeatures_seurat(seurat_obj = sub_obj, features = ggs.merged)
      
      dev.off()
      
      Check_someClusters = FALSE
      if(Check_someClusters){
        Idents(sub_obj) = sub_obj$seurat_clusters
        
        for(cluster in c(6)){
          
          # cluster = 4
          cluster.markers <- FindMarkers(sub_obj, only.pos = TRUE, ident.1 = cluster, ident.2 = 2, 
                                         min.pct = 0.2, logfc.threshold = 0.2, )
          head(cluster.markers, n = 10)
          
          FeaturePlot(sub_obj, features = c(rownames(cluster.markers)[1:10], 'Foxa2', 'Pax6'))
          
          ggsave(filename = paste0(outDir_cc, 
                                   '/RAtargets_chipseq_sparseMerged_FeaturePlots_markers_cluster.',
                                   cluster, '.pdf'), 
                 width = 14, height = 8)
          
        }
        
      }
      
    }
    
    
    ### save results 
    Save_table_plots = FALSE
    if(Save_table_plots){
      
      saveVersion = paste0(cc, collapse = '_')
      saveDir = paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/scRNAseq_mNT/',
                       'RA_symetryBreaking/sparse_featureSelection_timePoints/', saveVersion,
                       '/')
      system(paste0('mkdir -p ', saveDir))
      
      source('functions_utility.R')
      ggs.merged = merge_sparseFeatures(list(genes.dub, genes.basis, genes.smd, ggs.addition), 
                                        manual_addtion = TRUE)
      
      # test seurat umap using selected sparse features
      ggs = rownames(ggs.merged)
      # ggs = unique(ggs, res$genes)
      sub_obj = subset(aa, features = ggs)
      sub_obj <- RunPCA(sub_obj, features = rownames(sub_obj), verbose = FALSE, weight.by.var = FALSE, 
                        npcs = 50)
      ElbowPlot(sub_obj, ndims = 50)
      
      sub_obj <- RunUMAP(sub_obj, dims = 1:10, n.neighbors = 20, min.dist = 0.05)
      DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
      
      sub_obj <- FindNeighbors(sub_obj, dims = 1:10)
      sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.7)
      
      p1 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
      p11 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
      p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
      p3 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
      (p1 + p11) / (p2 + p3) 
      
      ggsave(filename = paste0(saveDir, '/updatedUMAP_condition_clustering_cellcyclePhase_',
                               saveVersion, '.pdf'), 
             width = 20, height = 12)
      
      write.csv2(ggs.merged, file = paste0(saveDir, 'merged_sparse_features_withRanks_', saveVersion, '.csv'), 
                 row.names = TRUE, quote = FALSE)
      
      if(!is.null(dev.list())) dev.off()
      pdf(paste0(saveDir, '/updatedUMAP_merged_sparseFeatures_', saveVersion, '.pdf'),
          width =16, height = 10, useDingbats = FALSE)
      plot_manyFeatures_seurat(seurat_obj = sub_obj, features = ggs)
      
      dev.off()
      
      p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
      p11 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
      p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
      p3 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
      (p1 + p11) / (p2 + p3) 
      
      ggsave(filename = paste0(saveDir, '/originalUMAP_condition_clustering_cellcyclePhase_', saveVersion, '.pdf'), 
             width = 20, height = 12)
      
      if(!is.null(dev.list())) dev.off()
      pdf(paste0(saveDir, '/originalUMAP_merged_sparseFeatures_', saveVersion, '.pdf'),
          width =16, height = 10, useDingbats = FALSE)
      plot_manyFeatures_seurat(seurat_obj = aa, features = ggs)
      
      dev.off()
      
      
    }
    
  }
  
  ########################################################
  ########################################################
  # Section : test features of RA signaling pathways and RA targets
  # 
  ########################################################
  ########################################################
  Check_RA.genes_RAtargets = FALSE
  if(Check_RA.genes_RAtargets){
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
    
    Only_RA_genes = TRUE
    if(!Only_RA_genes){
      #genes.baisis = readRDS(file = paste0(outDir_cc, 
      #                                    'sparse_features_test_geneBasis_TFs.SPs_ntop', ntop_sparseFeatures, '.rds')) 
      #ggs = unique(c(ggs, genes.baisis$gene))
      
      genes.alltime.merged = readRDS(file = paste0(outDir_cc, 'sparseFeatures_merged_dubstepR.geneBasis.SMD_v2.rds'))
      ggs = unique(c(ggs, genes.alltime.merged))
    }
    
    ggs = intersect(ggs, rownames(aa))
    cat(length(ggs), ' RA signaling bias features \n')
    
    sub_obj <- RunPCA(aa, features = ggs, verbose = FALSE, weight.by.var = FALSE, 
                      npcs = 30)
    ElbowPlot(sub_obj, ndims = 30)
    
    sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
    DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    
    #p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    p1 =  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
    p1 + p2
    
    ggsave(filename = 
             paste0(outDir_cc, '/Umap_comparison_origine_sparseFeatures_RA.pathways.pdf'), 
           width = 10, height = 6) 
    
    FeaturePlot(sub_obj, features = c('Pax6', 'Foxa2', 'Zfp42', 'Tcf15', 'Skil', 'Lef1', 
                                      'Tfap2c', 'Cebpb', 'Sox17', 'Prdm1', 'Cdh1', 'Irx3', 'Sox1',
                                      "Cyp26a1", "Cyp26b1"))
    
    ggsave(filename = paste0(outDir,
                             '/check_representiveFeature_umapEmbedding_mergedFeatures_RA.pathways.pdf'), 
           width = 14, height = 10)
    
    sub_obj <- FindNeighbors(sub_obj, dims = 1:20)
    sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.7)
    
    p1 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
    p1 + p2
    
    Find_cluster_specific_markers = FALSE
    if(Find_cluster_specific_markers){
      Idents(sub_obj) = sub_obj$seurat_clusters
      
      cluster4.markers <- FindMarkers(sub_obj, only.pos = TRUE, ident.1 = 4, ident.2 = 5, 
                                      min.pct = 0.25, logfc.threshold = 0.1)
      head(cluster4.markers, n = 10)
      
      FeaturePlot(sub_obj, features = c(rownames(cluster4.markers)[1:10],'Foxa2', 'Pax6'))
      ggsave(filename = paste0(outDir, '/FeaturePlots_markers_cluster4.pdf'), 
             width = 14, height = 8)
      
      cluster5.markers <- FindMarkers(sub_obj, only.pos = TRUE, ident.1 = 5, ident.2 = 4, 
                                      min.pct = 0.25, logfc.threshold = 0.1)
      head(cluster5.markers, n = 10)
      
      FeaturePlot(sub_obj, features = c(rownames(cluster5.markers)[1:10], 'Foxa2', 'Pax6'))
      ggsave(filename = paste0(outDir, '/FeaturePlots_markers_cluster5.pdf'), 
             width = 14, height = 8)
      
      
    }
    
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
    ggs = read.csv2(paste0(RARtargetDir, 'RARtarget_intersection/DESeq2_DEgenes_pairwiseComparison_RA_d2.18h.vs.noRA_d2.18h.csv'),
                    row.names = c(1))
    ggs = rownames(ggs)[which(abs(ggs$log2FoldChange_RA_d2.18hvs.noRA_d2.18h)>2 & 
                                ggs$padj_RA_d2.18hvs.noRA_d2.18h<0.01)]
    
    #saveRDS(ggs, file = paste0(RdataDir, '/RA_targets_smartseq2_chipseq_strigent.rds'))
    genes.merged = readRDS(file = paste0(RdataDir, 'sparseFeatures_merged_dubstepR.geneBasis.SMD_v1.rds'))
    ggs = unique(c(ggs, genes.merged))
    
    ggs = intersect(ggs, rownames(aa))
    cat(length(ggs), ' genes selected \n')
    sub_obj <- RunPCA(aa, features = ggs, verbose = FALSE, weight.by.var = FALSE, 
                      npcs = 30)
    ElbowPlot(sub_obj, ndims = 30)
    
    sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
    DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    
    p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    p1 =  DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    
    p0 + p1
    
    ggsave(filename = 
             paste0(outDir, '/Umap_comparison_origine_selectedFeatures_RAtargets_chipseq_sparseMerged.pdf'), 
           width = 16, height = 6) 
    
    sub_obj <- FindNeighbors(sub_obj, dims = 1:20)
    sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.5)
    p1 = DimPlot(sub_obj, label = TRUE, repel = TRUE,  raster=FALSE)
    p2 = DimPlot(sub_obj, group.by = 'Phase')
    p1 + p2
  }
  
  
}

##########################################
# save the workspace  
##########################################
save.image(paste0(outDir_cc, "saved_Image.RData"))
