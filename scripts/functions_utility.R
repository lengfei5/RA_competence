##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: utility functions 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Oct 25 10:37:05 2022
##########################################################################
##########################################################################
cal_sample_means = function(cpm, conds = c("mUA", "mLA", "mHand"))
{
  sample.means = c()
  for(n in 1:length(conds)) 
  {
    kk = grep(conds[n], colnames(cpm))
    if(length(kk)>1) {
      sample.means = cbind(sample.means, apply(cpm[, kk], 1, mean))
    }else{
      sample.means = cbind(sample.means, cpm[, kk])
    }
    
  }
  colnames(sample.means) = conds
  return(sample.means)
  
}


##########################################
# plot many features in one pdf 
##########################################
plot_manyFeatures_seurat = function(seurat_obj, features)
{
  #if(!is.null(dev.list())) dev.off()
  #while(!is.null(dev.list())){dev.off()}
  
  nn = ceiling(length(features)/10)
  for(n in 1:nn)
  {
    kk = c(((n-1)*10+1):(n*10))
    cat(' -- ', kk, '-- \n')
    p = FeaturePlot(seurat_obj, features = features[kk])
    plot(p)                                         
  }
}


merge_sparseFeatures = function(gene.list = list(genes.dub, genes.basis, genes.smd, ggs.addtions), 
                                manual_addtion = TRUE)
{
  gg.uniq = c()
  for(n in 1:length(gene.list)){
    gg.uniq = unique(c(gg.uniq, gene.list[[n]]))
  }
  
  gg.merged = matrix(NA, ncol = length(gene.list), nrow = length(gg.uniq))
  rownames(gg.merged) = gg.uniq
  colnames(gg.merged) = c('dubstep', 'gene.basis', 'smd', 'manual.added')
  
  for(n in 1:length(gene.list)){
    gg.merged[,n] = match(rownames(gg.merged), gene.list[[n]])
  }
  
  if(manual_addtion){
    gg.merged[which(!is.na(gg.merged[, ncol(gg.merged)])), ncol(gg.merged)] = 1  
    ranks = gg.merged[, c(1:(ncol(gg.merged)-1))]
  }else{ranks = gg.merged}
  
  ranks[which(is.na(ranks))] = nrow(ranks)
  ss = apply(ranks, 1, sum)
  ss_order = order(ss)
  
  gg.merged = gg.merged[ss_order, ]
  
  
  return(gg.merged) 
  
}

##########################################
# 
##########################################
## original functions from Monocle
## https://github.com/cole-trapnell-lab/monocle-release/blob/master/R/plotting.R
plot_genes_branched_heatmap <- function(seuratObj = aa,
                                        gene_subset = candidates,
                                        nbCell_condition = 50,
                                        Get.Smooth.Curve = TRUE, 
                                        scale_max=3, 
                                        scale_min=-3, 
                                        hmcols = NULL, 
                                        hclust_method = "ward.D2", 
                                        num_clusters = 6,
                                        cluster_rows = TRUE,
                                        add_annotation_row = NULL,
                                        show_rownames = TRUE
                                        ) 
{
  # seuratObj = aa; nbCell_condition = 100;scale_max=3; scale_min=-3;hmcols = NULL; Get.Smooth.Curve = TRUE;
  # gene_subset = candidates;hclust_method = "ward.D2";num_clusters = 6
  library(VGAM) # an example code from https://online.stat.psu.edu/stat504/lesson/8/8.2/8.2.2
  library(MASS)
  library(tidyr)
  require(colorRamps)
  require(pheatmap)
  
  cell_beforeRA = colnames(seuratObj)[which(!is.na(seuratObj$pseudot) & 
                                              seuratObj$condition == 'day2_beforeRA')]
  cell_noRA = colnames(seuratObj)[which(!is.na(seuratObj$pseudot) & grepl('_noRA', seuratObj$condition))]
  cell_RA = colnames(seuratObj)[which(!is.na(seuratObj$pseudot) & grepl('_RA', seuratObj$condition))]
  
  cat('subsampling ', nbCell_condition, ' cells\n')
  cell.sels = c()
  cc = unique(seuratObj$condition)
  cc = cc[which(cc != "day3_RA.rep2")]
  
  for(n in 1:length(cc))
  {
    cell.sels = c(cell.sels, sample(colnames(seuratObj)[which(seuratObj$condition == cc[n] & 
                                                                !is.na(seuratObj$pseudot))], 
                                    size = nbCell_condition, 
                                    replace = FALSE))
  }
  
  subs = subset(seuratObj, cells = cell.sels)
  
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
  
  if(Get.Smooth.Curve){
    cds <- seuratObj@assays$RNA@scale.data
    cds = cds[which(!is.na(match(rownames(cds), gene_subset))), ]
    cat(' -- smoothing the single cell data for subsampled cells -- \n')  
    
    # before RA
    jj = match(cell_beforeRA, colnames(seuratObj))
    jj = jj[which(!is.na(seuratObj$pseudot[jj]))]
    Pseudotime = as.numeric(seuratObj$pseudot[jj])
    kk_common = which(subs$condition == 'day2_beforeRA')
    kk_common = kk_common[order(subs$pseudot[kk_common])]
    pseudot_comomon = subs$pseudot[kk_common]
    
    common_ancestor_cells = t(apply(cds[ ,jj], 1, get_smooth_curve_spline, t = Pseudotime, newt = pseudot_comomon))
    
    jj = grep('_RA$|_RA.rep1', seuratObj$condition)
    jj = jj[which(!is.na(seuratObj$pseudot[jj]))]
    Pseudotime = as.numeric(seuratObj$pseudot[jj])
    
    kk_BrachA = grep('_RA$|_RA.rep1', subs$condition)
    kk_BrachA = kk_BrachA[order(subs$pseudot[kk_BrachA])]
    pseudot_BrachA = subs$pseudot[kk_BrachA]
    
    BranchA_exprs <- t(apply(cds[ ,jj], 1, get_smooth_curve_spline, t = Pseudotime, newt = pseudot_BrachA))
    
    jj = grep('_noRA', seuratObj$condition)
    jj = jj[which(!is.na(seuratObj$pseudot[jj]))]
    Pseudotime = as.numeric(seuratObj$pseudot[jj])
    
    kk_BrachB = grep('_noRA', subs$condition)
    kk_BrachB = kk_BrachB[order(subs$pseudot[kk_BrachB])]
    pseudot_BrachB = subs$pseudot[kk_BrachB]
    
    BranchB_exprs <- t(apply(cds[ ,jj], 1, get_smooth_curve_spline, t = Pseudotime, newt = pseudot_BrachB))
    
  }
  
  col_gap_ind <- c(length(kk_BrachB), length(kk_BrachB) + length(kk_common), 
                   length(kk_BrachB) + 2*length(kk_common))
  
  heatmap_matrix <- cbind(BranchB_exprs[, ncol(BranchB_exprs):1], 
                          common_ancestor_cells[, ncol(common_ancestor_cells):1],
                          common_ancestor_cells,
                          BranchA_exprs)
  
  indexs = c(kk_BrachB[ncol(BranchB_exprs):1], 
             kk_common[ncol(common_ancestor_cells):1], 
             kk_common,
             kk_BrachA)
  
  heatmap_matrix=heatmap_matrix[!apply(heatmap_matrix, 1, sd)==0,]
  heatmap_matrix=Matrix::t(scale(Matrix::t(heatmap_matrix), center=TRUE))
  heatmap_matrix=heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE, ]
  heatmap_matrix[is.nan(heatmap_matrix)] = 0
  heatmap_matrix[heatmap_matrix>scale_max] = scale_max
  heatmap_matrix[heatmap_matrix<scale_min] = scale_min
  
  saveRDS(heatmap_matrix, file = paste0(outDir, "/heatmap_matrix_forPlot.rds"))
  #heatmap_matrix_ori <- heatmap_matrix
  #heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & is.finite(heatmap_matrix[, col_gap_ind]), ] #remove the NA fitting failure genes for each branch 
  
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  
  exp_rng <- range(heatmap_matrix) #bks is based on the expression range
  
  bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
  if(is.null(hmcols)) {
    hmcols <- blue2green2red(length(bks) - 1)
  }
  
  # prin  t(hmcols)
  ph <- pheatmap(heatmap_matrix, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=TRUE, 
                 show_rownames=F, 
                 show_colnames=F, 
                 #scale="row",
                 clustering_distance_rows=row_dist, #row_dist
                 clustering_method = hclust_method,
                 cutree_rows=num_clusters,
                 silent=TRUE,
                 #filename=NA,
                 breaks=bks,
                 color=hmcols
                 #color=hmcols#
  )
  
  annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)),
                               pseudot = subs$pseudot[indexs],
                               condition = subs$condition[indexs])
  
  write.csv2(annotation_row, file = paste0(outDir, 'gene_clusters.csv'), row.names = TRUE, quote = FALSE)
  
  #colnames(annotation_col) <- "Cell Type"  
  #if(!is.null(add_annotation_col)) {
  #  annotation_col <- cbind(annotation_col, add_annotation_col[fData(cds[row.names(annotation_col), ])$gene_short_name, 1])  
  #}
  branch_colors = cols_sel[grep("day3_RA.rep2", names(cols_sel), invert = TRUE)]
  annotation_colors=list("condition"=branch_colors)
  #names(annotation_colors$`condition`) = c('Pre-branch', branch_labels)
  
  
  #names(branch_colors) <- c("Pre-branch", branch_labels[1], branch_labels[2])
  #annotation_colors=list("Cell Type"=branch_colors)
  #names(annotation_colors$`Cell Type`) = c('Pre-branch', branch_labels)
  feature_label <- row.names(heatmap_matrix)
  row_ann_labels <- row.names(annotation_row)
  
  row.names(heatmap_matrix) <- feature_label
  row.names(annotation_row) <- row_ann_labels
  
  ##########################################
  # all DE genes
  ##########################################
  pheatmap(heatmap_matrix[, ], #ph$tree_row$order
           useRaster = T,
           cluster_cols=FALSE, 
           cluster_rows=TRUE, 
           show_rownames=FALSE,
           show_colnames=FALSE, 
           scale='none',
           clustering_distance_rows=row_dist, #row_dist
           clustering_method = hclust_method, #ward.D2
           cutree_rows=num_clusters,
           # cutree_cols = 2,
           annotation_row=annotation_row,
           annotation_col=annotation_col,
           annotation_colors=annotation_colors,
           gaps_col = col_gap_ind,
           treeheight_row = 30, 
           breaks=bks,
           fontsize = 6,
           color=hmcols, 
           border_color = NA,
           silent=TRUE, 
           filename=paste0(outDir, "/expression_pseudotime_pheatmap_allDEgenes.pdf"),
           width = 6, height = 12
  )
  
  pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                     useRaster = T,
                     cluster_cols=FALSE, 
                     cluster_rows=TRUE, 
                     show_rownames=TRUE,
                     show_colnames=FALSE, 
                     scale='none',
                     clustering_distance_rows=row_dist, #row_dist
                     clustering_method = hclust_method, #ward.D2
                     cutree_rows=num_clusters,
                     # cutree_cols = 2,
                     annotation_row=annotation_row,
                     annotation_col=annotation_col,
                     annotation_colors=annotation_colors,
                     gaps_col = col_gap_ind,
                     treeheight_row = 50, 
                     breaks=bks,
                     fontsize = 2,
                     color=hmcols, 
                     border_color = NA,
                     silent=TRUE, 
                     filename=paste0(outDir, "/expression_pseudotime_pheatmap_allDEgenes_with.geneNames.pdf"),
                     width = 12, height = 60
  )
  
  
  ##########################################
  # all DE genes intersected with RAR target
  ##########################################
  targets = readRDS('../results/RA_targets_L118404_smartseq3_20221117/Rdata/RAR_targets_chip_chiapet.rds')
  
  sels = which(!is.na(match(rownames(heatmap_matrix), targets)))
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix[sels,])))/2)
  row_dist[is.na(row_dist)] <- 1
  annotation_rowSel = data.frame(Cluster = annotation_row[sels, ])
  rownames(annotation_rowSel) = rownames(annotation_row)[sels]
  
  
  pheatmap(heatmap_matrix[sels, ], #ph$tree_row$order
           useRaster = T,
           cluster_cols=FALSE, 
           cluster_rows=TRUE, 
           show_rownames=FALSE,
           show_colnames=FALSE, 
           scale='none',
           clustering_distance_rows=row_dist, #row_dist
           clustering_method = hclust_method, #ward.D2
           cutree_rows=num_clusters,
           # cutree_cols = 2,
           annotation_row=annotation_rowSel,
           annotation_col=annotation_col,
           annotation_colors=annotation_colors,
           gaps_col = col_gap_ind,
           treeheight_row = 30, 
           breaks=bks,
           fontsize = 6,
           color=hmcols, 
           border_color = NA,
           silent=TRUE, 
           filename=paste0(outDir, "/expression_pseudotime_pheatmap_allDEgenes_intersectedRARtarget.pdf"),
           width = 6, height = 12
  )
  
  pheatmap(heatmap_matrix[sels, ], #ph$tree_row$order
           useRaster = T,
           cluster_cols=FALSE, 
           cluster_rows=TRUE, 
           show_rownames=TRUE,
           show_colnames=FALSE, 
           scale='none',
           clustering_distance_rows=row_dist, #row_dist
           clustering_method = hclust_method, #ward.D2
           cutree_rows=num_clusters,
           # cutree_cols = 2,
           annotation_row=annotation_rowSel,
           annotation_col=annotation_col,
           annotation_colors=annotation_colors,
           gaps_col = col_gap_ind,
           treeheight_row = 30, 
           breaks=bks,
           fontsize = 2,
           color=hmcols, 
           border_color = NA,
           silent=TRUE, 
           filename=paste0(outDir, "/expression_pseudotime_pheatmap_allDEgenes_intersectedRARtarget",
                           "_withgeneNames.pdf"),
           width = 8, height = 20
  )
  
  
  
  ##########################################
  # DE TFs and signaling pathways 
  ##########################################
  tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
  tfs = unique(tfs$`HGNC symbol`)
  tfs = as.character(unlist(sapply(tfs, firstup)))
  sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v2.rds'))
  targets = unique(c(tfs, sps$gene))
  xx = read.table('../data/annotations/GO_term_summary_RAR_signalingPathway.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  targets = unique(c(targets, xx[,2]))
  xx = read.table('../data/annotations/GO_term_summary_TGFb.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  targets = unique(c(targets, xx[,2]))
  #sps = toupper(unique(sps$gene))
  #sps = setdiff(sps, tfs)
  
  sels = which(!is.na(match(rownames(heatmap_matrix), targets)))
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix[sels,])))/2)
  row_dist[is.na(row_dist)] <- 1
  annotation_rowSel = data.frame(Cluster = annotation_row[sels, ])
  rownames(annotation_rowSel) = rownames(annotation_row)[sels]
  
  pheatmap(heatmap_matrix[sels, ], #ph$tree_row$order
           useRaster = T,
           cluster_cols=FALSE, 
           cluster_rows=TRUE, 
           show_rownames=FALSE,
           show_colnames=FALSE, 
           scale='none',
           clustering_distance_rows=row_dist, #row_dist
           clustering_method = hclust_method, #ward.D2
           cutree_rows=num_clusters,
           # cutree_cols = 2,
           annotation_row=annotation_rowSel,
           annotation_col=annotation_col,
           annotation_colors=annotation_colors,
           gaps_col = col_gap_ind,
           treeheight_row = 30, 
           breaks=bks,
           fontsize = 6,
           color=hmcols, 
           border_color = NA,
           silent=TRUE, 
           filename=paste0(outDir, "/expression_pseudotime_pheatmap_allDEgenes_TF.SP.pdf"),
           width = 6, height = 12
  )
  
  pheatmap(heatmap_matrix[sels, ], #ph$tree_row$order
           useRaster = T,
           cluster_cols=FALSE, 
           cluster_rows=TRUE, 
           show_rownames=TRUE,
           show_colnames=FALSE, 
           scale='none',
           clustering_distance_rows=row_dist, #row_dist
           clustering_method = hclust_method, #ward.D2
           cutree_rows=num_clusters,
           # cutree_cols = 2,
           annotation_row=annotation_rowSel,
           annotation_col=annotation_col,
           annotation_colors=annotation_colors,
           gaps_col = col_gap_ind,
           treeheight_row = 30, 
           breaks=bks,
           fontsize_row = 4,
           #fontsize = 4,
           color=hmcols, 
           border_color = NA,
           silent=TRUE, 
           filename=paste0(outDir, "/expression_pseudotime_pheatmap_allDEgenes_TF.SP",
                           "_withgeneNames.pdf"),
           width = 8, height = 20
  )
  
  
  
  ##########################################
  # DE TFs and signaling pathways intersected with RAR target
  ##########################################
  targets = readRDS('../results/RA_targets_L118404_smartseq3_20221117/Rdata/RAR_targets_chip_chiapet.rds')
  
  tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
  tfs = unique(tfs$`HGNC symbol`)
  tfs = as.character(unlist(sapply(tfs, firstup)))
  sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v2.rds'))
  tf.sp = unique(c(tfs, sps$gene))
  xx = read.table('../data/annotations/GO_term_summary_RAR_signalingPathway.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  tf.sp = unique(c(tf.sp, xx[,2]))
  xx = read.table('../data/annotations/GO_term_summary_TGFb.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  tf.sp = unique(c(tf.sp, xx[,2]))
  
  targets = intersect(targets, tf.sp)
    
  sels = which(!is.na(match(rownames(heatmap_matrix), targets)))
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix[sels,])))/2)
  row_dist[is.na(row_dist)] <- 1
  annotation_rowSel = data.frame(Cluster = annotation_row[sels, ])
  rownames(annotation_rowSel) = rownames(annotation_row)[sels]
  
  pheatmap(heatmap_matrix[sels, ], #ph$tree_row$order
           useRaster = T,
           cluster_cols=FALSE, 
           cluster_rows=TRUE, 
           show_rownames=TRUE,
           show_colnames=FALSE, 
           scale='none',
           clustering_distance_rows=row_dist, #row_dist
           clustering_method = hclust_method, #ward.D2
           cutree_rows=num_clusters,
           # cutree_cols = 2,
           annotation_row=annotation_rowSel,
           annotation_col=annotation_col,
           annotation_colors=annotation_colors,
           gaps_col = col_gap_ind,
           treeheight_row = 30, 
           breaks=bks,
           fontsize = 6,
           color=hmcols, 
           border_color = NA,
           silent=TRUE, 
           filename=paste0(outDir, "/expression_pseudotime_pheatmap_all.TF.SP_intersectedRARtargets",
                           "_withgeneNames.pdf"),
           width = 6, height = 12
  )
  
}


do_scatter <- function(umap_use, meta_data, label_name, no_guides = TRUE,
                       do_labels = TRUE, nice_names, 
                       palette_use = colors_use,
                       pt_size = 4, point_size = .5, base_size = 12, 
                       do_points = TRUE, do_density = FALSE, h = 6, w = 8) {
  # umap_use = harmony_embeddings; label_name = 'dataset'
  library(ggplot2)
  #colors_use <- c(`jurkat` = '#810F7C', `t293` = '#D09E2D',`half` = '#006D2C')
  
  umap_use <- umap_use[, 1:2]
  colnames(umap_use) <- c('X1', 'X2')
  plt_df <- umap_use %>% data.frame() %>% 
    cbind(meta_data) %>% 
    dplyr::sample_frac(1L) 
  plt_df$given_name <- plt_df[[label_name]]
  
  if (!missing(nice_names)) {
    plt_df %<>%
      dplyr::inner_join(nice_names, by = "given_name") %>% 
      subset(nice_name != "" & !is.na(nice_name))
    
    plt_df[[label_name]] <- plt_df$nice_name        
  }
  
  plt <- plt_df %>% 
    ggplot2::ggplot(aes_string("X1", "X2", col = label_name, fill = label_name)) + 
    theme_test(base_size = base_size) + 
    theme(panel.background = element_rect(fill = NA, color = "black")) + 
    guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1,
                                                    shape = 16, size = 4)), 
           alpha = FALSE) +
    scale_color_manual(values = palette_use) + 
    scale_fill_manual(values = palette_use) +    
    theme(plot.title = element_text(hjust = .5)) + 
    labs(x = "PC 1", y = "PC 2") 
  
  if (do_points) 
    plt <- plt + geom_point(shape = '.')
  if (do_density) 
    plt <- plt + geom_density_2d()    
  
  
  if (no_guides)
    plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)
  
  if (do_labels) {
    data_labels <- plt_df %>% 
      dplyr::group_by(label_name) %>% 
      dplyr::summarise(X1 = mean(X1), X2 = mean(X2)) %>% 
      dplyr::ungroup()
    
    plt <- plt + geom_label(data = data_labels, label.size = NA,
                            aes_string(label = label_name), 
                            color = "white", size = pt_size, alpha = 1,
                            segment.size = 0) +
      guides(col = FALSE, fill = FALSE)
  }
  
  return(plt)
}


########################################################
########################################################
# Section IV : Genie3 random-forest-based method to predict TF regulators for FoxA2
# 
########################################################
########################################################
Test_GENIE3 = FALSE
if(Test_GENIE3){
  motif.oc = readRDS(file = paste0(outDir, '/motif_oc_fimo_jaspar2022_pval.0.0001_Foxa2.enhancers.promoters.rds'))
  motif.oc = t(motif.oc)
  mapping = readRDS(paste0('../results/scRNAseq_R13547_10x_mNT_20220813/motif_analysis/', 
                           'JASPAR2022_CORE_UNVALIDED_vertebrates_nonRedundant_metadata_manual_rmRedundantUNVALIDED.rds'))
  
  mapping = mapping[which(!is.na(match(mapping$name, rownames(motif.oc)))), ]
  mapping$gene = firstup(mapping$gene)
  
  Idents(aa) = factor(aa$condition)
  
  levels_sels = c(
    "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA")
  
  data_version = "d2.5.d3,d3.5,d4"
  
  outDir_cc = paste0(outDir, '/Genie3_FoxA2', data_version)
  system(paste0('mkdir -p ', outDir_cc))
  
  sce = subset(aa, idents = levels_sels)
  sce$condition = droplevels(sce$condition)
  #sce = as.SingleCellExperiment(sce)
  
  E = sce@assays$RNA@data
  gg_sels = intersect(c(mapping$gene, 'Foxa2', 'Lhx1'), rownames(E))
  
  E = E[match(gg_sels, rownames(E)), ]
  ss = rowSums(E)
  E = E[ss>0, ]
  
  E = as.matrix(E)
  
  saveRDS(E, file = paste0('analysis_examples/Genie3/ExprMatrix_4FoxA2.rds'))
  
  ##########################################
  # run the GENIE3 on the laptop due to dynamical.loading issue
  ##########################################
  # tic()
  # source('myGENIE3.R')
  # wtm = GENIE3_withSpecificTargets(expr.matrix = E, priorTargets = c('Foxa2'), ncore = 1)
  # 
  # saveRDS(wtm, file = paste0(RdataDir, '/first_test_Genie3_v3.rds'))
  # 
  # toc()
  
  
  wtm = readRDS(file = paste0('analysis_examples/Genie3/ranked_predictedRegulators_FoxA2.rds'))
  
  
  pdf(paste0(outDir_cc, '/GENIE3_regulators_expressionPattern.pdf'),
      width =16, height = 10, useDingbats = FALSE)
  
  features = c('Foxa2', 'Pax6', wtm$regulator[which(wtm$Foxa2>10^-3)])
  
  plot_manyFeatures_seurat(seurat_obj = aa, features = unique(features))
  
  dev.off()
  
  FeaturePlot(aa, features = c('Pax6', 'Foxa2', 'Rarg', 'Pou5f1', 'Zfp42', 'Lhx1', 'Jun', 'Lef1', 'Gata3'))
  
}

########################################################
########################################################
# Section : functions of quantifying cell heterogeneity
# 
########################################################
########################################################
calc_heterogeneity_RA.noRA = function(seuratObj, method = 'pairwiseDist_pca', 
                                      subsample.cells = 200,
                                      subsample.sketch = FALSE,
                                      nb_features = 500, nb_pcs = 30)
{
  # seuratObj = aa;method = 'pairwiseDist_pca';subsample.cells = 200;nb_features = 1000;nb_pcs = 30;
  
  seuratObj = subset(seuratObj, cells = colnames(seuratObj)[which(!is.na(seuratObj$pst_group))])
  
  Idents(seuratObj) = as.factor(seuratObj$pst_group)
  
  cc = unique(seuratObj$pst_group)
  
  if(method == 'pairwiseDist_pca'){
    
    library(scater)
    library(SingleCellExperiment)
    library(scran)
    
    set.seed(2024)
    
    if(subsample.sketch){
      seuratObj = subset(seuratObj, cells = colnames(seuratObj)[which(!is.na(seuratObj$sketch))])
    }else{
      seuratObj = subset(seuratObj, downsample = subsample.cells)
    }
    #table(seuratObj$pst_group)
    
    #Idents(seuratObj) = factor(seuratObj$condition, levels = levels_sels)
    #cc = unique(seuratObj$condition)
    
    Run_HVGs_perTimePoint = FALSE
    if(Run_HVGs_perTimePoint){
      hete = c()
      #nb_features = 1000
      for(n in 1:length(levels_sels))
      {
        # n = 1
        cat(n, ' -- ', levels_sels[n], '\n')
        sce = as.SingleCellExperiment(subset(seuratObj, idents = levels_sels[n]))
        
        dec <- modelGeneVar(sce)
        top.hvgs <- getTopHVGs(dec, n=nb_features)
        sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 30)
        # reducedDimNames(sce)
        #ll.pca = reducedDim(sce, 'PCA')[, c(1:30)]
        ll.pca = logcounts(sce)
        ll.pca = as.matrix(ll.pca[match(top.hvgs, rownames(ll.pca)), ])
        dists = sqrt((1-cor((ll.pca), method = 'pearson'))/2)
        
        #dists = dists[which(dists>0)]
        
        #hete = rbind(hete, data.frame(condition = rep(levels_sels[n], length(dists)), dists))
        hete = rbind(hete, data.frame(condition = rep(levels_sels[n], length(dists)), dists))
        
      }
    }else{
      
      sce = as.SingleCellExperiment(seuratObj)
      
      dec <- modelGeneVar(sce)
      top.hvgs <- getTopHVGs(dec, n=nb_features)
      sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = nb_pcs)
      
      # reducedDimNames(sce)
      
      ll.pca = reducedDim(sce, 'PCA')[, c(1:nb_pcs)]
      
      #ll.pca = logcounts(sce)
      #ll.pca = t(as.matrix(ll.pca[match(top.hvgs, rownames(ll.pca)), ]))
      
      hete = c()
      
      for(n in 1:length(cc))
      {
        # n = 1
        kk = match(colnames(sce)[which(sce$pst_group == cc[n])], rownames(ll.pca))
        
        cat(n, ' -- ', cc[n], '--', length(kk), 'cells \n')
        
        dists =  dist(ll.pca[kk, ], method = "euclidean", diag = FALSE, upper = TRUE) 
        #sqrt((1-cor((ll
        
        #dists = sqrt((1-cor(ll.pca[kk, ], method = 'pearson'))/2)
        #dists = dists[which(dists>0)]
        
        hete = rbind(hete, data.frame(nb_cells = rep(length(kk), length(dists)),
                                      condition = rep(cc[n], length(dists)),
                                      as.numeric(dists)))
        
      }
      
    }
    
    #hete$dists = log10(hete$dists)
    colnames(hete)[ncol(hete)] = 'dists'
    hete = as.data.frame(hete)
    
  }
  
  if(method == 'entropy_metric'){
    suppressMessages(library(ROGUE))
    suppressMessages(library(ggplot2))
    suppressMessages(library(tidyverse))
    
    hete = c()
    
    for(n in 1:length(cc))
    {
      # n = 1
      #cat(n, ' -- ', cc[n], '\n')
      
      kk = which(seuratObj$pst_group == cc[n])
      expr = seuratObj@assays$RNA@counts[, kk]
      
      #meta <- seuratObj@meta.data
      ent.res <- SE_fun(expr)
      ent.res <- SE_fun(expr, span = 0.1, r = 1, mt.method = "fdr")
      #head(ent.res)
      
      SEplot(ent.res)
      
      rogue.value <- CalculateRogue(ent.res, platform = "UMI")
      cat(n, ' -- ', cc[n], '--', rogue.value, '\n') 
      
      #rogue.res <- rogue(expr, labels = meta$ct, samples = meta$Patient, platform = "UMI", span = 0.6)
      #rogue.res
      #ent.res <- SE_fun(expr, span = 0.1, r = 1, mt.method = "fdr")
      #CalculateRogue(ent.res, platform = "UMI")
      #CalculateRogue(ent.res, k = 30)
      
      hete = rbind(hete,  c(cc[n], (1.0 - rogue.value)))
      
    }
    
    colnames(hete) = c('condition', 'hete')
    
    hete = as.data.frame(hete)
    
    
  }
  
  
  
  return(hete)
  
}

##########################################
# plot the gene-gene scatterplot from the output scFates 
##########################################
make_scatterplot_scFates = function(geneA = 'Foxa2', geneB = 'Pax6', fit, assign, pst,
                                    win_keep=c(2, 3, 4))
{
  library(gridExtra)
  
  # geneA = 'Foxa2'; geneB = 'Pax6'; win_keep=c(2, 3, 4, 5)
  if(length(!is.na(match(geneA, colnames(fit)))) != 1){
    stop('-- geneA not found --')
  }
  
  if(length(!is.na(match(geneB, colnames(fit)))) != 1){
    stop('-- geneB not found --')
  }
  
  jj1 = which(colnames(fit) == geneA)
  jj2 = which(colnames(fit) == geneB)
  nt = nrow(assign)
  cat('there are ', nt, 'non-intersecting time windowns \n')
  assignment = apply(assign, 2, which.max) 
  
  
  cmd = c()
  for(n in win_keep)
  {
    # n = 3
    cells = names(which(assignment == n))
    cells = gsub('[.]', '-', cells)
    mm = match(cells, rownames(fit))
    mm = mm[which(!is.na(mm))]
    
    kk = match(rownames(fit)[mm], rownames(pst))
    
    dat = data.frame(fit[mm, c(jj1, jj2)], pst = pst$t[kk])
    colnames(dat) = c('geneA', 'geneB', 'pseudotime')
    title = paste0('pst-window : ', n)
    eval(parse(text= paste0('p', n,  '=  ggplot(data = dat, aes(x = geneA, y=geneB)) +
      geom_point(size = 0.5) + 
      theme_classic() +
      labs(x = geneA, y = geneB ) +
      theme(axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10)) + 
      ggtitle(title)')))
    cmd = c(cmd, paste0('p', n))
           
  }
  
  figure = eval(parse(text= paste0(cmd, collapse = '+')))
  
  plot(figure)
  
}


