##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: analyze the signaling pathways for semmetry breaking
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Jan 13 11:01:15 2023
##########################################################################
##########################################################################
# all 16 samples
aa = readRDS(file = paste0(RdataDir, 
              'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
              'cellCycleScoring_annot.v1_savedUMAP.v1_', species, version.analysis, '.rds'))

# only RA samples incl. dya2_beforeRA
aa = readRDS(file = paste0(RdataDir, 
              'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
              'cellCycleScoring_annot.v1_savedUMAP.subs.v2_', species, version.analysis, '.rds'))

Idents(aa) = aa$condition
DimPlot(aa, cols = cols, group.by = 'condition', label = TRUE, repel = TRUE, raster = FALSE)

levels_sels = c("day3_RA.rep1", "day3.5_RA", "day4_RA")
data_version = "_d3_d3.5_d4"

#levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1")
#data_version = "_d2_d2.5_d3"

#levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA")
#data_version = "_d2_d2.5_d3_d3.5_d4"

#levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA")
#data_version = "_d2_d2.5_d3_d3.5_d4"

levels_sels = c("day2_beforeRA", "day2.5_RA", "day3_RA.rep1",  "day3_RA.rep2",
                "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
data_version = "_d2_d2.5_d3_d3.5_d4_d5_d6"

names(cols) = levels
cols_sel = cols[match(levels_sels, names(cols))]

outDir = paste0(resDir, '/RA_symetryBreaking/signaling_pathway_activity', data_version)
system(paste0('mkdir -p ', outDir))

Idents(aa) = factor(aa$condition, levels = levels)
aa = subset(aa, idents = levels_sels)

########################################################
########################################################
# Section 0: prepare the scRNA-seq samples  
# select only day3, day3.5 and day4
########################################################
########################################################
Prepare_selected_conditions = FALSE
if(Prepare_selected_conditions){
  # aa = readRDS(file = paste0(RdataDir, 
  #                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
  #                           'cellCycleScoring_annot.v2_newUMAP_clusters_noRAday6_noRAday2.5_',
  #                           species, version.analysis, '.rds'))
  
  
  ##########################################
  # subset the time points
  ##########################################
  Idents(aa) = factor(aa$condition, levels = levels)
  aa = subset(aa, idents = levels_sels)
  
  #Idents(aa) = aa$condition
  # aa = subset(aa, cells = colnames(aa)[which(aa$condition != 'day5_RA')])
  
  # rerun the umap 
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 50, min.dist = 0.1)
  
  
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
                            'cellCycleScoring_annot.v2_newUMAP_clusters', data_version, '_',
                            species, version.analysis, '.rds'))
  
}

########################################################
########################################################
# Section I: test global pathway analysis 
#  - e.g. decoupleR, Pagoda2 (optional)
# decoupleR original code from https://saezlab.github.io/decoupleR/articles/pw_sc.html
# and https://www.sc-best-practices.org/conditions/gsea_pathway.html
########################################################
########################################################
library(Seurat)
library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters', data_version, '_',
                           species, version.analysis, '.rds'))

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_', 
                           "d3.d3.5.d4_mNT_scRNAseq_R13547_10x_mNT_20220813.rds"
                           ))

# downsample for each time point 
Idents(aa) = aa$condition
aa = subset(aa, downsample = 2500)

DimPlot(aa, group.by = 'clusters', label = TRUE)

##########################################
# test progeny 
##########################################
# PROGENy is a comprehensive resource containing a curated collection of pathways and their target genes, 
# with weights for each interaction
net <- get_progeny(organism = 'mouse', top = 100)
net

# Activity inference with Weighted Mean
# Extract the normalized log-transformed counts
mat <- as.matrix(aa@assays$RNA@data)

# Run wmean
tic()
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',
                  .mor='weight', times = 100, minsize = 5)
toc()

acts

## Visualization 
# Extract norm_wmean and store it in pathwayswmean in data
aa[['pathwayswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = aa) <- "pathwayswmean"

# Scale the data
aa <- ScaleData(aa)
aa@assays$pathwayswmean@data <- aa@assays$pathwayswmean@scale.data

p1 <- DimPlot(aa, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(aa, features = c("Trail")) & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('Trail activity')
p1 | p2

FeaturePlot(aa, features = rownames(aa)) & 
  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')

ggsave(filename = paste0(outDir, '/progeny_14_signalingPathways.pdf'), 
       width = 20, height = 12)


## Exploration 
# Extract activities from object as a long dataframe
Idents(aa) = aa$condition
DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)

df <- t(as.matrix(aa@assays$pathwayswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(aa)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

# Plot
library(pheatmap)
pheatmap::pheatmap(t(top_acts_mat), 
                   border_color = NA, 
                   cluster_cols = FALSE,
                   color=my_color, breaks = my_breaks, 
                   filename = paste0(outDir, '/progeny_14_signalingPathways_summary_time.pdf'), 
                   width = 6, height = 6) 


##########################################
# test TF activity with decoupleR (original code https://saezlab.github.io/decoupleR/articles/tf_sc.html)  
##########################################
net <- get_dorothea(organism='mouse', levels=c('A', 'B', 'C'))
net

xx = (data.frame(net[grep('Smad|Tead', net$source), ]))
write.csv(xx, file = paste0(outDir, '/targets_decoupleR_for_Smad_Tead_activtiyInference.csv'), 
          quote = FALSE, row.names = FALSE)

# Extract the normalized log-transformed counts
mat <- as.matrix(aa@assays$RNA@data)

# Run wmean
tic()
acts <- run_wmean(mat=mat, net=net, .source='source', .target='target',
                  .mor='mor', times = 100, minsize = 5)
toc()

acts

# Extract norm_wmean and store it in tfswmean in pbmc
aa[['tfswmean']] <- acts %>%
  filter(statistic == 'norm_wmean') %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = aa) <- "tfswmean"

# Scale the data
aa <- ScaleData(aa)
aa@assays$tfswmean@data <- aa@assays$tfswmean@scale.data


DimPlot(aa, reduction = "umap", group.by =  "Phase",  label = TRUE, pt.size = 0.5, label.size = 6) + 
  ggtitle('time points')
ggsave(filename = paste0(outDir,  '/cellcycle_umap.pdf'), 
       width = 10, height = 6)

p0 = DimPlot(aa, reduction = "umap", group.by = 'condition',  label = TRUE, pt.size = 0.5, label.size = 6) + 
  ggtitle('time points')

p1 <- DimPlot(aa, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6) + 
  NoLegend() + ggtitle('clusters')
p0 | p1

ggsave(filename = paste0(outDir,  '/timePoints_clusters.pdf'), 
       width = 12, height = 6)


DefaultAssay(object = aa) <- "tfswmean"
n_tfs <- 70
# Extract activities from object as a long dataframe
Idents(aa) = aa$clusters

df <- t(as.matrix(aa@assays$tfswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(aa)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap::pheatmap(t(top_acts_mat), border_color = NA, color=my_color, breaks = my_breaks,
         filename = paste0(outDir, '/decoupleR_TF.activity_summary_clusters.pdf'), 
         width = 6, height = 16) 


Idents(aa) = aa$condition

df <- t(as.matrix(aa@assays$tfswmean@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(aa)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Get top tfs with more variable means across clusters
tfs <- df %>%
  group_by(source) %>%
  summarise(std = sd(mean)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)

# Subset long data frame to top tfs and transform to wide matrix
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pheatmap::pheatmap(t(top_acts_mat), border_color = NA, color=my_color, breaks = my_breaks,
                   cluster_cols = FALSE,
                   filename = paste0(outDir, '/decoupleR_TF.activity_summary_time.pdf'), 
                   width = 4, height = 12) 

#for (gene in c('Foxa2', 'Pax6', 'Sox2', 'Sox1'))
#for(gene in colnames(top_acts_mat))
for(gene in rownames(aa))
{
  cat('gene -- ', gene, '\n')
  # gene = 'Smad3'
  p2 <- (FeaturePlot(aa, features = gene) & 
           scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) + 
    ggtitle(paste0(gene, ' activity'))
  DefaultAssay(object = aa) <- "RNA" 
  p3 <- FeaturePlot(aa, features = c(gene)) + ggtitle(paste0(gene, ' expression'))
  DefaultAssay(object = aa) <- "tfswmean"
  p2 | p3
  
  ggsave(filename = paste0(outDir,  '/decoupleR_TF_activity_expression_', gene, '.pdf'), 
         width = 12, height = 6)
  
  
}

saveRDS(aa, file = paste0(RdataDir, 
                      'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                      'cellCycleScoring_annot.v1_savedUMAP.subs.v2_savedTFactivities.decoupleR.rds'))

##########################################
# try ReactomeGSA for global pathway analysis 
# https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/analysing-scRNAseq.html
# 
##########################################
use_ReactomeGSA = FALSE
if(use_ReactomeGSA){
  require(ReactomeGSA)
  require(tictoc)
  
  cell_ids  = 'condition'
  
  if(cell_ids == 'condition')  {
    Idents(aa) = aa$condition
  }else{
    Idents(aa) = aa$clusters
  }
  
  tic()
  gsva_result <- analyse_sc_clusters(object = aa, verbose = TRUE)
  toc()
  
  # The resulting object is a standard ReactomeAnalysisResult object.
  gsva_result
  
  # pathways returns the pathway-level expression values per cell cluster
  pathway_expression <- pathways(gsva_result)
  
  # simplify the column names by removing the default dataset identifier
  colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))
  pathway_expression[1:3, ]
  
  # A simple approach to find the most relevant pathways is to assess the maximum difference 
  # in expression for every pathway:
  
  # find the maximum differently expressed pathway
  max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(abs(values)), max = max(abs(values))))
  }))
  
  max_difference$diff <- max_difference$max - max_difference$min
  
  # sort based on the difference
  max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]
  head(max_difference)
  
  plot(max_difference$max, max_difference$diff, cex = 0.5)
  abline(v = 0.2)
  abline(h = 0.02)
  sels = which(max_difference$max>=0.25 |max_difference$diff>0.05)
  
  cat(length(sels), ' pathways passing thresholds \n')
  
  #relevant_pathways <- c("R-HSA-983170", "R-HSA-388841", "R-HSA-2132295", "R-HSA-983705", "R-HSA-5690714")
  jj = intersect(grep('FGF|WNT|BMP|TGF|Retinoic Acid|Hippo|EGFR|NOTCH|GLI|SHH|Hh', max_difference$name),
                 grep('Signaling|signaling|Hippo|GLI|SHH|Hh|ERK', max_difference$name))
  
  max_difference_jj = max_difference[unique(c(jj)), ]
  max_difference = max_difference[unique(c(sels, jj)), ]
  
  gsva_mat = pathway_expression[match(max_difference$name, pathway_expression$Name), ]
  rownames(gsva_mat) = gsva_mat$Name
  gsva_mat = gsva_mat[, -c(1)]
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  pheatmap::pheatmap(gsva_mat, border_color = NA, color=my_color, 
           #breaks = my_breaks,
           scale = 'row',
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize = 7,
           filename = paste0(outDir, '/ReactomeGSA_global_pathwayActivity_', cell_ids, '_many.pdf'), 
           width = 10, height = 28) 
  
  # selected pathways
  gsva_mat = pathway_expression[match(max_difference_jj$name, pathway_expression$Name), ]
  rownames(gsva_mat) = gsva_mat$Name
  gsva_mat = gsva_mat[, -c(1)]
  
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  pheatmap::pheatmap(gsva_mat, border_color = NA, color=my_color, 
           #breaks = my_breaks,
           scale = 'row',
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           fontsize = 10,
           filename = paste0(outDir, '/ReactomeGSA_global_pathwayActivity_', cell_ids, '_selected.pdf'), 
           width = 10, height = 14) 
  
}

##########################################
# try different method for signaling pathways with decoupleR
# original code from https://www.sc-best-practices.org/conditions/gsea_pathway.html
##########################################
Test_DecoupleR_reactome = FALSE
if(Test_DecoupleR_reactome){
  #reactome = gmt_to_decoupler("c2.cp.reactome.v7.5.1.symbols.gmt")
  msigdb = get_resource("MSigDB")
  
  # Get reactome pathways
  reactome = msigdb %>% filter(collection == 'reactome_pathways')
  
  # Filter duplicates
  reactome = reactome %>% dplyr::select(geneset, genesymbol) %>% 
    distinct()
  
  # Filtering genesets to match behaviour of fgsea
  geneset_size = table(as.data.frame(reactome)[,1])
  geneset_sel = names(geneset_size)[which((geneset_size > 15) & (geneset_size < 500))]
  
  #res_gsea <- run_fgsea(mat, network, .source='source', .target='target', nproc=1, minsize = 0)
  
  library(tictoc)
  
  # run aucel for each cell
  mat <- as.matrix(aa@assays$RNA@data)
  rownames(mat) = toupper(rownames(mat))
  
  tic() # run aucell for each cell
  aucel = run_aucell(mat = mat, network = reactome, 
                     .source="geneset",
                     .target="genesymbol",
                     minsize = 5
  )
  toc()
  
  
  # Extract norm_wmean and store it in tfswmean in pbmc
  aa[['aucell']] <- aucel %>%
    filter(statistic == 'aucell') %>%
    pivot_wider(id_cols = 'source', names_from = 'condition',
                values_from = 'score') %>%
    column_to_rownames('source') %>%
    Seurat::CreateAssayObject(.)
  
  # Change assay
  DefaultAssay(object = aa) <- "aucell"
  
  # Scale the data
  aa <- ScaleData(aa)
  aa@assays$aucell@data <- aa@assays$aucell@scale.data
  
  Idents(aa) = aa$clusters
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)
  
  df <- t(as.matrix(aa@assays$aucell@data)) %>%
    as.data.frame() %>%
    mutate(cluster = Idents(aa)) %>%
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
    group_by(cluster, source) %>%
    summarise(mean = mean(score))
  
  # Get top tfs with more variable means across clusters
  ntop = 50
  top_SP <- df %>%
    group_by(source) %>%
    summarise(std = sd(mean)) %>%
    arrange(-abs(std)) %>%
    #summarise(mean_mean = mean(mean)) %>%
    #arrange(-abs(mean_mean)) %>%
    head(ntop) %>%
    pull(source)
  
  # Subset long data frame to top tfs and transform to wide matrix
  top_aucell_mat <- df %>%
    filter(source %in% top_SP) %>%
    pivot_wider(id_cols = 'cluster', names_from = 'source',
                values_from = 'mean') %>%
    column_to_rownames('cluster') %>%
    as.matrix()
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 3, length.out=floor(palette_length/2)))
  
  # Plot
  pheatmap(t(top_aucell_mat), border_color = NA, color=my_color, breaks = my_breaks,
           filename = paste0(outDir, '/decoupleR_REACTOME_summary_clusters.pdf'), 
           width = 20, height = 20) 
  
  # Transform to wide matrix
  top_accell_mat <- df %>%
    pivot_wider(id_cols = 'cluster', names_from = 'source',
                values_from = 'mean') %>%
    column_to_rownames('cluster') %>%
    as.matrix()
  
  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 2, length.out=floor(palette_length/2)))
  
  # Plot
  pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks, 
           filename = paste0(outDir, '/progeny_14_signalingPathways_summary_clusters.pdf'), 
           width = 10, height = 6) 
  
  sp = "REACTOME-SIGNALING-BY-RETINOIC-ACID"
  sp = "REACTOME-SIGNALING-BY-FGFR" 
  (FeaturePlot(aa, features = c(sp)) &
      scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red'))
  
  
}

########################################################
########################################################
# Section II : ligand-receptor anlaysis
# e.g. LIANA and NicheNet with defined subpopulation
# test LIANA and NicheNet using clusters 
########################################################
########################################################

## test LIANA
source("functions_ligandReceptor_analysis.R")
aa$celltypes = aa$clusters

## test NicheNet
source("functions_ligandReceptor_analysis.R")


########################################################
########################################################
# Section III: native way of checking gene examples of pathways of interest
# (RA, WNT, BMP, TGF-beta, FGF, Noch, Yap, SHH, PI3K)
########################################################
########################################################

##########################################
# plot gene examples of selected signaling pathways
##########################################
curate.geneList.signalingPathways = FALSE
if(curate.geneList.signalingPathways){
  sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v2.rds'))
  #sps = unique(sps$gene)
  xx = read.table('../data/annotations/GO_term_summary_RAR_signalingPathway.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  xx = unique(xx[,2])
  xx = data.frame(gene = xx, pathway = rep('RA', length(xx)))
  sps = rbind(sps, xx)
  
  
  xx = read.table('../data/annotations/GO_term_summary_TGFb.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  xx = unique(xx[,2])
  xx = data.frame(gene = xx, pathway = rep('TGF_beta', length(xx)))
  sps = rbind(sps, xx)
  
  xx = read.table('../data/annotations/GO_term_summary_NOTCH.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  xx = unique(xx[,2])
  xx = data.frame(gene = xx, pathway = rep('NOTCH', length(xx)))
  sps = rbind(sps, xx)
  
  xx = read.table('../data/annotations/GO_term_summary_PI3K.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  xx = unique(xx[,2])
  xx = data.frame(gene = xx, pathway = rep('PI3K', length(xx)))
  sps = rbind(sps, xx)
  
  xx = read.table('../data/annotations/GO_term_summary_Hippo.txt', header = TRUE, sep = '\t', 
                  row.names = NULL)
  xx = unique(xx[,2])
  xx = data.frame(gene = xx, pathway = rep('Hippo', length(xx)))
  sps = rbind(sps, xx)
  
  saveRDS(sps, file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
  
  
}


saveDir = paste0(outDir, '/geneExamples_pathways/')
system(paste0('mkdir -p ', saveDir))

DefaultAssay(aa) = 'RNA'

cell_ids  = 'clusters'

sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

#sps = intersect(sps, markers$gene)
markers_saved = c()

for(cell_ids in c('condition', 'clusters'))
{
  # cell_ids = 'condition'
  if(cell_ids == 'condition')  {
    Idents(aa) = aa$condition
  }else{
    Idents(aa) = aa$clusters
  }
  
  markers = FindAllMarkers(aa, only.pos = TRUE,  
                           #min.pct = 0.25, 
                           logfc.threshold = 0.25)
  
  markers_saved = unique(c(markers_saved, markers$gene))
  markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top10
  DoHeatmap(aa, features = top10$gene) + NoLegend()
  
  ggsave(filename = paste0(saveDir,  'heatmap_markerGenes_top20_by.', cell_ids, '.pdf'), 
         width = 12, height = 20)
  
  DoHeatmap(aa, features = intersect(sps, markers$gene)) + NoLegend()
  ggsave(filename = paste0(saveDir,  'heatmap_markerGenes_signalingGenes_by.', cell_ids, '.pdf'), 
         width = 12, height = 25)
  
  tfs_sel = intersect(tfs, markers$gene)
  DoHeatmap(aa, features = tfs_sel) + NoLegend()
  ggsave(filename = paste0(saveDir,  'heatmap_markerGenes_TFs_by.', cell_ids, '.pdf'), 
         width = 12, height = 16)
  
}

sps_sels = intersect(sps, markers_saved)
# feature plots 
pdfname = paste0(saveDir, "gene_examples_signalingPathways.pdf")
pdf(pdfname, width=20, height = 6)
par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

#for(n in 1:20)
for(n in 1:length(sps_sels))
{
  # n = 1
  gg = sps_sels[n];
  # gg = "Foxa2"
  cat(n, '--', gg, '\n')
  
  p1 = FeaturePlot(aa, features = gg)
  p2 = VlnPlot(aa, features = gg, group.by = 'condition', pt.size = 0.01)  + NoLegend()  
    #scale_fill_manual(values=c("blue", 'green', "red",  "blue", 
    #                           'red', 'green', 'blue', 'green', 'red', 'green', 'green'))
  p3 = VlnPlot(aa, features = gg, group.by = 'clusters', pt.size = 0.01) + NoLegend()
  plot(p1 | p2| p3)
  
}

dev.off()

tfs = unique(c(tfs, "Zfp703"))
tfs_sel = intersect(tfs, markers_saved)

pdfname = paste0(saveDir, "gene_examples_TFs.pdf")
pdf(pdfname, width=24, height = 6)
par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

#for(n in 1:20)
for(n in 1:length(tfs_sel))
{
  # n = 1
  gg = sps_sels[n];
  # gg = "Foxa2"
  cat(n, '--', gg, '\n')
  
  p1 = FeaturePlot(aa, features = gg)
  p2 = VlnPlot(aa, features = gg, group.by = 'condition', pt.size = 0.01)  + NoLegend()  
  #scale_fill_manual(values=c("blue", 'green', "red",  "blue", 
  #                           'red', 'green', 'blue', 'green', 'red', 'green', 'green'))
  p3 = VlnPlot(aa, features = gg, group.by = 'clusters', pt.size = 0.01) + NoLegend()
  plot(p1 | p2| p3)
  
}

dev.off()

