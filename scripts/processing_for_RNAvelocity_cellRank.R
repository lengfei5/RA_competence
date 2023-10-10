##########################################################################
##########################################################################
# Project: RA competence  
# Script purpose: prepare files for RNA velocity using alevin 
# https://github.com/csoneson/rna_velocity_quant
# https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Oct  7 14:09:17 2022
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
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA", "day4_RA", "day5_RA", "day6_RA")
cols_sel = cols[match(levels_sels, names(cols))]

#version.analysis = paste0(version.analysis, '_ES.beforeRA.and.RA')

outDir = paste0(resDir, '/RA_symetryBreaking/cellRank_scvelo_d2.beforeRA.to.d6_rm.weired.clusters/')
system(paste0('mkdir -p ', outDir))


##########################################
# preapre the intron and transcript reference fasta files
##########################################
# gtf and fa files of mouse using eisaR
# we load the eisaR package and extract a GRanges object 
# containing the genomic coordinates of each annotated transcript and intron. 
# In this example, we use the ‘separate’ approach to define introns separately for each transcript, 
# and add a flank length of 90nt to each intron.
genomeDir = '/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/'
gtf <- "/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/Mus_musculus.GRCm38.87.gtf"
genome.file = "/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/Mus_musculus.GRCm38.dna.toplevel.fa"

grl <- eisaR::getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "intron"), 
  intronType = "separate", 
  flankLength = 90L, 
  joinOverlappingIntrons = FALSE, 
  verbose = TRUE
)

grl[4:6]


genome <- Biostrings::readDNAStringSet(genome.file)

names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = grl
)

Biostrings::writeXStringSet(
  seqs, filepath = paste0(genomeDir,  "Mus.GRCm38.87.annotation.expanded.fa")
)

eisaR::exportToGtf(
  grl, 
  filepath = paste0(genomeDir, "Mus.GRCm38.87.annotation.expanded.gtf")
)

head(metadata(grl)$corrgene)

write.table(
  metadata(grl)$corrgene, 
  file = paste0(genomeDir, "Mus.GRCm38.87.annotation.expanded.features.tsv"),
  row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)


df <- eisaR::getTx2Gene(
  grl, filepath = paste0(genomeDir, "Mus.GRCm38.87.annotation.expanded.tx2gene.tsv")
)


##########################################
# Index the reference features 
# https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
##########################################
## run the salmon index 
## bash script

## we also create a linked transcriptome with tximeta. 
## This allows tximeta to recognize the reference annotation when reading the alevin quantification, 
## and automatically annotate the resulting SummarizedExperiment object.
tximeta::makeLinkedTxome(
  indexDir = "/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/salmon_mm10_annotation.expanded.sidx", 
  source = "Ensembl", 
  genome = "GRCm38", 
  organism = "Mus musculus", 
  release = "87", 
  fasta = paste0(genomeDir, "alevin_velocity/Mus.GRCm38.87.annotation.expanded.fa"), 
  gtf =paste0(genomeDir, "alevin_velocity/Mus.GRCm38.87.annotation.expanded.gtf"),
  write = TRUE, 
  jsonFile = paste0(genomeDir, "alevin_velocity/Mus.GRCm38.87.annotation.expanded_v2.json")
)

rjson::fromJSON(file = paste0(genomeDir, "alevin_velocity/Mus.GRCm38.87.annotation.expanded_v2.json"))


########################################################
########################################################
# Section : process the alevin output and prepare the data for scvelo and cellrank
# 
########################################################
########################################################
aa =  readRDS(file = paste0(RdataDir, 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_',
                            'regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', species, version.analysis, '.rds'))

Idents(aa) = factor(aa$condition, levels = levels)

table(aa$condition)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)

ggsave(filename = paste0(outDir, 'UMAP_rmDoublet_rmRiboMT_regressed.nCounts_annot.v1_',
                         'beforeSubsetting_RAsymmetryBreaking.pdf'), width = 10, height = 8)

aa = subset(aa, idents = levels_sels)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)

ggsave(filename = paste0(outDir, 'UMAP_rmDoublet_rmRiboMT_regressed.nCounts_annot.v1_',
                         'subsetting.RAsymmetryBreaking.onlyday3rep1',
                         version.analysis, '.pdf'), 
       width = 10, height = 8)

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                          'cellCycleScoring_annot.v2_',
                          species, version.analysis, '.rds'))

# saveRDS(aa, file = paste0(RdataDir, 
#                           'seuratObject_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
#                           'cellCycleScoring_annot.v2_', 'RA.symmetry.breaking_d2.to.d6_',
#                           species, version.analysis, '.rds'))


##########################################
# Import abundances into R with tximeta
# https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
##########################################
library(Seurat)
library(tximport)

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_',
                           species, version.analysis, '.rds'))


design = read.csv(file = paste0(dataDir, '/sampleInfos.csv'))
colnames(design) = c('sampleID', 'condition')

design = design[, c(1:2)]
design$condition = gsub('day0', 'day0_beforeRA', design$condition)

design = design[match(levels_sels, design$condition), ]

design$condition = gsub('day0_beforeRA', 'day0', design$condition)

##########################################
# prepare gene annotation 
##########################################
ggs = read.csv(paste0("../R13547_10x/cellRanger_outs/196324/outs/filtered_feature_bc_matrix", "/features.tsv.gz"), 
               header = F, stringsAsFactors = F, sep = '\t')
## make unique gene names
g = ggs
rm(ggs)
g$name = g$V2
gg.counts = table(g$V2)
gg.dup = names(gg.counts)[which(gg.counts>1)]
index.dup = which(!is.na(match(g$V2, gg.dup)))
g$name[index.dup] = paste0(g$V2[index.dup], '_', g$V1[index.dup])

cg <- read.delim(paste0(genomeDir, "alevin_velocity/Mus.GRCm38.87.annotation.expanded.features.tsv"),
                 header = TRUE, as.is = TRUE)

## Rename the 'intron' column 'unspliced' to make assay names compatible with scVelo
colnames(cg)[colnames(cg) == "intron"] <- "unspliced"
cg$name = g$name[match(cg$spliced, g$V1)]

# import data from cellranger output
for(n in 1:nrow(design))
{
  # n = 1
  cat(n, ' : ', design$sampleID[n], '--', design$condition[n], '\n')
  
  topdir = paste0(dataDir, '/RNAvelocity_alevin_outs/alevin_out_',  design$sampleID[n], 
                  '/alevin/')
  
  # https://combine-lab.github.io/alevin-tutorial/2018/alevin-seurat/
  txi <- tximport(paste0(topdir, 'quants_mat.gz'), type="alevin")
  
  count.data = txi$counts
  rm(txi);
  
  meta = data.frame(cell.id = paste0(colnames(count.data), '-1', '_', design$condition[n], '_', 
                                     design$sampleID[n]), 
                    condition = design$condition[n])
  cell2keep = !is.na(match(meta$cell.id, aa$cell.id)) 
  meta$cell2keep = cell2keep
  rownames(meta) = colnames(count.data)
  
  mm1 = match(cg$spliced, rownames(count.data))
  kk1 = which(!is.na(mm1))
  xx_spliced = count.data[mm1[kk1], ]
  rownames(xx_spliced) = cg$name[kk1]
  
  mm2 = match(cg$unspliced, rownames(count.data))
  kk2 = which(!is.na(mm2))
  xx_unspliced = count.data[mm2[kk2], ]
  rownames(xx_unspliced) = cg$name[kk2]
  
  rm(count.data)
  
  # creat seurat object
  aa_spliced = CreateSeuratObject(counts = xx_spliced[, cell2keep],
                          meta.data = meta[cell2keep, ], 
                          min.cells = 5, min.features = 10)
  rm(xx_spliced)
  
  aa_unspliced = CreateSeuratObject(counts = xx_unspliced[, cell2keep],
                                  meta.data = meta[cell2keep, ], 
                                  min.cells = 5, min.features = 10)
  rm(xx_unspliced)
  rm(meta)
  
  if(n == 1) {
    spliced = aa_spliced
    unspliced = aa_unspliced
    
  }else{
    spliced = merge(spliced, aa_spliced)
    unspliced = merge(unspliced, aa_unspliced)
  }
  
  rm(aa_spliced)
  rm(aa_unspliced)
  
}

save(spliced, unspliced, 
     file = paste0(RdataDir, 'seuratObject_RNAvelocity_alevin_spliced_unspliced_', 
                   species, version.analysis, '.Rdata'))

##########################################
# intersect cells in the analysis with spliced and unspliced matrix
# save seurat object for scvelo
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_',
                           species, version.analysis, '.rds'))

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_',
                           species, version.analysis, '.rds'))


table(aa$condition)

# rerun the umap 
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(aa, ndims = 50)

Idents(aa) = aa$condition

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)

ggsave(filename = paste0(outDir, '/UMAP_RAsymmetryBreaking.onlyday3rep1_5000HVGs_noweighted.byvarPCA_',
                         '30pcs_30neighbors_minDist0.1', version.analysis, '.pdf'), 
       width = 10, height = 6)

# quickly run clustering
ElbowPlot(aa, ndims = 50)
aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p1 + p2

ggsave(filename = paste0(outDir, 'UMAP_RAsymmetryBreaking.onlyday3rep1_timePoints_clustering.res0.5',
                        version.analysis, '.pdf'), 
       width = 18, height = 8)

Discard_weired.clusters = FALSE
if(Discard_weired.clusters){
  
  aa = subset(aa, cells = colnames(aa)[which(aa$seurat_clusters != 11 & aa$seurat_clusters != 7)])
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p1 + p2
 
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 50, min.dist = 0.1)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)
  
  ElbowPlot(aa, ndims = 50)
  aa <- FindNeighbors(aa, dims = 1:30)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p1 + p2
  
  ggsave(filename = paste0(outDir, 'UMAP_RAsymmetryBreaking.onlyday3rep1_timePoints_clustering.res0.7',
                           version.analysis, '.pdf'), 
         width = 18, height = 8)
  
}

aa$clusters = aa$seurat_clusters

#jj = which(aa$celltypes == 'FP'|aa$celltypes == 'Neurons'|aa$celltypes == 'NP_RA')
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)

aa$celltypes = as.character(aa$clusters)

#aa$celltypes[which(aa$clusters == '6')] = 'FP'
#aa$celltypes[which(aa$clusters == '10')] = 'Neurons'
#aa$celltypes[which(aa$clusters == '8')] = 'NP'

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)

aa$time = aa$condition
aa$time = gsub('day','',aa$time)
aa$time = gsub('_RA','',aa$time)
aa$time = gsub('.rep1','',aa$time)
aa$time = gsub('_beforeRA', '', aa$time)

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                          'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d6_',
                          species, version.analysis, '.rds'))

##########################################
# ignore RA_d6 
##########################################
NotConsider.RA_day6 = FALSE
if(NotConsider.RA_day6){
  aa = readRDS(file = paste0(RdataDir, 
                             'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                             'cellCycleScoring_annot.v2_newUMAP_clusters_time_',
                             species, version.analysis, '.rds'))
  
  Idents(aa) = aa$condition
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)
  
  ## remove the neurons 
  aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != 'Neurons')])
  
  # rerun the umap
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.1)
  
  #DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols_sel, raster=FALSE)
  #ggsave(filename = paste0(outDir, 'UMAP_RAsymmetryBreaking.onlyday3rep1_3000HVGs_noweighted.byvarPCA_',
  #                         '30pcs_30neighbors_minDist0.1.pdf'), 
  #       width = 10, height = 8)
  
  # quickly run clustering
  #ElbowPlot(aa, ndims = 50)
  aa <- FindNeighbors(aa, dims = 1:30)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p1 + p2
  
  ggsave(filename = paste0(outDir, 
                           'UMAP_RAsymmetryBreaking.onlyday3rep1_timePoints_clustering.res0.7_noRAday6',
                           version.analysis, '.pdf'), 
         width = 18, height = 8)
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)
  
  aa$clusters = aa$seurat_clusters
  
  #xx = subset(aa, idents = "7", invert = TRUE)
  #aa = xx
  #xx = aa
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  aa$clusters[which(aa$clusters == "11")] = "4"
  aa$clusters[which(aa$clusters == "13")] = "1"
  
  aa$clusters[which(aa$clusters == "10")] = "5"
  aa$clusters[which(aa$clusters == "12")] = "5"
  
  aa$clusters = as.character(aa$clusters)
  aa$clusters[which(aa$clusters == "8")] = "FP"
  aa$clusters[which(aa$clusters == "6")] = "NP"
  
  aa$clusters = as.factor(aa$clusters)
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  
  # aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
  # 
  # ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  # aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  # ElbowPlot(aa, ndims = 50)
  # 
  # aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
  # 
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  p1 + p2
  
  ggsave(filename = paste0(outDir, 
                           'UMAP_RAsymmetryBreaking.onlyday3rep1_timePoints_clustering.res0.7_noNeurons',
                           '_clened', version.analysis, '.pdf'), 
         width = 18, height = 8)
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                            'cellCycleScoring_annot.v2_newUMAP_clusters_time_noNeurons_',
                            species, version.analysis, '.rds'))
  
  
  ##########################################
  # include only day2_before RA, otherwise it is too challenging due to the big gap between 
  # early time points
  ##########################################
  aa = readRDS(paste0(RdataDir, 
                      'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                      'cellCycleScoring_annot.v2_newUMAP_clusters_time_noNeurons_',
                      species, version.analysis, '.rds'))
  
  Idents(aa) = aa$condition
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)
  
  ## remove the neurons 
  aa = subset(aa, cells = colnames(aa)[which(aa$condition != 'day0_beforeRA' & 
                                               aa$condition != 'day1_beforeRA')])
  
  # rerun the umap
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.1)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  FeaturePlot(aa, features = c('Pax6', 'Foxa2'))
  
  aa <- FindNeighbors(aa, dims = 1:30)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  p1 + p2
  
  ggsave(filename = paste0(outDir, 
                           'UMAP_RAsymmetryBreaking.onlyday3rep1_timePoints_clustering.res0.7_',
                           'd2.beforRA_noRAday6',
                           version.analysis, '.pdf'), 
         width = 18, height = 8)
  
  
  aa$clusters = aa$seurat_clusters
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  FeaturePlot(aa, features = c('Pax6', 'Foxa2'))
  
  #aa$clusters[which(aa$clusters == "11")] = "2"
  #aa$clusters[which(aa$clusters == "13")] = "1"
  
  #aa$clusters[which(aa$clusters == "10")] = "5"
  #aa$clusters[which(aa$clusters == "12")] = "0"
  
  aa$clusters = as.character(aa$clusters)
  aa$clusters[which(aa$clusters == "6"|aa$clusters == '7')] = "FP"
  aa$clusters[which(aa$clusters == "5")] = "NP"
  
  aa$clusters = as.factor(aa$clusters)
  
  #p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  #p2 = FeaturePlot(aa, features = c('Pax6', 'Foxa2'))
  #VlnPlot(aa, features = c('Pax6', 'Foxa2'), group.by = 'clusters')
  # aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
  # 
  # ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  # aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  # ElbowPlot(aa, ndims = 50)
  # 
  # aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
  # 
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  p1 + p2
  
  ggsave(filename = paste0(outDir, 
                           'UMAP_RAsymmetryBreaking.onlyday3rep1_timePoints_clustering.res0.7_noNeurons',
                           '_d2.beforeRA_cleaned', version.analysis, '.pdf'), 
         width = 18, height = 8)
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                            'cellCycleScoring_annot.v2_newUMAP_clusters_time_noNeurons_d2.beforeRA_',
                            species, version.analysis, '.rds'))
  
  ## try to remove strange clusters 8, 9, 13
  xx = subset(aa, cells = colnames(aa)[which(aa$clusters != '8' & 
                                               aa$clusters != '9' &
                                               aa$clusters != '13')])
  aa = xx
  rm(xx)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
  p1 + p2
  
  ggsave(filename = paste0(outDir, 
                           'UMAP_RAsymmetryBreaking.onlyday3rep1_timePoints_clustering.res0.7_noNeurons',
                           '_d2.beforeRA_cleaned_rmSuspeciousClusters', version.analysis, '.pdf'), 
         width = 18, height = 8)
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                            'cellCycleScoring_annot.v2_newUMAP_clusters_time_noNeurons_d2.beforeRA_',
                            'rmSuspeciousClusters_',
                            species, version.analysis, '.rds'))
  
  
}

########################################################
########################################################
# Section : prepare the files for scVelo and cellRank
# 
########################################################
########################################################
# aa = readRDS(file = paste0(RdataDir, 
#                            'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
#                            'cellCycleScoring_annot.v2_newUMAP_clusters_time_noNeurons_d2.beforeRA_',
#                            'rmSuspeciousClusters_',
#                            species, version.analysis, '.rds'))
#aa = readRDS(file = paste0(RdataDir, 'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_',
#                           'regressout.nCounts_cellCycleScoring_annot.v2_mNT_scRNAseq_R13547_10x_mNT_20220813.rds'))
#aa = readRDS(file = paste0(RdataDir, 
#                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
#                           'cellCycleScoring_annot.v2_newUMAP_clusters_time_updated20231005_',
#                           species, version.analysis, '.rds'))
load(file = paste0(RdataDir, 'seuratObject_RNAvelocity_alevin_spliced_unspliced_', 
                   'mNT_scRNAseq_R13547_10x_mNT_20220813.Rdata'))

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_time_d2.to.d6_',
                           species, version.analysis, '.rds'))

load(file = paste0(RdataDir, 
                   'seuratObject_RNAvelocity_alevin_spliced_unspliced_',
                   'mNT_scRNAseq_R13547_10x_mNT_20220813_ES.beforeRA.and.RA.Rdata'))

FeaturePlot(aa, features = c('Pax6', 'Foxa2', 'Sox2'))


p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'clusters', raster=FALSE)
p1 + p2

p0 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)

p2 + p0


aa$clusters = aa$seurat_clusters
aa$celltypes = aa$clusters
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE)

Discard_cellCycle.corrrelatedGenes = FALSE
if(Discard_cellCycle.corrrelatedGenes){
  library(scater)
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
  
  print(intersect(genes_discard, gene_examples))
  
  aa = subset(aa, features = setdiff(rownames(aa), genes_discard))
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000) # find subset-specific HVGs
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.3)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
     
}

## select features shares by spliced and unspliced
features = intersect(rownames(aa), intersect(rownames(spliced), rownames(unspliced)))
spliced = subset(spliced, features = features)
unspliced = subset(unspliced, features = features)

# subsetting aa with the same cell as spliced and unspliced
cells.shared = intersect(spliced$cell.id, aa$cell.id)
mnt = subset(aa, cells = colnames(aa)[match(cells.shared, aa$cell.id)])

counts = GetAssayData(spliced, slot = 'counts')
counts = counts[, match(mnt$cell.id, spliced$cell.id)]
colnames(counts) = colnames(mnt)

# subsetting mnt with the same genes as spliced and unspliced
DefaultAssay(mnt) = 'RNA'
mnt = subset(mnt, features = rownames(counts))

counts = counts[match(rownames(mnt), rownames(counts)), ]
mnt[["spliced"]] <- CreateAssayObject(counts = counts)
rm(spliced)

counts = GetAssayData(unspliced, slot = 'counts')
counts = counts[, match(mnt$cell.id, unspliced$cell.id)]
colnames(counts) = colnames(mnt)
counts = counts[match(rownames(mnt), rownames(counts)), ]
mnt[["unspliced"]]<- CreateAssayObject(counts = counts)
rm(unspliced)

# try to save mulitple assays 
# https://github.com/mojaveazure/seurat-disk/issues/21

#library(SeuratData)
library(SeuratDisk)

VariableFeatures(mnt) = NULL
#mnt@assays$RNA@scale.data = NULL
#mnt@assays$RNA@data = NULL

DefaultAssay(mnt) = 'RNA'
mnt = DietSeurat(mnt, counts = TRUE, data = TRUE,
                 scale.data = FALSE,
                 features = rownames(mnt), 
                 assays = c('RNA', 'spliced', 'unspliced'), 
                 dimreducs = c('umap'), graphs = NULL, 
                 misc = TRUE
)


DefaultAssay(mnt) = 'RNA'
VariableFeatures(mnt)


saveRDS(mnt, file = paste0(outDir, '/SeuratObj_splice.unspliced_RA_day2_to_day6_all.rds'))

Idents(mnt) = mnt$condition
mnt = subset(mnt, downsample = 4000)
#saveDir = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/",
#                "results/scRNAseq_R13547_10x_mNT_20220813/RA_symetryBreaking/")

saveFile = '/RNAmatrix_umap_alevin.spliced_unspliced_RA_downsample_28k_rmCellCycle.correlatedGenes.h5Seurat'

SaveH5Seurat(mnt, filename = paste0(outDir, saveFile), 
             overwrite = TRUE)
Convert(paste0(outDir, saveFile), 
        dest = "h5ad", overwrite = TRUE)

