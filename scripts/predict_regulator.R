##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: predict Foxa2 and Pax6 regulators
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Aug 23 16:11:31 2023
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

outDir = paste0(resDir, '/RA_symetryBreaking/regulators_prediction')
system(paste0('mkdir -p ', outDir))

# tfs and sps annotations 
sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
sps = unique(sps$gene)

tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)
tfs = as.character(unlist(sapply(tfs, firstup)))

source('functions_utility.R')

##########################################
# import the data
##########################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
                           'cellCycleScoring_annot.v2_newUMAP_clusters_sparseFeatures', data_version, '_',
                           species, version.analysis, '.rds'))

# bb  = aa
p1 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
p2 = DimPlot(aa, group.by = 'clusters', label = TRUE, repel = TRUE)

p1 + p2


########################################################
########################################################
# Section I : make some guess
# 
########################################################
########################################################
FeaturePlot(aa, features = c('Pax6', 'Foxa2', "Tet2", "Tet3", "Tet1",
                             "Dnmt3b", "Dnmt1", "Dnmt3a", "Crabp2"))

ggsave(filename = paste0(outDir, '/FeaturePlots_Foxa2_pax6_regulatorGuesss.pdf'), 
       width = 14, height = 12)

FeaturePlot(aa, features = c('Pax6', 'Foxa2', 'Stra8', 'Halr1', 'Phc1', 'Aldh1a2'))


ggsave(filename = paste0(outDir, '/FeaturePlots_Foxa2_pax6_RAresponseGenesGuesss.pdf'), 
       width = 14, height = 12)

FeaturePlot(aa, features = c('Foxa2', 'Sox17'))

FeaturePlot(aa, features = c('Foxa2', 'Lhx1', 'Otx2', 'Eomes', 'Ldb1'))
FeaturePlot(aa, features = c('Foxa2', 'Smad2', 'Smad4', 'Amot', 'Dkk1'))

FeaturePlot(aa, features = c('Foxa2', 'Yap1', 'Tead1', 'Tead2', 'Tead3', 'Tead4'))

ggsave(filename = paste0(outDir, '/FeaturePlots_Foxa2_pax6_TeadGenes.pdf'), 
       width = 12, height = 8)


##########################################
# check the Yap targets
##########################################
FeaturePlot(aa, features = c('Foxa2', 'Ajuba', 'Limd1', 'Dlg5'))

targets = read.csv(file = paste0('/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/',
                                 'scRNAseq_mNT/RA_symetryBreaking/decoupleR_Smad_Tead_targets/', 
                                 'targets_decoupleR_for_Smad_Tead_activtiyInference.csv'))

targets = targets[grep('Tead', targets$source), ]

pdf(paste0(outDir, '/TransientGenes_Tead_targets.pdf'),
    width =16, height = 10, useDingbats = FALSE)

plot_manyFeatures_seurat(seurat_obj = aa, features = unique(c('Foxa2', targets$target)))

dev.off()


FeaturePlot(aa, features = c('Foxa2', 'Ctgf', 'Cyr61', 'Ax6', 'Birc5', 'Ankrd1', 'Amot',
                             'Amotl1', 'Amotl2', 'Axin2', 'Sox2', 'Id1', 'Id2'))

FeaturePlot(aa, features = c('Foxa2', 'Pax6', 'Etv6'))

FeaturePlot(aa, features = c('Foxa2', 'Pax6', 'Cdx1'))
FeaturePlot(aa, features = c('Foxa2', 'Pax6', "Mafb",  "Cflar",  "Bbx", 
                             "Prdm1",  "Prrx2",  "Sox17",  "Lhx1",  "Pdcd4",  
                             "Atp6ap2",  "Fgf3",  "Tmf1",  "Fzd2"))

##########################################
# transiently expressed gene at day2.5 and day3
##########################################
Idents(aa) = aa$condition
all.markers <- FindMarkers(aa, ident.1 = 'day2.5_RA', 
                           ident.2 = c('day2_beforeRA', 'day3.5_RA', 'day4_RA', 'day5_RA'),
                           only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)

markers2 <- FindMarkers(aa, ident.1 = 'day3_RA.rep1', 
                           ident.2 = c('day2_beforeRA', 'day3.5_RA', 'day4_RA', 'day5_RA'),
                           only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)

all.markers =rbind(all.markers, markers2[which(is.na(match(rownames(markers2), rownames(all.markers)))), ])

all.markers = all.markers[order(-all.markers$avg_log2FC), ]

saveRDS(all.markers, file = paste0(outDir, '/transiently_expressed_genes_day2.5_andDay3.rds'))

## import the RA and noRA samples
bb =  readRDS(file = paste0(RdataDir, 
                            'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                            'cellCycleScoring_annot.v1_', species, version.analysis, '.rds'))

levels_sels = c("day2_beforeRA",  
                "day2.5_RA", "day3_RA.rep1", "day3.5_RA",   "day4_RA", "day5_RA",
                "day2.5_noRA", "day3_noRA",  "day3.5_noRA", "day4_noRA", "day5_noRA")

Idents(bb) = factor(bb$condition, levels = levels)

bb = subset(bb, idents = levels_sels)

bb <- FindVariableFeatures(bb, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs

## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
bb <- RunPCA(bb, features = VariableFeatures(object = bb), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(bb, ndims = 50)

Idents(bb) = bb$condition
bb <- RunUMAP(bb, dims = 1:30, n.neighbors = 100, min.dist = 0.2)
DimPlot(bb, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

p1 = DimPlot(bb, group.by = 'condition', label = TRUE, repel = TRUE) + NoLegend()
p2 = FeaturePlot(bb, features = c('Foxa2', 'Pax6'), ncol = 1)

p1 + p2

ggsave(filename = paste0(outDir, '/UMAP_RA_noRA_samples.pdf'), 
       width = 16, height = 8)

pdf(paste0(outDir, '/TransientGenes_day2.5_day3_top124_onlyTFs.SPs_umap.RAnoRA.pdf'),
    width =16, height = 10, useDingbats = FALSE)

all.markers = readRDS(file = paste0(outDir, '/transiently_expressed_genes_day2.5_andDay3.rds'))

xx = all.markers[which(!is.na(match(rownames(all.markers), c(tfs, sps)))), ]

plot_manyFeatures_seurat(seurat_obj = bb, features = rownames(xx))

dev.off()


########################################################
########################################################
# Section II : scan Foxa2 enhancer and promoter
# 
########################################################
########################################################
Dir_mouse_annot = "/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/"

promoters = read.delim(paste0(Dir_mouse_annot, 'ens_BioMart_GRCm38.p6.txt'))
promoters = promoters[which(promoters$Gene.type == 'protein_coding'), ]
promoters = promoters[which(promoters$GENCODE.basic.annotation == 'GENCODE basic'), ]

promoters = promoters[which(!is.na(match(promoters$Chromosome.scaffold.name, c(1:19)))), ]

# select only one promoter for each gene symbols
Make_background = FALSE

if(Make_background){
  ggs = unique(promoters$Gene.name)
  promoters = promoters[match(ggs, promoters$Gene.name), ]
  
  # randomly select 2000 promoters as background
  set.seed(2023)
  sels = sample(c(1:nrow(promoters)), size = 2000)
}else{
  sels = which(promoters$Gene.name == 'Foxa2'| promoters$Gene.name == 'Pax6')
}

promoters = promoters[sels, ]

bgs = data.frame(promoters$Chromosome.scaffold.name, promoters$Transcription.start.site..TSS.,
                 promoters$Transcription.start.site..TSS., promoters$Transcript.name, rep(0, nrow(promoters)), 
                 promoters$Strand, stringsAsFactors = FALSE)
colnames(bgs) = c('chr', 'start', 'end', 'gene', 'score', 'strand')
bgs$strand[which(bgs$strand == '1')] = '+'
bgs$strand[which(bgs$strand == '-1')] = '-'

jj = which(bgs$strand == '+')
bgs$start[jj] = bgs$end[jj] - 2000

jj = which(bgs$strand == '-')
bgs$end[jj] = bgs$start[jj] + 2000

write.table(bgs, file = paste0(Dir_mouse_annot, 'mm10_Foxa2_Pax6_promoters.bed'), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')


##########################################
# double check the potential TFs after FoxA2 promoter and enhancer scanning
##########################################
# test AME enriched motifs
pdf(paste0(outDir, '/AME_enrichedMotifs.pdf'),
    width =16, height = 10, useDingbats = FALSE)

features = c('Foxa2', 'Pax6', 
             rownames(aa)[grep('Onecut', rownames(aa))], 
             rownames(aa)[grep('Zfp384', rownames(aa))],
             'RREB1',
             rownames(aa)[grep('Egr', rownames(aa))],
             rownames(aa)[grep('Prdm', rownames(aa))],
             rownames(aa)[grep('Fox', rownames(aa))],
             rownames(aa)[grep('Pou', rownames(aa))],
             rownames(aa)[grep('Klf', rownames(aa))]
)
plot_manyFeatures_seurat(seurat_obj = aa, features = unique(features))

dev.off()

FeaturePlot(aa, features = c('Foxa2', 'Pax6', 'Klf4', 'Rreb1', 'Mzf1', 
                             'Zfp384','Zfp148', 'Zfp281','Patz1'
                             #rownames(aa)[grep('Zic', rownames(aa))],
                             #rownames(aa)[grep('Cdx', rownames(aa))]
                             ))

FeaturePlot(aa, features = c('Foxa2', 'Pax6', 
                             rownames(aa)[grep('Sp', rownames(aa))]
                             #rownames(aa)[grep('Cdx', rownames(aa))]
))

FeaturePlot(aa, features = c('Foxa2', 'Pax6', 
                             rownames(aa)[grep('Bptf', rownames(aa))]
                             #rownames(aa)[grep('Cdx', rownames(aa))]
                             ),
            max.cutoff = 'q99')

FeaturePlot(aa, features = c('Foxa2', 'Pax6', 
                             rownames(aa)[grep('^Rar|^Rxr', rownames(aa))]
                             #rownames(aa)[grep('Cdx', rownames(aa))]
))

FeaturePlot(aa, features = c('Foxa2', 'Pax6', 
                             rownames(aa)[grep('^Stat', rownames(aa))]
                             #rownames(aa)[grep('Cdx', rownames(aa))]
))

FeaturePlot(aa, features = c('Foxa2', 'Pax6', 
                             'Otx1', 'Lhx1'
))


##########################################
# make summary of motif occurrency in Foxa2 promoters and enhancers 
##########################################
library(data.table)

fimo.out = paste0('../results/scRNAseq_R13547_10x_mNT_20220813/motif_analysis/fimo/',
                  'FoxA2_promoters_enhancers/', 'fimo.tsv')

fimo = fread(fimo.out, header = TRUE)
fimo = fimo[which(fimo$`p-value` < 10^-4), ]

motif.oc = table(fimo$motif_id, fimo$sequence_name, useNA = 'ifany')
colnames(motif.oc) = c('enhancer.fp', 'enhancer.node', 'promoter_201', 'promoter_202')
motif.oc[grep('SMAD|TEAD|LEF|^TCF|FOXA|RAR|GLI|PAX', rownames(motif.oc)), ]

motif.oc = t(motif.oc)

print(head(rownames(motif.oc), 20))

saveRDS(motif.oc, file = paste0(outDir, '/motif_oc_fimo_jaspar2022_pval.0.0001_Foxa2.enhancers.promoters.rds'))

ss = apply(motif.oc, 2, sum)
motif.oc = motif.oc[, order(-ss)]

########################################################
########################################################
# Section III : Test scran correlaton analysis
# original code from 
# https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/misc.html
########################################################
########################################################
library(R.utils)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(ggrepel)

motif.oc = readRDS(file = paste0(outDir, '/motif_oc_fimo_jaspar2022_pval.0.0001_Foxa2.enhancers.promoters.rds'))
motif.oc = t(motif.oc)
mapping = readRDS(paste0('../results/scRNAseq_R13547_10x_mNT_20220813/motif_analysis/', 
                         'JASPAR2022_CORE_UNVALIDED_vertebrates_nonRedundant_metadata_manual_rmRedundantUNVALIDED.rds'))

mapping = mapping[which(!is.na(match(mapping$name, rownames(motif.oc)))), ]
mapping$gene = firstup(mapping$gene)

Idents(aa) = factor(aa$condition)


levels_sels = c("day3_RA.rep1")

data_version = "_d3"

outDir_cc = paste0(outDir, '/scran_pairCorrelation_scran', data_version)
system(paste0('mkdir -p ', outDir_cc))

sce = subset(aa, idents = levels_sels)
sce = as.SingleCellExperiment(sce)

set.seed(100)
mm = match(rownames(sce), c(tfs, sps))
subset = which(!is.na(mm))
var.cor <- correlatePairs(sce, pairings = list(rep(which(rownames(sce) == 'Foxa2'), length(subset)), 
                                               subset))

head(var.cor)

saveRDS(var.cor, file = paste0(outDir_cc, '/correlatePairs_FoxA2_vs_tfs.sps.rds'))

var.cor = readRDS(file = paste0(outDir_cc, '/correlatePairs_FoxA2_vs_tfs.sps.rds'))

sig.cor <- var.cor$FDR <= 0.01
summary(sig.cor)

var.cor = var.cor[which(sig.cor == TRUE), ]

ave.counts <- calculateAverage(sce)

var.cor$findMotifs = !is.na(match(var.cor$gene2, mapping$gene))
var.cor$fdr = -log10(var.cor$FDR+10^-10)
var.cor$ave.counts = log10(ave.counts[match(var.cor$gene2, names(ave.counts))])

var.cor = as.data.frame(var.cor)

ggplot(data = var.cor, aes(x=rho, y=ave.counts, label = gene2)) +
  geom_point(size = 0.2) + 
  theme_classic() + 
  labs(x = "Correlation", y = "log10(ave.counts)") +
  geom_text_repel(data= var.cor[which(var.cor$findMotifs == TRUE & abs(var.cor$rho) >0.1), ], 
                  aes(x=rho, y=ave.counts),
                  size = 4,
                  #color = "darkgreen",
                  #family = 'Times',
                  #fontface = 'bold',
                  # Add extra padding around each text label.
                  #box.padding = unit(0.3, 'lines'),
                  # Add extra padding around each data point.
                  #point.padding = unit(1.6, 'lines')
                  )
ggsave(filename = paste0(outDir_cc, '/scran_correlation.with.FoxA2.pdf'), width = 20, height = 16)

dev.off()

pdf(paste0(outDir_cc, '/Compare_expressionPattern.pdf'),
    width =16, height = 10, useDingbats = FALSE)

features = c('Foxa2', 'Pax6', 
             var.cor$gene2[which(var.cor$findMotifs == TRUE & abs(var.cor$rho) >0.1)]
)
plot_manyFeatures_seurat(seurat_obj = aa, features = unique(features))

dev.off()

########################################################
########################################################
# Section IV : Genie3 random-forest-based method to predict TF regulators for FoxA2
# 
########################################################
########################################################
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

features = c('Foxa2', 'Pax6', wtm$regulator[which(wtm$Foxa2>10^-)])

plot_manyFeatures_seurat(seurat_obj = aa, features = unique(features))

dev.off()







