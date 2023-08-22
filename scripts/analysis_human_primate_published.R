##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: process the published human/primate embryo scRNA-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Feb 27 15:11:00 2023
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_human.primate.embryo.scRNAseq/'

resDir = paste0("../results/dataset_scRNAseq", version.analysis)
RdataDir = paste0(resDir, 'Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)


functionDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts'
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(Seurat)
library(DropletUtils)
library(edgeR)
library(future)
library(data.table)
library(tidyverse)
options(future.globals.maxSize = 80000 * 1024^2)
mem_used()

species = 'human.primate_scRNAseq'


########################################################
########################################################
# Section I : process Zhai et al. monkey data
# 
########################################################
########################################################
topdir = paste0('../published_human.primate.embryo/Zhai_monkey_CS8_11/filtered_feature_bc_matrix/')

exp = Matrix::readMM(paste0(topdir, "matrix.mtx.gz")) #read matrix

bc = read.csv(paste0(topdir, "/barcodes.tsv.gz"), header = F, stringsAsFactors = F)
g = read.csv(paste0(topdir, "/features.tsv.gz"), header = F, stringsAsFactors = F, sep = '\t')

genesymbol.mapping = openxlsx::read.xlsx(paste0('../published_human.primate.embryo/Tyser_et_al_2021/',
                                                '41586_2016_BFnature19096_MOESM297_ESM.xlsx'),
                                         sheet = 1, startRow = 2)

g$name = genesymbol.mapping$hg19_gene_symbol[match(g$V2, genesymbol.mapping$macFas5_gene_symbol)]

mm = which(is.na(g$name))
g$name[mm] = paste0(g$V2[mm], '.macFas5')

gg.counts = table(g$name)
gg.dup = names(gg.counts)[which(gg.counts>1)]

index.dup = which(!is.na(match(g$V2, gg.dup)))
g$name[index.dup] = paste0(g$V2[index.dup], '_', g$V1[index.dup])

colnames(exp) = bc$V1
rownames(exp) = g$name

count.data = exp
rm(exp);

meta = data.frame(row.names = colnames(count.data), dataset = rep('Zhai.et.al.', ncol(count.data)))
cellAnnot = read.csv('../published_human.primate.embryo/Zhai_monkey_CS8_11/MFE56636-meta.csv')
colnames(cellAnnot)[1] = 'cell.ids'
bcs = sapply(cellAnnot$cell.ids, function(x){unlist(strsplit(as.character(x), '_'))[2]})
cellAnnot$bc = bcs

meta$filtered = FALSE
mm = match(rownames(meta), cellAnnot$bc)
meta$filtered[which(is.na(mm))] = TRUE

meta = data.frame(meta, cellAnnot[mm, ], stringsAsFactors = FALSE)

# use defaultDrop to select cells.
hmp = CreateSeuratObject(counts = count.data,
                        meta.data = meta, 
                        min.cells = 10, min.features = 100)

hmp = subset(hmp, cells = colnames(hmp)[which(hmp$filtered == FALSE)])

embedding = as.matrix(cbind(hmp$UMAP_1, hmp$UMAP_2))
colnames(embedding) <- paste0("UMAP_", 1:2)
  
hmp[["umap"]] <- CreateDimReducObject(embeddings = embedding, key = "UMAP_", 
                                      assay = DefaultAssay(hmp))

p1 = DimPlot(hmp, group.by = 'stage', label = TRUE, repel = FALSE)
p2 = DimPlot(hmp, group.by = 'theiler_stage', label = TRUE, repel = FALSE)

p1 + p2 

ggsave(filename = paste0(resDir, '/Zhai_etal_stages.pdf'), width = 20, height = 8)

p1 = DimPlot(hmp, group.by = 'cell_cluster', label = TRUE, repel = FALSE)
p2 = DimPlot(hmp, group.by = 'cell_type', label = TRUE, repel = FALSE)

p1 + p2
ggsave(filename = paste0(resDir, '/Zhai_etal_celltypes.pdf'), width = 26, height = 8)

hmp$cell_type[which(hmp$cell_type == 'AI')] = 'Al'
hmp$cell_type[which(hmp$cell_type == 'Cardi. ')] = 'Cardi.'

celltypes = openxlsx::read.xlsx(paste0("../published_human.primate.embryo/Zhai_monkey_CS8_11/",
                                       "41586_2022_5526_MOESM3_ESM.xlsx"))
meta = hmp@meta.data
mm = match(hmp$cell_type, celltypes$Abbreviation)
cat(length(which(is.na(mm))), ' cell missing labels\n')

hmp$celltype_fullName = celltypes$Cell.types[mm]

p1 = DimPlot(hmp, group.by = 'cell_type', label = TRUE, repel = FALSE)
p2 = DimPlot(hmp, group.by = 'celltype_fullName', label = TRUE, repel = FALSE)

p1 + p2
ggsave(filename = paste0(resDir, '/Zhai_etal_celltypes_full.name.pdf'), width = 30, height = 8)

saveRDS(hmp, file = paste0(RdataDir, 'seuratObject_Zhai.et.al.2022_geneName.remapped.rds'))

##########################################
# double check some markers 
##########################################
hmp = readRDS(file = paste0(RdataDir, 'seuratObject_Zhai.et.al.2022_geneName.remapped.rds'))

FeaturePlot(hmp, features = c('FOXA2', 'POU5F1', 'OTX2', 'T'))

ggsave(filename = paste0(resDir, '/Zhai_etal_celltypes_featurePlots_FOXA2_OCT4.pdf'),
       width = 12, height = 10)

p1 = Seurat::VlnPlot(hmp, features = c('FOXA2'), group.by = 'cell_type')
p2 = Seurat::VlnPlot(hmp, features = c('POU5F1'), group.by = 'cell_type')

p1/p2

########################################################
########################################################
# Section II : import the data share by Maria
# 
########################################################
########################################################
dataDir = '../published_human.primate.embryo/Data_for_Hannah_fromMaria/Raw_data_integration/'

proteinCoding = read.delim(paste0(dataDir, 'Protein.coding.genes.filtering.txt'), sep = '\t', 
                           header = TRUE)

##########################################
# process first the full dataset from Tyser et al. (humanGastrula)  
##########################################
yy = readRDS(file = paste0(dataDir, 'raw_matrix_humanGastrula_Tyser2021.rds'))
yy1 = readRDS(file = paste0(dataDir, 'annot_umap_Tyser2021.rds'))
colnames(yy1)[1:3] = c('cell.index', 'umap_1', 'umap_2') 
yy = t(yy)
colnames(yy) = yy1$cell_name
mm = match(rownames(yy), proteinCoding$Feature)
yy = yy[which(!is.na(mm)), ]

meta = data.frame(row.names = yy1$cell_name, yy1)

# use defaultDrop to select cells.
tyser = CreateSeuratObject(counts = yy,
                         meta.data = meta, 
                         min.cells = 10, min.features = 100)

embedding = as.matrix(cbind(tyser$umap_1, tyser$umap_2))
colnames(embedding) <- paste0("UMAP_", 1:2)

tyser[["umap"]] <- CreateDimReducObject(embeddings = embedding, key = "UMAP_", 
                                      assay = DefaultAssay(tyser))

p1 = DimPlot(tyser, group.by = 'cluster_id', label = TRUE, repel = FALSE)
p2 = DimPlot(tyser, group.by = 'sub_cluster', label = TRUE, repel = FALSE)

p1 + p2 

ggsave(filename = paste0(resDir, '/Tyser_et_al_clusters_subclusters.pdf'), width = 20, height = 6)

tyser <- NormalizeData(tyser, normalization.method = "LogNormalize", scale.factor = 10000)
FeaturePlot(tyser, features = c('FOXA2', 'POU5F1', 'OTX2', 'T', 'NANOG', 'SOX11', 'SOX2', 'PAX6'))

ggsave(filename = paste0(resDir, '/Tyser.et.al_celltypes_featurePlots_FOXA2_OCT4.pdf'),
       width = 12, height = 10)

saveRDS(tyser, file = paste0(RdataDir, 'seuratObject_Tyser.2021_1195cells.rds'))

##########################################
# process Ma et al. 2019 data (monkey)
##########################################
meta = read.table(file = paste0(dataDir, 'Xiang.Ma.Tyser.integrated.dataset.annotation.txt'), header = TRUE)
genesymbol.mapping = openxlsx::read.xlsx(paste0('../published_human.primate.embryo/Tyser_et_al_2021/',
                                                '41586_2016_BFnature19096_MOESM297_ESM.xlsx'),
                                         sheet = 1, startRow = 2)

xx = fread(paste0(dataDir, "GSE130114_MF1453_Ma.csv.gz")) %>%
  as_tibble() %>%
  dplyr::rename(Gene = V1)

xx = as.data.frame(xx)

ggs = genesymbol.mapping$hg19_gene_symbol[match(xx$Gene, genesymbol.mapping$macFas5_gene_symbol)]
mm = which(is.na(ggs))
ggs[mm] = paste0(xx$Gene[mm], '.macFas5')

xx$Gene = ggs

#gg.counts = table(g$name)
#gg.dup = names(gg.counts)[which(gg.counts>1)]
#index.dup = which(!is.na(match(g$V2, gg.dup)))
#g$name[index.dup] = paste0(g$V2[index.dup], '_', g$V1[index.dup])

rownames(xx) = xx$Gene
xx = as.matrix(xx[, -1])

mm = match(colnames(xx), meta$cell.ID)
cat(sum(is.na(mm)), ' cell missing metadata \n')

meta = meta[mm, ]
rownames(meta) = meta$cell.ID

ma = CreateSeuratObject(counts = xx,
                           meta.data = meta, 
                           min.cells = 10, min.features = 100)


print(rownames(ma)[grep('FOXA|PAX|POU5F', rownames(ma))])

VlnPlot(ma, features = c("nCount_RNA", "nFeature_RNA"))

ma <- NormalizeData(ma, normalization.method = "LogNormalize", scale.factor = 10000)

ma <- FindVariableFeatures(ma, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(ma)

ma <- ScaleData(ma, features = all.genes)
ma <- RunPCA(ma, features = VariableFeatures(object = ma), verbose = FALSE)

ElbowPlot(ma, ndims = 30)

ma <- RunUMAP(ma, dims = 1:20, n.neighbors = 30, min.dist = 0.3)

DimPlot(ma, label = TRUE, repel = TRUE, group.by = 'day', raster=FALSE)
ggsave(filename = paste0(resDir, '/Ma_et_al.2019_umap_days.pdf'), width = 8, height = 6)

p1 = DimPlot(ma, label = TRUE, repel = TRUE, group.by = 'original.annotation', raster=FALSE)
p2 = DimPlot(ma, label = TRUE, repel = TRUE, group.by = 'new.identity', raster=FALSE)

p1 + p2

ggsave(filename = paste0(resDir, '/Ma_et_al.2019_umap_old_newAnnotation.pdf'), width = 16, height = 6)

FeaturePlot(ma, features = c('FOXA2', 'POU5F1', 'OTX2', 'T', 'NANOG', 'SOX1', 'SOX2', 'PAX6'))

ggsave(filename = paste0(resDir, '/Ma_et_al.2019_celltypes_featurePlots_FOXA2_OC4_markers.pdf'),
       width = 12, height = 10)

saveRDS(ma, file = paste0(RdataDir, 'seuratObject_Ma.2019_1453cells.rds'))


##########################################
# process Xiang et al. 2020 (Human)
##########################################
rm(genesymbol.mapping) # this is human data, no need for gene symbol mapping with monkey

meta = read.table(file = paste0(dataDir, 'Xiang.Ma.Tyser.integrated.dataset.annotation.txt'), header = TRUE)

xx = fread(paste0(dataDir, "merged_counts_Xiang.txt")) %>% 
  as_tibble() %>%
  dplyr::rename(Gene = gene)

xx = as.data.frame(xx)
ggs = proteinCoding$Feature[match(xx$Gene, proteinCoding$ID)] 
#ggs = genesymbol.mapping$hg19_gene_symbol[match(xx$Gene, genesymbol.mapping$macFas5_gene_symbol)]

xx$Gene = ggs
xx = xx[!is.na(xx$Gene), ]

rownames(xx) = make.unique(xx$Gene)
xx = as.matrix(xx[, -1])

cells = colnames(xx)
cells = sapply(cells, function(x) {test = unlist(strsplit(as.character(x), '_')); 
paste0(test[-c(1:3)], collapse = '_')})
colnames(xx) = cells

mm = match(colnames(xx), meta$cell.ID)
cat(sum(is.na(mm)), ' cell missing metadata \n')

meta = meta[mm, ]
rownames(meta) = meta$cell.ID

xx = CreateSeuratObject(counts = xx,
                        meta.data = meta, 
                        min.cells = 10, min.features = 100)

print(rownames(xx)[grep('FOXA|PAX|POU5F', rownames(xx))])

VlnPlot(xx, features = c("nCount_RNA", "nFeature_RNA"))

xx <- NormalizeData(xx, normalization.method = "LogNormalize", scale.factor = 10000)

xx <- FindVariableFeatures(xx, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(xx)

xx <- ScaleData(xx, features = all.genes)
xx <- RunPCA(xx, features = VariableFeatures(object = xx), verbose = FALSE)

ElbowPlot(xx, ndims = 30)

xx <- RunUMAP(xx, dims = 1:30, n.neighbors = 30, min.dist = 0.3)

DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'day', raster=FALSE)
ggsave(filename = paste0(resDir, '/Xiang.2020_umap_days.pdf'), width = 8, height = 6)

p1 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'original.annotation', raster=FALSE)
p2 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'new.identity', raster=FALSE)

p1 + p2

ggsave(filename = paste0(resDir, '/Xiang.2020_umap_old_newAnnotation.pdf'), width = 16, height = 6)

FeaturePlot(xx, features = c('FOXA2', 'POU5F1', 'OTX2', 'T', 'NANOG', 'SOX1', 'SOX2', 'PAX6'))

ggsave(filename = paste0(resDir, '/Xiang.2020_celltypes_featurePlots_FOXA2_OC4_markers.pdf'),
       width = 12, height = 10)


saveRDS(xx, file = paste0(RdataDir, 'seuratObject_Xiang.2020_555cells.rds'))

##########################################
# import the processed R object shared from Maria Rostovskya 
##########################################
emb = readRDS(file = paste0(dataDir, 'Xiang.Ma.Tyser.human.monkey.embryos.rds'))

emb <- NormalizeData(emb, normalization.method = "LogNormalize", scale.factor = 10000)

emb <- FindVariableFeatures(emb, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(emb)

emb <- ScaleData(emb, features = all.genes)
emb <- RunPCA(emb, features = VariableFeatures(object = emb), verbose = FALSE)

ElbowPlot(emb, ndims = 30)

emb <- RunUMAP(emb, dims = 1:30, n.neighbors = 50, min.dist = 0.3)

p1 = DimPlot(emb, label = TRUE, repel = TRUE, group.by = 'dataset', raster=FALSE)
p2 = DimPlot(emb, label = TRUE, repel = TRUE, group.by = 'new.identity', raster=FALSE)

p1
ggsave(filename = paste0(resDir, '/Xiang.Ma.Tyser.human.monkey.embryos_mergedMaria_dataset.pdf'), 
       width = 12, height = 8)

p2 
ggsave(filename = paste0(resDir, '/Xiang.Ma.Tyser.human.monkey.embryos_mergedMaria_celltypes_new.pdf'), 
       width = 12, height = 8)


