########################################################
########################################################
# cell type annotations
# 
########################################################
########################################################
system(paste0('mkdir -p ', resDir, '/manual_annotation/'))

outDir = paste0(resDir, '/manual_annotation/')


aa = readRDS(file = paste0(RdataDir, 
                'seuratObject_merged_cellFiltered_doubletFinderOut.v2_geneFiltered.15kGene_regressed.nCounts_', 
                           species, version.analysis, '.rds'))

aa$condition = factor(aa$condition, levels = levels)
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)

# remove the doublet cells
aa = subset(aa, DF_out == 'Singlet')

##########################################
# transfer annotated for from old version of analysis to the new one 
##########################################
Transfer.cell.labels.to.new.analysis = FALSE
if(Transfer.cell.labels.to.new.analysis){
  refs = readRDS(file = paste0(RdataDir, 
                               'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene_regressed.nCounts_clusterd0.7_', 
                               'manualAnnot_v1_', species, version.analysis, '.rds'))
  
  annots = refs@meta.data
  rm(refs)
  
  mm = match(colnames(aa), rownames(annots))
  length(which(is.na(mm)))
  length(which(!is.na(mm)))
  
  aa$celltypes = annots$celltypes[mm] 
  aa = subset(aa, cells = rownames(annots)) # discard blood cells
  
  # rerun the umap
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 50)
  
  Idents(aa) = aa$condition
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.2)
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', cols = cols, raster=FALSE)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes')
  
  saveRDS(aa, file = paste0(RdataDir, 
                'seuratObject_merged_cellFiltered_doublet.rm_mt.ribo.geneFiltered_regressout.nCounts_',
                'cellCycleScoring_annot.v1_', species, version.analysis, '.rds'))
  
}

##########################################
# rerun the PCA and umap 
# run clustering and save the cluster index
##########################################
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 100, min.dist = 0.2)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7) 

DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)
ggsave(filename = paste0(resDir, '/manual_annotation/UMAP_clusters.pdf'), width = 14, height = 10)


## save the first clusters
aa$clusters = aa$seurat_clusters
aa$celltypes = NA

saveRDS(aa, file = paste0(RdataDir, 
                     'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene_regressed.nCounts_clusterd0.7', 
                     species, version.analysis, '.rds'))

## discard the weird blood cells 
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene_regressed.nCounts_clusterd0.7', 
                           species, version.analysis, '.rds'))
aa$celltypes = NA
aa$subtypes = NA

DimPlot(aa, cells.highlight = colnames(aa)[which(aa$clusters == 17)])

cell.sels = colnames(aa)[which(aa$clusters == 17)]
mm = which(!is.na(match(colnames(aa), cell.sels)))
cat(length(cell.sels), '--', length(mm), '\n')
aa$celltypes[mm] = "Blood"

aa = subset(aa, cells = colnames(aa)[which(aa$clusters != 17)])

saveRDS(aa, file = paste0(RdataDir, 
          'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene_regressed.nCounts_clusterd0.7_rmBloodCells', 
                          species, version.analysis, '.rds'))

# reload the cleaned object
aa = readRDS(file = paste0(RdataDir, 
          'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene_regressed.nCounts_clusterd0.7_rmBloodCells', 
                           species, version.analysis, '.rds'))

aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 50)

aa <- RunUMAP(aa, dims = 1:50, n.neighbors = 50, min.dist = 0.1)

DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)


p1 = DimPlot(aa, label = TRUE, group.by = 'condition', cols = cols, repel = TRUE, raster=FALSE)
p2 = DimPlot(aa, label = TRUE, group.by = 'seurat_clusters', repel = TRUE, raster=FALSE)

p1 | p2

ggsave(filename = paste0(resDir, '/manual_annotation/UMAP_timepoints_clusters_v2.pdf'), 
       width = 20, height = 10)

# double check if clusters are due to sequencing depth and pert.mt
p1 = FeaturePlot(aa, features = 'nFeature_RNA')
p2 = FeaturePlot(aa, features = 'percent.mt')

p1 | p2

ggsave(filename = paste0(resDir, '/manual_annotation/FeaturesPlot_nFeatures_pertMT.pdf'), width = 14, height = 6)

##########################################
# check the cluster-specific markers
##########################################
markers = readRDS(file = paste0(RdataDir, 
        'seuratObject_merged_cellFiltered_doubletRemoved_geneFiltered.15kGene_regressed.nCounts_allMarkers_v1', 
                                species, version.analysis, '.rds'))
markers %>%
  filter(!str_detect(gene, '^(AMEX|LOC)')) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) -> top10

xx = subset(aa, downsample = 200)
DoHeatmap(xx, features = top10$gene) + NoLegend()
ggsave(filename = paste0(resDir, '/manual_annotation/MakrerGenes_20_perCluster_v2.pdf'), width = 45, height = 40)


FeaturePlot(aa, features = 'Sry', cols = c('gray', 'red'))

##########################################
##########################################
# Corase cluster annotation: ES, FP, Neuron, NP_RA and NP_noRA
##########################################
##########################################
Idents(aa) = aa$seurat_clusters

table(aa$seurat_clusters[which(aa$condition == 'day2_beforeRA')])

cluster_sels = c("7", '16', '4', '18')
sub.obj = subset(x = aa, idents = cluster_sels)

#sub.obj = SCTransform(aa, ncells = 3000, assay = "RNA", verbose = FALSE, 
#            variable.features.n = 3000, return.only.var.genes = TRUE, vst.flavor = "v2")
# sub.obj = NormalizeData(sub.obj, normalization.method = "LogNormalize", scale.factor = 10000)
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
sub.obj = ScaleData(sub.obj, vars.to.regress = 'nCount_RNA')

sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 30)

nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
n.neighbors = 20; min.dist = 0.1;
sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                   dims = 1:nb.pcs, 
                   n.neighbors = n.neighbors,
                   min.dist = min.dist)

sub.obj$clusters = sub.obj$seurat_clusters

#features = c('Pou5f1', 'Sox2', 'Lef1', 'Otx2', 'Zfp703', 'Pax6', 'Foxa2', 'Shh', 'Nkx6-1', 'Nkx2-2', 'Olig2', 
#             'Sox1', 'Tubb3', 'Bmp4', 'Bmp7', 'Nog', 'Pax3', 'Pax7', 'Arx')
features = c('Pou5f1', 'Sox2', 'Nanog', 'Zfp42', 'Klf4', 'Esrrb', 'Nr0b1', 'Dazl', 
             'Pou3f1', 'Otx2', 'Klf2', 'Klf5', 'Etv4', 'Etv5')

features = rownames(sub.obj)[!is.na(match(rownames(sub.obj), features))]
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))

VlnPlot(sub.obj, features = features)

#sub.obj <- FindNeighbors(sub.obj, dims = 1:20)
#sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.3)

p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) 

DimPlot(sub.obj, group.by = 'condition', reduction = 'umap', label = TRUE, label.size = 5) 

p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 6) +
  ggtitle(paste0(celltype.sels, ' -- subclusters'))

Idents(sub.obj) = sub.obj$clusters
markers = FindMarkers(sub.obj, ident.1 = 16, ident.2 = 7,
                      only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)

#save(markers, sub.obj, file = paste0(RdataDir, 'seuratObj_subset_day0.day1_markers.Rdata')) 
markers %>%
  #filter(!str_detect(gene, '^(AMEX|LOC|N/A)')) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

markers = FindMarkers(sub.obj, ident.1 = 18, ident.2 = 4,
                      only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)

xx = subset(sub.obj, downsample = 500)
#tops.common = c(markers.coarse$gene[!is.na(match(markers.coarse$cluster, c(1, 2, 3, 6, 13, 14, 15)))])

p3 = DoHeatmap(xx, features = top10$gene) + NoLegend()
p3
#ggsave(filename = paste0(resDir, '/first_test_clusterMarkers_v2.pdf'), width = 10, height = 30)

pdfname = paste0(resDir, '/subclustering_associatedMarkerGenes_', celltype.sels, '_v2.pdf')
pdf(pdfname, width=16, height = 16)

p1
p2
p3 

dev.off()

## save the subclustering labels 
cell.sels = colnames(aa)[which(aa$clusters == 7| aa$clusters == 16)]
mm = which(!is.na(match(colnames(aa), cell.sels)))
cat(length(cell.sels), '--', length(mm), '\n')
aa$celltypes[mm] = "ESC"

DimPlot(aa, group.by = 'celltypes')

##########################################
# annotate cells for RA day6 or cluster 8, 12, 13
##########################################
Idents(aa) = aa$seurat_clusters
aa$clusters = aa$seurat_clusters
#table(aa$seurat_clusters[which(aa$condition == 'day2_beforeRA')])

cluster_sels = c("8", '12', '13')
sub.obj = subset(x = aa, idents = cluster_sels)

#sub.obj = SCTransform(aa, ncells = 3000, assay = "RNA", verbose = FALSE, 
#            variable.features.n = 3000, return.only.var.genes = TRUE, vst.flavor = "v2")
# sub.obj = NormalizeData(sub.obj, normalization.method = "LogNormalize", scale.factor = 10000)
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
sub.obj = ScaleData(sub.obj, vars.to.regress = 'nCount_RNA')

sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 30)

nb.pcs = 30 # nb of pcs depends on the considered clusters or ids
n.neighbors = 30; min.dist = 0.1;

sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                   dims = 1:nb.pcs, 
                   n.neighbors = n.neighbors,
                   min.dist = min.dist)

sub.obj <- FindNeighbors(sub.obj, dims = 1:30)
sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.7)

pdfname = paste0(resDir, '/subclustering_annotation_', paste0(cluster_sels, collapse = "_"),
                 '_FloorPlate.pdf')
pdf(pdfname, width=16, height = 12)

DimPlot(sub.obj, group.by = 'condition', reduction = 'umap', label = TRUE, label.size = 5) 

p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) 
p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 6) +
  ggtitle(paste0(' subclusters'))

p1 | p2

#features = c('Pou5f1', 'Sox2', 'Lef1', 'Otx2', 'Zfp703', 'Pax6', 'Foxa2', 'Shh', 'Nkx6-1', 'Nkx2-2', 'Olig2', 
#             'Sox1', 'Tubb3', 'Bmp4', 'Bmp7', 'Nog', 'Pax3', 'Pax7', 'Arx')
#features = c('Pou5f1', 'Sox2', 'Nanog', 'Zfp42', 'Klf4', 'Esrrb', 'Nr0b1', 'Dazl', 
#             'Pou3f1', 'Otx2', 'Klf2', 'Klf5', 'Etv4', 'Etv5')

features = c('Foxa2', 'Shh', 'Arx', 'Pax6', 'Sox1', 'Nkx2-2', 'Nkx6-1', 'Olig2', 'Tubb3', 
             'Map2', 'Rbfox3', 'Irx3', 'Irx5', 'Sp8', 'Dbx2', 'Msx2', 'Lmx1b', 'Pax3', 'Pax7', 
             'Sox2', 'Elavl3', 'Sox10', 'Tlx2', 'Six1', 'Bmp4', 'Bmp7', 'Wnt1', 'Olig3')

features = c('Foxa2', 'Shh', 'Arx', 'Ferd3l', 'Lmx1b', 'Sox2', 'Nkx6-1', 
             'Tubb3', 'Elavl3') # floor plate markers
features = rownames(sub.obj)[!is.na(match(rownames(sub.obj), features))]
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))

VlnPlot(sub.obj, features = features)


dev.off()

## save the subclustering labels 
cell.sels = colnames(aa)[which(aa$clusters == 8)]
mm = which(!is.na(match(colnames(aa), cell.sels)))
cat(length(cell.sels), '--', length(mm), '\n')
aa$celltypes[mm] = "FP"

cell.sels = colnames(aa)[which(aa$clusters == 13)]
mm = which(!is.na(match(colnames(aa), cell.sels)))
cat(length(cell.sels), '--', length(mm), '\n')
aa$celltypes[mm] = "Neurons"

cell.sels = colnames(aa)[which(aa$clusters == 12)]
mm = which(!is.na(match(colnames(aa), cell.sels)))
cat(length(cell.sels), '--', length(mm), '\n')
aa$celltypes[mm] = "NP_RA"

cell.sels = colnames(aa)[which(aa$clusters == 11)]
mm = which(!is.na(match(colnames(aa), cell.sels)))
cat(length(cell.sels), '--', length(mm), '\n')
aa$celltypes[mm] = "NP_noRA"

DimPlot(aa, group.by = 'celltypes')

saveRDS(aa, file = paste0(RdataDir, 
'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene_regressed.nCounts_clusterd0.7_', 
'manualAnnot_v1_', species, version.analysis, '.rds'))


##########################################
# annotate all terminal cells (clusters 8, 11, 12, 13)
# too many cluster and also the cluster 8 is too different from others
##########################################
Idents(aa) = aa$seurat_clusters
aa$clusters = aa$seurat_clusters
#table(aa$seurat_clusters[which(aa$condition == 'day2_beforeRA')])

cluster_sels = c("8", '12', '13', '11')
sub.obj = subset(x = aa, idents = cluster_sels)

#sub.obj = SCTransform(aa, ncells = 3000, assay = "RNA", verbose = FALSE, 
#            variable.features.n = 3000, return.only.var.genes = TRUE, vst.flavor = "v2")
# sub.obj = NormalizeData(sub.obj, normalization.method = "LogNormalize", scale.factor = 10000)
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 2000)
sub.obj = ScaleData(sub.obj, vars.to.regress = 'nCount_RNA')

sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 30)

nb.pcs = 30 # nb of pcs depends on the considered clusters or ids
n.neighbors = 30; min.dist = 0.1;

sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                   dims = 1:nb.pcs, 
                   n.neighbors = n.neighbors,
                   min.dist = min.dist)

sub.obj <- FindNeighbors(sub.obj, dims = 1:30)
sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.7)

pdfname = paste0(resDir, '/subclustering_annotation_', paste0(cluster_sels, collapse = "_"),
                 '_FloorPlate.pdf')
pdf(pdfname, width=16, height = 12)

DimPlot(sub.obj, group.by = 'condition', reduction = 'umap', label = TRUE, label.size = 5) 

p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) 
p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 6) +
  ggtitle(paste0(' subclusters'))

p1 | p2

#features = c('Pou5f1', 'Sox2', 'Lef1', 'Otx2', 'Zfp703', 'Pax6', 'Foxa2', 'Shh', 'Nkx6-1', 'Nkx2-2', 'Olig2', 
#             'Sox1', 'Tubb3', 'Bmp4', 'Bmp7', 'Nog', 'Pax3', 'Pax7', 'Arx')
#features = c('Pou5f1', 'Sox2', 'Nanog', 'Zfp42', 'Klf4', 'Esrrb', 'Nr0b1', 'Dazl', 
#             'Pou3f1', 'Otx2', 'Klf2', 'Klf5', 'Etv4', 'Etv5')

features = c('Foxa2', 'Shh', 'Arx', 'Pax6', 'Sox1', 'Nkx2-2', 'Nkx6-1', 'Olig2', 'Tubb3', 
             'Map2', 'Rbfox3', 'Irx3', 'Irx5', 'Sp8', 'Dbx2', 'Msx2', 'Lmx1b', 'Pax3', 'Pax7', 
             'Sox2', 'Elavl3', 'Sox10', 'Tlx2', 'Six1', 'Bmp4', 'Bmp7', 'Wnt1', 'Olig3')

features = c('Foxa2', 'Shh', 'Arx', 'Ferd3l', 'Lmx1b', 'Sox2', 'Nkx6-1', 
             'Tubb3', 'Elavl3') # floor plate markers
features = rownames(sub.obj)[!is.na(match(rownames(sub.obj), features))]
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))

VlnPlot(sub.obj, features = features)



dev.off()

##########################################
##########################################
# Subclustering the corase clusters 
# continue to annotate NP and neurons with Hannah
# from clusters 11, 12, 13
##########################################
##########################################
#aa = readRDS(file = paste0(RdataDir, 
#                           'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene_regressed.nCounts_clusterd0.7_manualAnnot_v1_', 
#                           species, version.analysis, '.rds'))

cluster_sels = c('11', '12', '13')
sub.obj = subset(x = aa, idents = cluster_sels)

# sub.obj = NormalizeData(sub.obj, normalization.method = "LogNormalize", scale.factor = 10000)
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
sub.obj = ScaleData(sub.obj, vars.to.regress = 'nCount_RNA')

sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 30)

nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
n.neighbors = 30; min.dist = 0.1;

sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                   dims = 1:nb.pcs, 
                   n.neighbors = n.neighbors,
                   min.dist = min.dist)

sub.obj <- FindNeighbors(sub.obj, dims = 1:nb.pcs)
sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 1.0)

pdfname = paste0(resDir, '/subclustering_annotation_v2_moreHoxMarkers_', paste0(cluster_sels, collapse = "_"),
                 '_FloorPlate.pdf')
pdf(pdfname, width=16, height = 12)

p1 = DimPlot(sub.obj, group.by = 'condition', reduction = 'umap', label = TRUE, label.size = 5)
p2 = DimPlot(sub.obj, group.by = 'celltypes', reduction = 'umap', label = TRUE, label.size = 5) 
p1 | p2

ggsave(filename = paste0(outDir, 'subcluster_vhs_condition_celltype',  
                         paste0(cluster_sels, collapse = "_"), '.pdf'), width = 16, height = 8)

p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) 
p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 6) +
  ggtitle(paste0(' subclusters'))

p1 | p2

ggsave(filename = paste0(outDir, 'subcluster_vhs_',  paste0(cluster_sels, collapse = "_"), '.pdf'), width = 16, height = 8)

#features = c('Pou5f1', 'Sox2', 'Lef1', 'Otx2', 'Zfp703', 'Pax6', 'Foxa2', 'Shh', 'Nkx6-1', 'Nkx2-2', 'Olig2', 
#             'Sox1', 'Tubb3', 'Bmp4', 'Bmp7', 'Nog', 'Pax3', 'Pax7', 'Arx')
#features = c('Pou5f1', 'Sox2', 'Nanog', 'Zfp42', 'Klf4', 'Esrrb', 'Nr0b1', 'Dazl', 
#             'Pou3f1', 'Otx2', 'Klf2', 'Klf5', 'Etv4', 'Etv5')

features = c('Foxa2', 'Shh', 'Arx', 'Pax6', 'Sox1', 'Nkx2-2', 'Nkx6-1', 'Olig2', 'Tubb3', 
             'Map2', 'Rbfox3', 'Irx3', 'Irx5', 'Sp8', 'Dbx2', 'Msx2', 'Lmx1b', 'Pax3', 'Pax7', 
             'Sox2', 'Elavl3', 'Sox10', 'Tlx2', 'Six1', 'Bmp4', 'Bmp7', 'Wnt1', 'Olig3')

features = c('Foxa2', 'Shh', 'Arx', 'Ferd3l', 'Lmx1b', 'Sox2', 'Nkx6-1', 
             'Tubb3', 'Elavl3') # floor plate markers
features = rownames(sub.obj)[!is.na(match(rownames(sub.obj), features))]

FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))

features = unique(c('Sox2', 'Tubb3', 'Elavl3', 
                    'Hoxa2', 'Hoxa3', 'Hoxa4',
                    'Hoxb1','Hoxb2', 'Hoxb3','Hoxb4', 
                    'Hoxc4',    
                    'Hoxd3', 'Hoxd4', 'Hoxc6', 'Hoxc8', 'Phox2b', 'Hoxa9', 'Hoxc9',
                    'Cdx2', 'Hoxd13', 'Hoxb13', 'Hoxb9', 'Hoxb4', 'Hoxc4'))

FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))


features = c('Sox2', 'Tubb3', 'Elavl3')

features = c('Sox2', 'Foxa2', 'Shh', 'Arx', 'Pax6', 'Tubb3') # general markers
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))  

features = unique(c('Sox2', 'Sox1', 'Tubb3', 'Elavl3', 
                    'Irx3', 'Irx5', 'Pax3', 'Pax7',
              'Pax6', 'Olig2', 'Nkx2-9', 'Nkx2-2', 
              'Nkx6-1', 'Foxa2', 'Arx', 'Shh'
             )) # DV overview
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))  

ggsave(filename = paste0(outDir, 'subcluster_v_HS_DVaxis_markers_',  paste0(cluster_sels, collapse = "_"), '.pdf'), 
       width = 20, height = 16)

VlnPlot(sub.obj, features = features)


features = c('Sox2', 'Msx1', 'Msx2', 'Pax3', 'Wnt1', 'Olig3', 'Lmx1a'
            ) # RP
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))
ggsave(filename = paste0(outDir, 'subcluster_v_HS_RP_markers_',  paste0(cluster_sels, collapse = "_"), '.pdf'), 
       width = 10, height = 8)

features = c('Pax3', 'Pax7', 'Foxn4', 'Olig2', 'Nkx6-1', 'Nkx6-2')
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))

features = unique(c( 'Nkx2-2', 'Nkx6-1', 'Olig2',  
                     'Sox1', 'Msx1', 'Msx2', 'Pax3', 'Wnt1', 'Olig3',
                     'Irx3', 'Irx5', 'Pax6', 'Nkx6-1', 'Olig2', 'Nkx2-9', 'Nkx2-2',
                     'Map2', 'Rbfox3', 'Irx3', 'Irx5', 'Msx2', 'Lmx1b', 'Pax3', 'Pax7', 
                     'Sox2', 'Elavl3', 'Wnt1', 'Olig3', 'Nkx6-2',
                     'Ascl1', 'Gbx2', 'Dbx2', 'Dbx1','Sp8', 'Prdm12', 'Foxn4', 'Foxa2'
))

DotPlot(sub.obj, features = features, group.by = 'annot') + RotatedAxis()

VlnPlot(sub.obj, features = features)

VlnPlot(sub.obj, features = features)

# marker of all 22 clusters
markers = FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)

xx = subset(sub.obj, downsample = 500)
markers %>%
  #filter(!str_detect(gene, '^(AMEX|LOC|N/A)')) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

p3 = DoHeatmap(xx, features = top10$gene) + NoLegend()
p3                     
                      
ggsave(filename = paste0(outDir, '/subcluster_Markers_clusters.11.12,13.pdf'), width = 10, height = 30)

dev.off()

# ## save the subclustering labels 
# cell.sels = colnames(sub.obj)[which(!is.na(match(Idents(sub.obj), c(10, 11, 15, 16, 21))))]
# mm = which(!is.na(match(colnames(aa), cell.sels)))
# cat(length(cell.sels), '--', length(mm), '\n')
# aa$celltypes[mm] = "Neurons"
# 
# cell.sels = colnames(sub.obj)[which(!is.na(match(Idents(sub.obj), c(0, 3, 5, 8, 13, 14, 18))))]
# mm = which(!is.na(match(colnames(aa), cell.sels)))
# cat(length(cell.sels), '--', length(mm), '\n')
# aa$celltypes[mm] = "NP.RA"
# 
# cell.sels = colnames(sub.obj)[which(!is.na(match(Idents(sub.obj), c(1, 4, 7, 6, 12, 2, 17, 19, 9, 20))))]
# mm = which(!is.na(match(colnames(aa), cell.sels)))
# cat(length(cell.sels), '--', length(mm), '\n')
# aa$celltypes[mm] = "NP.noRA"

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene_', 
                          'regressed.nCounts_clusterd0.7_manualAnnot_v2_', 
                          species, version.analysis, '.rds'))




############################
## annotation only cluster 12: NP with RA
# v_HS
############################
#terminals = readRDS(file = paste0(RdataDir, 'annotation_ref.scmap.rds'))
#DimPlot(terminals, group.by = 'scmap.pred.id.50', reduction = 'umap', label = TRUE, label.size = 5) 
#sub.obj = terminals
# FeaturePlot(termianls, features = 'Hoxb4')

cell.sels = which(aa$clusters == 12)
sub.obj = subset(x = aa, cells = cell.sels)

sub.obj = NormalizeData(sub.obj, normalization.method = "LogNormalize", scale.factor = 10000)
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
sub.obj = ScaleData(sub.obj, vars.to.regress = 'nCount_RNA')

sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 30)

nb.pcs = 30 # nb of pcs depends on the considered clusters or ids
n.neighbors = 30; min.dist = 0.1;

sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                   dims = 1:nb.pcs, 
                   n.neighbors = n.neighbors,
                   min.dist = min.dist)

sub.obj <- FindNeighbors(sub.obj, dims = 1:nb.pcs)
sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 1.0)

DimPlot(sub.obj, group.by = 'celltypes', reduction = 'umap', label = TRUE, label.size = 5) 


DimPlot(sub.obj, group.by = 'condition', reduction = 'umap', label = TRUE, label.size = 5)
ggsave(filename = paste0(outDir, 'subcluster_Cluster12_vHS_condition.pdf'), width = 12, height = 8)

p0 = DimPlot(sub.obj, group.by = 'scmap.pred.id.30', reduction = 'umap', label = TRUE, label.size = 5) 

p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) 
p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 6) +
  ggtitle(paste0(' subclusters'))
p1 | p2
ggsave(filename = paste0(outDir, 'subcluster_Cluster12_vHS.pdf'), width = 20, height = 8)


p0|p2
ggsave(filename = paste0(outDir, 'subcluster_scmap_Neurons.pdf'), width = 20, height = 8)

pdfname = paste0(resDir, '/subclustering_annotation_', paste0(cluster_sels, collapse = "_"),
                 '_FloorPlate.pdf')
pdf(pdfname, width=16, height = 12)


features = unique(c('Sox2', 'Sox1', 'Tubb3', 'Elavl3', 
                    'Irx3', 'Irx5', 'Pax3', 'Pax7',
                    'Pax6', 'Olig2', 'Nkx2-9', 'Nkx2-2', 
                    'Nkx6-1', 'Foxa2', 'Arx', 'Shh'
)) # DV overview
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))  

ggsave(filename = paste0(outDir, 'subcluster_12_v_HS_DVaxis_markers_',  paste0(cluster_sels, collapse = "_"), '.pdf'), 
       width = 20, height = 16)

features = c('Tubb3', 'Elavl3', 'Isl1', 'Lhx3', 'Isl2', 'Mnx1', 'Slc10a4', 'Olig2', 
             'Aldn1a2', 'Arhgap36', 'Nkx2-2', 'Sim1', 'Pou4f1') # MN and V3

features = c('Tubb3', 'Elavl3', 'Pou4f1', 'Lhx2', 'Lhx9', 'Barhl1', 'Barhl2', 'Atoh1', 
             'Foxd3', 'Lhx1', 'Lhx5', 'Nkx2-2', 
             'Lhx3', 'Vsx2', 'Sox21', 'Foxn4', 'Vsx1', 'Tal1', 'Gata2', 'Gata3', 'En1', 'Evx1', 'Evx2') # dl1 and dl2

FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))

features = unique(c('Sox2', 'Tubb3', 'Elavl3', ## NP with and without RA 
                    'Hoxa2', 'Hoxa3', 'Hoxa4',
                    'Hoxb1','Hoxb2', 'Hoxb3','Hoxb4', 
                    'Hoxc4',    
                    'Hoxd3', 'Hoxd4'))


features = c('Sox2', 'Foxa2', 'Shh', 'Arx', 'Pax6', 'Tubb3') # general markers

features = c('Pax6', 'Sox1', 'Foxa2',
             'Irx3', 'Irx5', 'Pax6', 'Pax3', 'Sp8', 'Nkx6-1', 'Olig2', 'Nkx2-9', 'Nkx2-2') # p3, pMN
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))  

VlnPlot(sub.obj, features = features)


features = c('Sox1', 'Msx1', 'Msx2', 'Pax3', 'Wnt1', 'Olig3',
             'Irx3', 'Irx5', 'Pax6', 'Sp8', 'Nkx6-1', 'Olig2', 'Nkx2-9', 'Nkx2-2') # RP
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))
features = c('Pax3', 'Pax7', 'Foxn4', 'Olig2', 'Nkx6-1', 'Nkx6-2')
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))

features = unique(c( 'Nkx2-2', 'Nkx6-1', 'Olig2',  
                     'Sox1', 'Msx1', 'Msx2', 'Pax3', 'Wnt1', 'Olig3',
                     'Irx3', 'Irx5', 'Pax6', 'Nkx6-1', 'Olig2', 'Nkx2-9', 'Nkx2-2',
                     'Map2', 'Rbfox3', 'Irx3', 'Irx5', 'Msx2', 'Lmx1b', 'Pax3', 'Pax7', 
                     'Sox2', 'Elavl3', 'Wnt1', 'Olig3', 'Nkx6-2',
                     'Ascl1', 'Gbx2', 'Dbx2', 'Dbx1','Sp8', 'Prdm12', 'Foxn4', 'Foxa2'
))

DotPlot(sub.obj, features = features, group.by = 'annot') + RotatedAxis()


VlnPlot(sub.obj, features = features)

dev.off()

## save the subclustering labels 
#aa$subtypes = aa$celltypes
sub.obj$annot = as.character(sub.obj$RNA_snn_res.1)

jj = which(!is.na(match(Idents(sub.obj), c(9))))
subtypes = 'MN' # Isl1, Olig2 and Mnx1
sub.obj$annot[jj] = subtypes
mm = which(!is.na(match(colnames(aa), colnames(sub.obj)[jj])))
cat(length(jj), '--', length(mm), '\n')
aa$subtypes[mm] = subtypes

jj = which(!is.na(match(Idents(sub.obj), c(8)))) # Nkx2-2, Nkx6-1, no Pax6
subtypes = 'pMN' 
sub.obj$annot[jj] = subtypes
mm = which(!is.na(match(colnames(aa), colnames(sub.obj)[jj])))
cat(length(jj), '--', length(mm), '\n')
aa$subtypes[mm] = subtypes

jj = which(!is.na(match(Idents(sub.obj), c(12)))) # Nkx2-2, Nkx6-1, no Pax6
subtypes = 'RP' 
sub.obj$annot[jj] = subtypes
mm = which(!is.na(match(colnames(aa), colnames(sub.obj)[jj])))
cat(length(jj), '--', length(mm), '\n')
aa$subtypes[mm] = subtypes


saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene_',
                          'regressed.nCounts_clusterd0.7_manualAnnot_v2.7_', 
                          species, version.analysis, '.rds'))

############################
## annotation only cluster 8 (FP) to check the DV axis stuff
# v_HS
############################
cell.sels = which(aa$clusters == 8)
sub.obj = subset(x = aa, cells = cell.sels)

sub.obj = NormalizeData(sub.obj, normalization.method = "LogNormalize", scale.factor = 10000)
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
sub.obj = ScaleData(sub.obj, vars.to.regress = 'nCount_RNA')

sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 30)

nb.pcs = 30 # nb of pcs depends on the considered clusters or ids
n.neighbors = 30; min.dist = 0.1;

sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                   dims = 1:nb.pcs, 
                   n.neighbors = n.neighbors,
                   min.dist = min.dist)

sub.obj <- FindNeighbors(sub.obj, dims = 1:nb.pcs)
sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 1.0)

DimPlot(sub.obj, group.by = 'celltypes', reduction = 'umap', label = TRUE, label.size = 5) 

DimPlot(sub.obj, group.by = 'condition', reduction = 'umap', label = TRUE, label.size = 5)
ggsave(filename = paste0(outDir, 'subcluster_Cluster8.FP_vHS_condition.pdf'), width = 12, height = 8)


p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) 
p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 6) +
  ggtitle(paste0(' subclusters'))
p1 | p2
ggsave(filename = paste0(outDir, 'subcluster_Cluster8.FP_vHS.pdf'), width = 20, height = 8)


p0|p2
ggsave(filename = paste0(outDir, 'subcluster_scmap_Neurons.pdf'), width = 20, height = 8)

pdfname = paste0(resDir, '/subclustering_annotation_', paste0(cluster_sels, collapse = "_"),
                 '_FloorPlate.pdf')
pdf(pdfname, width=16, height = 12)


features = unique(c('Sox2', 'Sox1', 'Tubb3', 'Elavl3', 
                    'Irx3', 'Irx5', 'Pax3', 'Pax7',
                    'Pax6', 'Olig2', 'Nkx2-9', 'Nkx2-2', 
                    'Nkx6-1', 'Foxa2', 'Arx', 'Shh'
)) # DV overview
FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))  

ggsave(filename = paste0(outDir, 'subcluster_8FP_v_HS_DVaxis_markers.pdf'), 
       width = 20, height = 16)


DotPlot(sub.obj, features = features, group.by = 'annot') + RotatedAxis()

VlnPlot(sub.obj, features = features)

dev.off()


saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObject_merged_cellFiltered_doubletRm_geneFiltered.15kGene_',
                          'regressed.nCounts_clusterd0.7_manualAnnot_v2.7_', 
                          species, version.analysis, '.rds'))







