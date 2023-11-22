Dir_saveObjs = paste0('/groups/tanaka/People/current/jiwang/projects/RA_competence/',
                      'results/scRNAseq_R13547_10x_mNT_20220813/mapping_to_MouseGastrulationData/',
                      'dataMapping_subsettingRef_mNT.noRA.RA.d2_d5_SeuratRPCA/')

ref.combined = readRDS(paste0(Dir_saveObjs, 'integrated_mNT_mouseGastrulation_SeuratRPCA.rds'))

# make output directory if it doesn't exist
outDir = "/groups/tanaka/Collaborations/Jingkui-Hannah/RA_competence/scRNAseq_mNT/Hannahs_analysis2"
system(paste0('mkdir -p ', outDir))

cols_used = readRDS(paste0(Dir_saveObjs, 'integrated_mNT_mouseGastrulation_colorsUsed.rds'))

# Visualization
DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
        repel = TRUE, raster=FALSE, cols = cols_used)

ggsave(paste0(outDir, '/Integration_dataset_celltypes.pdf'), 
       width = 14, height = 8)

DimPlot(ref.combined, reduction = "umap", group.by = "dataset", raster=FALSE)
ggsave(paste0(outDir, '/Integration_dataset.pdf'), 
       width = 10, height = 8)

pdf(paste0(outDir, '/FeaturePlot_Markers.pdf'),
    width =10, height = 8, useDingbats = FALSE)


ggs = c('Pax6', 'Foxa2', 'Pou5f1', 'Sox17', 'Sox1', 'Sox2')
for(n in 1:length(ggs))
{
  p1 = FeaturePlot(ref.combined, features = ggs[n], min.cutoff = 'q5')
  #FeaturePlot(ref.combined, features = 'Foxa2', min.cutoff = 'q5')
  #FeaturePlot(ref.combined, features = 'Sox17', min.cutoff = 'q5')
  plot(p1)
}

dev.off()

DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE, split.by = 'dataset',
        repel = TRUE, raster=FALSE) + NoLegend()

ggsave(paste0(outDir, '/Integration_celltypes_split.dataset.pdf'), 
       width = 24, height = 8)

DimPlot(ref.combined, reduction = "umap", group.by = "stage", label = TRUE,
        repel = TRUE, raster=FALSE)

ggsave(paste0(outDir, '/Integration_stage.pdf'), 
       width = 16, height = 8)
