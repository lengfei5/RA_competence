inputDir = '../images_data/test_WT/'

xx = read.csv(file = paste0(inputDir, 'cyst_size_genotype_FoxA2_Pax6.csv'), row.names = c(1))

xx = xx[which(xx$pct_genotype > 0.3 & xx$pct_genotype < 0.9 & xx$pct_wt > 0.05 &   xx$cyst_size > 10^4), ]

xx$either = xx$pct_foxa2 + xx$pct_pax6
plot(xx$pct_wt, xx$pct_foxa2)

plot(xx$pct_wt, xx$either)


p = seq(0, 1.0, by = 0.1)

for(i in 1:(length(p)-1))
{
  ii = which(xx$pct_wt >= p[i] & xx$pct_wt < p[i+1])
  cat(p[i] + 0.05, '--', median(xx$pct_foxa2[ii]), '\n')
}


library(ggplot2)
xx = data.frame(xx, stringsAsFactors = FALSE)

ggplot(xx, aes(x=pct_wt, y=pct_foxa2)) + 
  geom_point()+
  geom_smooth(method=lm)
# Remove the confidence interval


require(graphics)

utils::example(nhtemp)
myNHT <- as.vector(nhtemp)
myNHT[20] <- 2 * nhtemp[20]
plot(myNHT, type = "b", ylim = c(48, 60), main = "Running Medians Example")
lines(runmed(myNHT, 7), col = "red")

jj = order(xx$pct_wt)
xx = xx[jj, ]
plot(xx$pct_wt, xx$pct_foxa2)
lines(xx$pct_wt, runmed(xx$pct_foxa2, 5), col = "red")



########################################################
########################################################
# Section : test SH output to summarize with UMAP
# 
########################################################
########################################################
library(CATALYST)
library(dplyr)
library(flowCore)
library(flowWorkspace)
library(ggcyto)
library(ggplot2)
library(mvtnorm)
library(openCyto)
library(scDataviz)
library(cowplot)
library(SingleCellExperiment)
library(scater)

figureDir = "../results/figures_tables_R13547_10x_mNT_20240522/"

#sce = readRDS(file = '../data/image_SHout_sce.rds')
#sce = readRDS(file = '../data/image_SHout_l_sce_wt_conditions_metadata.rds')
#sce = readRDS(file = '../data/image_SHout_dlogl_sce_wt_conditions_metadata.rds')
sce = readRDS(file = "../data/image_SHout_dlogl_sce_wt_conditions_metadata_filteringFoxA2SD_v3.rds")

#xx = readRDS(file = '../data/image_SHout_dlogl_collectedFeatures_sce_wt_conditions_metadata.rds')
#xx = xx[, !is.na(match(colnames(xx), colnames(sce)))]

#sce = xx

sce$size_log = log10(sce$cyst_size)
sce$r2_log = log10(sce$cyst_r^2)

sce$genotype_foxa2[which(sce$condition == 'WT')] = 0.3

nb_l = 20

#sce = sce[1:20, ]

#sce = sce[1:(nb_l +1), which(sce$time != 'd3.5' & sce$time != 'd5')]
#sce = sce[2:(nb_l + 1), which(sce$condition == "WT")]
#sce = sce[2:(nb_l + 1), ]
sce = sce[2:(nb_l + 1), ]

#sce = sce[, which(sce$time != 'd3.5' & sce$time != 'd5')]
#sce = sce[2:(nb_l + 1), which(sce$condition != "WT" & sce$condition != "KO_KO")]
sce = sce[ ,which(sce$time != 'd3.5' & sce$time != 'd5')]
sce = sce[ ,which(colnames(sce) != "211210_d6_RAd2_D2_49_01_isotropic_cyst_7")]

sce$condition = factor(sce$condition, levels = c("WT", "KO_KO", "TetOn_TetON_RA", "TetOn_TetON_dox"))

table(sce$condition, sce$time)

##########################################
# data transformation: tricky
##########################################
y <- assay(sce, "counts")

#y = y[which(rownames(y) != "label_foxa2"), ]

print(range(log10(y)))

y = log10(y) + 16
#y = log10(y )
#for(n in 1:nrow(y)) y[n, ] = y[n, ]/n 
#y <- asinh(sweep(y, 1, cf, "/"))

assay(sce, "exprs", FALSE) <- y

# run t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(1234)
#sce <- runDR(sce, "TSNE", cells = 500, features = "type")
sce <- runDR(sce, "PCA", ncomponents = 5,  scale = FALSE, assay = "exprs",
             cells = NULL, features = "type")


plotDR(sce, "PCA", dims = c(1, 2), color_by = "time", facet_by = "condition", ncol = 2
       #k_pal = c("lightgrey", "cornflowerblue", "navy")
       ) + coord_fixed(2) + 
  theme_classic() + 
  geom_point(size=0.8) +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size=14)) 

ggsave(paste0(figureDir, 'image_WT_vsConditions_RA_SHout_PCA_time_v2.pdf'), width=10, height = 6) 


plotDR(sce, "PCA", dims = c(1, 2), color_by = "r2_log", facet_by = "condition", ncol = 2
       #k_pal = c("lightgrey", "cornflowerblue", "navy")
) + coord_fixed(2) + 
  theme_classic() + 
  geom_point(size=0.8) +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12), 
        legend.text = element_text(size=12), 
        legend.title = element_text(size=12)) 

ggsave(paste0(figureDir, 'image_WT_vsConditions_RA_SHout_PCA_surfaceSize_v2.pdf'), width=10, height = 6) 


plotDR(sce, "PCA", dims = c(1, 2), color_by = "genotype_foxa2", facet_by = "condition", ncol = 2
       #k_pal = c("lightgrey", "cornflowerblue", "navy")
) + coord_fixed(2) + 
  geom_point(size=0.8) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12), 
        legend.text = element_text(size=12), 
        legend.title = element_text(size=12)) +
  ggplot2::scale_colour_gradient2(limits = c(0, 1), breaks = c(0, 0.1, 0.2,  0.3, 0.5, 0.7, 0.9),
                                  low = "#313695", mid = "#FFFFBF", high = "#C61010", midpoint = 0.3)

ggsave(paste0(figureDir, 'image_WT_vsConditions_RA_SHout_PCA_genotype_v2.pdf'), width=10, height = 6) 


p1 = plotDR(sce, "PCA", color_by = "condition")

p2 = plotDR(sce, "PCA", color_by = "condition", dims = c(3, 4))

p1 + p2

### filter outliers cells
pcs = reducedDim(sce, 'PCA')

jj = which(pcs[, 1] > 10 & sce$condition == "TetOn_TetON_dox" & sce$time == 'd3')

#jj = which(pcs[, 1] < 10 & sce$condition == "WT" & sce$time == 'd3')
#sels = which(pcs[ ,1] > (-10) & abs(pcs[, 2]) < (2))
#sce = sce[, sels]


# set.seed(1234)
# sce <- runDR(sce, "UMAP", cells = NULL, 
#              features = "type",
#              n_neighbors = 20, scale = FALSE,
#              min_dist = 0.05, metric = "cosine"
# )
# 
# p1 = plotDR(sce, "UMAP", color_by = "condition")
# p1
# 

set.seed(1234)
sce <- runDR(sce, "UMAP", #pca = 5,
             cells = NULL, 
             features = "type",
             n_neighbors = 20, scale = FALSE,
             min_dist = 0.2, 
             metric = "euclidean"
             #metric = "cosine"
             #metric = "correlation"
) 


plotDR(sce, "UMAP", color_by = "time", facet_by = "condition") +
  theme_classic() + 
  geom_point(size=1.0) +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) 



sce = sce[, which(sce$condition == 'TetOn_TetON_RA'| sce$condition == 'TetOn_TetON_dox')]
sce$condition = droplevels(sce$condition)

plotDR(sce, "UMAP", dims = c(1, 2), color_by = "time", facet_by = "condition", ncol = 2
       #k_pal = c("lightgrey", "cornflowerblue", "navy")
) + coord_fixed(1) + 
  theme_classic() + 
  geom_point(size=1.0) +
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14), 
        legend.text = element_text(size=16), 
        legend.title = element_text(size=20)) 

ggsave(paste0(figureDir, 'image_WT_vsConditions_RA_SHout_UMAPs_time_v3.pdf'), width=10, height = 4) 


plotDR(sce, "UMAP", dims = c(1, 2), color_by = "r2_log", facet_by = "condition", ncol = 2
       #k_pal = c("lightgrey", "cornflowerblue", "navy")
) + coord_fixed(1) + 
  theme_classic() + 
  geom_point(size=1.0) +
  theme(axis.text.x = element_text(angle = 0, size = 14, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 14), 
        legend.text = element_text(size=14), 
        legend.title = element_text(size=20)) 

ggsave(paste0(figureDir, 'image_WT_vsConditions_RA_SHout_UMAP_surfaceSize_v3.pdf'), width=10, height = 4) 


plotDR(sce, "UMAP", dims = c(1, 2), color_by = "genotype_foxa2", facet_by = "condition", ncol = 2
       #k_pal = c("lightgrey", "cornflowerblue", "navy")
) + coord_fixed(1) + 
  geom_point(size=0.8) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12), 
        legend.text = element_text(size=12), 
        legend.title = element_text(size=12)) +
  ggplot2::scale_colour_gradient2(limits = c(0, 1), breaks = c(0, 0.1, 0.2,  0.3, 0.5, 0.7, 0.9),
                                  low = "#313695", mid = "#FFFFBF", high = "#C61010", midpoint = 0.3)

ggsave(paste0(figureDir, 'image_WT_vsConditions_RA_SHout_UMAP_genotype_v2.pdf'), width=10, height = 6) 



#ggsave(paste0(figureDir, 'image_WT_SHout_dlogl_umap_cystSurface_v4.pdf'), width=16, height = 8) 

ggsave(paste0(figureDir, 'image_WT_vsConditions_RA_SHout_umap_cystSurface_v5.pdf'), width=30, height = 10) 

saveRDS(sce, file = paste0("../data/image_SHout_dlogl_sce_wt_conditions_metadata_umap_v5.rds"))

rds = reducedDim(sce, 'UMAP')
head(rds)

which(rds[, 1]<0 & sce$condition == "WT")

sce = sce[, which(rds[, 1] < 5)]
set.seed(1234)
sce <- runDR(sce, "UMAP", #pca = 5,
             cells = NULL, 
             features = "type",
             n_neighbors = 30, scale = FALSE,
             min_dist = 0.1, 
             metric = "euclidean"
             #metric = "cosine"
             #metric = "correlation"
) 


p2 = plotDR(sce, "UMAP", color_by = "r2_log") +
  theme_classic() + 
  geom_point(size=1.5) +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) 

p1 / p2

ggsave(paste0(figureDir, 'image_WT_vsConditions_RA_SHout_umap_cystSurface_v2.pdf'), width=10, height = 8) 

sce <- runDR(sce, "DiffusionMap", ncomponents = 5, cells = NULL, features = "type")

plotDR(sce, dims = c(1, 3),  "DiffusionMap", color_by = "time", facet_by = 'condition')

p1 + 
  theme_classic() + 
  geom_point(size=0.7) +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.4),
        axis.text.y = element_text(angle = 0, size = 12)) 

ggsave(paste0(figureDir, 'image_WT_RA_SHout_test_diffusionMap.pdf'), width=8, height = 6) 

