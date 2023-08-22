run_LIANA_defined_celltype = function(subref,
                                      celltypes,
                                      receiver_cells = NULL,
                                      ntop = 100, 
                                      additionalLabel = '_fixedCelltypes')
{
  require(liana)
  require(Seurat)
  require(scater)
  require(scran)
  # subref = aa; additionalLabel = '_fixedCelltypes'; celltypes = unique(aa$celltypes)
  
  system(paste0('mkdir -p ', paste0(outDir, '/LR_analysis_LIANA_v2')))
  # source('functions_scRNAseq.R') 
 
  sce <- as.SingleCellExperiment(subref)
  colLabels(sce) = as.factor(sce$celltypes)
  rownames(sce) = toupper(rownames(sce))
  
  ave.counts <- calculateAverage(sce, assay.type = "counts")
  
  #hist(log10(ave.counts), breaks=100, main="", col="grey80",
  #     xlab=expression(Log[10]~"average count"))
  
  num.cells <- nexprs(sce, byrow=TRUE)
  #smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
  #              xlab=expression(Log[10]~"average count"))
  
  # detected in >= 5 cells, ave.counts >=5 but not too high
  genes.to.keep <- num.cells > 5 & ave.counts >= 10^-4  & ave.counts <10^2  
  summary(genes.to.keep)
  
  sce <- sce[genes.to.keep, ]
  
  ## run the liana wrap function by specifying resource and methods
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_resources()
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_methods()
  
  liana_test <- liana_wrap(sce,  
                           method = c("natmi", "connectome", "logfc", "sca", "cytotalk"  
                                      #'cellphonedb'
                           ),
                           resource = c("Consensus", 'CellPhoneDB', "OmniPath", "LRdb",
                                        "CellChatDB",  "CellTalkDB"), 
                           assay.type = "logcounts", 
                           idents_col = 'celltypes')
  
  # Liana returns a list of results, each element of which corresponds to a method
  liana_test %>% glimpse
  
  # We can aggregate these results into a tibble with consensus ranks
  saveRDS(liana_test, file = paste0(outDir, '/LR_analysis_LIANA_v2/res_lianaTest_Consensus_v2', 
                                    additionalLabel, '.rds'))
  
  liana_test = readRDS(file = paste0(outDir, '/LR_analysis_LIANA_v2/res_lianaTest_Consensus_v2', 
                                     additionalLabel, '.rds'))
  
  liana_test <- liana_test %>%
    liana_aggregate(resource = 'Consensus')
  
  if(is.na(receiver_cells)){ # loop over all cell type candidates
    celltypes = as.character(celltypes)
    receiver_cells = as.character(celltypes)
  }
  
  ntop = 200
  manually_specifying_sender_receiver = FALSE
  if(manually_specifying_sender_receiver){
    ntop = 100
    sender_cells = c('3', '4')
    receiver_cells = sender_cells
    
    liana_test %>%
      liana_dotplot(source_groups = sender_cells,
                    target_groups = receiver_cells,
                    ntop = ntop)
    ggsave(filename = paste0(outDir, '/LR_analysis_LIANA_v2/liana_LR_prediction_cluster3.vs.cluster4', 
                             additionalLabel, 
                             #'_receiverCells.', receiver_cells[m], 
                             '_ntop.', ntop, '.pdf'), 
           width = 15, height = 0.25*ntop, limitsize = FALSE)
    
    ntop = 500
    sender_cells = c('0', '6', '2', '8')
    receiver_cells = sender_cells
    
    liana_test %>%
      liana_dotplot(source_groups = sender_cells,
                    target_groups = receiver_cells,
                    ntop = ntop)
    ggsave(filename = paste0(outDir, '/LR_analysis_LIANA_v2/liana_LR_prediction_cluster_0_6_2_8', 
                             additionalLabel, 
                             #'_receiverCells.', receiver_cells[m], 
                             '_ntop.', ntop, '.pdf'), 
           width = 25, height = 0.25*ntop, limitsize = FALSE)
    
    
  }
  
  for(m in 1:length(receiver_cells))
  {
    # m = 1
    cat(m, '-- receiver cells : ', receiver_cells[m], '\n')
    #liana_test %>%
    #  liana_dotplot(source_groups = celltypes[n],
    #                target_groups = celltypes,
    #                ntop = ntop)
    liana_test %>%
      liana_dotplot(source_groups = celltypes,
                    target_groups = receiver_cells[m],
                    ntop = ntop)
    #liana_test_save =  liana_test %>% filter()
    #  liana_dotplot(source_groups = celltypes,
    #                target_groups = receiver_cells[m],
    #                ntop = ntop)
    
    ggsave(filename = paste0(outDir, '/LR_analysis_LIANA/liana_LR_prediction_recieveCell', 
                             additionalLabel, 
                             '_receiverCells.', receiver_cells[m], 
                             '_ntop.', ntop, '.pdf'), 
           width = 30, height = 0.25*ntop, limitsize = FALSE)
  }
  
  
  pdfname = paste0(outDir, '/LR_analysis_LIANA/liana_celltype_communication_freqHeatmap_', ntop, 
                   additionalLabel, '.pdf')
  pdf(pdfname, width=20, height = 8)
  
  liana_trunc <- liana_test %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01) # this can be FDR-corr if n is too high
  
  heat_freq(liana_trunc)
  
  dev.off()
  
  
}


##########################################
# NichetNet analysis
##########################################
Run_NicheNet_Seurat = function()
{
  library(pryr) # monitor the memory usage
  require(ggplot2)
  library(nichenetr)
  library(Seurat) # please update to Seurat V4
  library(tidyverse)
  library(circlize)
  library(RColorBrewer)
  require(scran)
  require(scater)
  library(nichenetr)
  library(Seurat)
  library(tidyverse)
  library(knitr)
  library(tidyr)
  library(reshape2)
  library(gridExtra)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(plotly)
  library(ggalluvial)
  library(stringr)
  
  options(stringsAsFactors = F)
  
  dataPath_nichenet = '/users/jingkui.wang/projects/heart_regeneration/data/NicheNet/'
  out_Nichenet = paste0(outDir, '/LR_analysis_NicheNet')
  system(paste0('mkdir -p ', out_Nichenet))
  
  
  mem_used()
  
  ##########################################
  # load processed scRNA-seq and visium data
  ##########################################
  refs = aa
  refs$celltypes = refs$clusters
  
  table(refs$celltypes)
  Idents(refs) = as.factor(refs$celltypes)
  
  ########################################################
  ########################################################
  # Section : test default NicheNet with Seurat object 
  # 
  ########################################################
  ########################################################
  
  ##########################################
  # step 0):  load the Nichenet data
  ##########################################
  # NicheNet’s ligand-target prior model
  ligand_target_matrix = readRDS(paste0(dataPath_nichenet,  "ligand_target_matrix.rds"))
  ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
  
  # ligand-receptor network, Putative ligand-receptor links from NicheNet
  lr_network = readRDS(paste0(dataPath_nichenet, "lr_network.rds"))
  lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
    distinct(ligand, receptor, bonafide)
  
  # If wanted, users can remove ligand-receptor interactions that were predicted based on 
  # protein-protein interactions and 
  # only keep ligand-receptor interactions that are described in curated databases.
  # lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  head(lr_network)
  
  ligands = lr_network %>% pull(ligand) %>% unique()
  receptors = lr_network %>% pull(receptor) %>% unique()
  
  # ## weighted integrated networks 
  weighted_networks = readRDS(paste0(dataPath_nichenet,  "weighted_networks.rds"))
  head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
  head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
  
  weighted_networks_lr = weighted_networks$lr_sig %>% 
    dplyr::rename(ligand = from, receptor = to) %>%
    inner_join(lr_network %>% 
                 distinct(ligand,receptor), by = c("ligand","receptor"))
  
  ##########################################
  # we directly consider time-specific 
  # loop over the cell types for each time points
  # here we are using the basic function predict_ligand_activities in NicheNet
  # Shoval is also using the same one, more flexible to adapt for data 
  # original code from https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
  ##########################################
  celltypes = levels(refs$celltypes)
  cat(celltypes, '\n')
  
  celltypes = c('0', "3", "6", "2", "4", "8")
  control_cells = c('1')
  receiver_cells = celltypes
  
  celltypes_sel = unique(c(celltypes, receiver_cells, control_cells)) 
  
  subref = subset(refs, cells = colnames(refs)[!is.na(match(refs$celltypes, celltypes_sel))])
  subref$celltypes = droplevels(as.factor(subref$celltypes))
  table(subref$celltypes)
  
  #subref = subset(x = subref, downsample = 1000)
  cat('celltype to consider -- ', names(table(subref$celltypes)), '\n')
  Idents(subref) = subref$celltypes
  
  # process and clean the gene names
  Gene.filtering.preprocessing = TRUE
  if(Gene.filtering.preprocessing){
    sce <- as.SingleCellExperiment(subref)
    ggs = rownames(sce)
    
    colLabels(sce) = as.factor(sce$celltypes)
    rownames(sce) = toupper(rownames(sce))
    
    ave.counts <- calculateAverage(sce, assay.type = "counts")
    
    #hist(log10(ave.counts), breaks=100, main="", col="grey80",
    #     xlab=expression(Log[10]~"average count"))
    
    num.cells <- nexprs(sce, byrow=TRUE)
    #smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
    #              xlab=expression(Log[10]~"average count"))
    
    # detected in >= 10 cells, ave.counts cutoff are both very lose
    genes.to.keep <- num.cells > 10 & ave.counts >= 10^-5  & ave.counts <10^3  
    summary(genes.to.keep)
    sce <- sce[genes.to.keep, ]
    ggs = ggs[genes.to.keep]
    
    rownames(sce) = make.names(rownames(sce), unique = TRUE)
    
    #geneDup = data.frame(gene = ggs, gg.uniq = rownames(sce), stringsAsFactors = FALSE)
    #geneDup$geneSymbol = get_geneName(geneDup$gene)
    
    #kk = which(geneDup$geneSymbol != geneDup$gg.uniq)
    
    subref = as.Seurat(sce, counts = "counts", data = "logcounts")
    rm(sce)
    Idents(subref) = as.factor(subref$celltypes)
    
    #saveRDS(geneDup, paste0(outDir, '/geneSymbol_duplication_inLRanalysis.rds'))
    
  }
  
  # step 1:  Define a “sender/niche” cell population and a “receiver/target” cell population 
  # present in your expression data and determine which genes are expressed in both populations
  if(is.na(receiver_cells)){ # loop over all cell type candidates
    receiver_cells = celltypes
  }
  
  # to save the result for each time point
  expressed_genes = c()
  res = tibble(sample=character(), sender=character(), receiver=character(),  
               test_ligand=character(), auroc=double(), aupr=double(), 
               pearson = double(),  rank=integer())
  
  for(receiver in receiver_cells) # loop over receiver cells 
  {
    # specify receiver
    # receiver = receiver_cells[1]
    cat('-- receiver : ', receiver, '-- \n')
    expressed_genes_receiver = get_expressed_genes(receiver, subref, pct = 0.10)
    background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    
    # step 2: Define a gene set of interest: these are the genes in the “receiver/target” cell population 
    # that are potentially affected by ligands expressed by interacting cells 
    # (e.g. genes differentially expressed upon cell-cell interaction)
    DE_table_receiver = FindMarkers(object = subref, 
                                    ident.1 = receiver, 
                                    ident.2 = control_cells, 
                                    min.pct = 0.10) %>% rownames_to_column("gene")
    
    geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
    geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
    
    ## sender
    sender_celltypes = celltypes
    #loop over the same receiver for sender cells, i.e. incl. self-activation (shoval's example)
    for (sender in sender_celltypes){
      # sender = sender_celltypes[1]
      cat('-- sender : ', sender, '-- \n')
      expressed_genes_sender = get_expressed_genes(sender, subref, pct=0.1)
      
      # Step 3: Define a set of potential ligands
      # these are ligands that are expressed by the “sender/niche” cell population and 
      # bind a (putative) receptor expressed by the “receiver/target” population
      # Get potential ligands in our datasets
      expressed_ligands = intersect(ligands, expressed_genes_sender)
      expressed_receptors = intersect(receptors,expressed_genes_receiver)
      
      potential_ligands = lr_network %>% filter(ligand %in% expressed_ligands & receptor %in% expressed_receptors) %>%
        pull(ligand) %>% unique()
      
      # Step 4: Perform NicheNet’s ligand activity analysis on the gene set of interest
      # rank the potential ligands based on the presence of their target genes in the gene set of interest 
      # (compared to the background set of genes)     
      # the main function used : predict_ligand_activities
      ligand_activities = predict_ligand_activities(geneset = geneset_oi,
                                                    background_expressed_genes = background_expressed_genes,
                                                    ligand_target_matrix = ligand_target_matrix, 
                                                    potential_ligands = potential_ligands) 
      ligand_activities = ligand_activities %>% 
        arrange(-pearson) %>% dplyr::mutate(sender=sender, receiver=receiver, 
                                            rank = rank(dplyr::desc(pearson)))
      
      # pearson correlation between the target expresson and prediciton as ligand activity scores 
      res = rbind(res, ligand_activities) 
    }
    
    # list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
    # expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  }
  
  # processing after the ligand activity prediction
  res = res %>% dplyr::mutate(interaction=paste0(sender,"-",receiver))
  #res = res %>% filter(test_ligand %in% growth.factors) %>% dplyr::mutate(interaction=paste0(sender,"-",receiver))
  
  res %>% DT::datatable()
  expressed_genes = unique(expressed_genes)
  
  ##########################################
  # ligand visualization
  ##########################################
  res = res %>% dplyr::mutate(interaction=paste0(sender,"-",receiver))
  
  # Plot ligand activities scores in general
  p1 = res %>% ggplot(aes(x=pearson)) + 
    geom_histogram(color="black", fill="darkorange") + 
    geom_vline(aes(xintercept=0.10), color="red", linetype="dashed", size=1) +  
    labs(x="ligand activity (Pearson Correlation)", y = "# ligands") +  
    theme_classic() + ggtitle("Ligands activity scores")
  plot(p1)
  
  ggsave(filename = paste0(out_Nichenet, '/ligand_pearson_correlation.pdf'), 
         width = 10, height = 8)
  
  # Plot ligand activities scores per interaction
  p2 = res %>% ggplot(aes(x=pearson)) + 
    geom_histogram(color="black", fill="darkorange") + 
    geom_vline(aes(xintercept=0.1), color="red", linetype="dashed", size=1) +  
    labs(x="ligand activity (Pearson Correlation)", y = "# ligands") +  
    theme_classic() + 
    ggtitle("Ligands activity scores") + 
    facet_wrap(.~interaction, )
  plot(p2)
  ggsave(filename = paste0(out_Nichenet, '/ligand_pearson_correlation_individualInteraction.pdf'), 
         width = 12, height = 10)
  
  # Visualize ligand expression along different days
  ligand.expression.df = AverageExpression(subref, features =unique(res$test_ligand))[['RNA']]
  pseudo_signal = ceiling(log2(range(ligand.expression.df[ligand.expression.df>0]))[1])-1
  ligand.expression.df = log2(2^pseudo_signal +ligand.expression.df)[sort(rownames(ligand.expression.df)), 
                                                                     sort(colnames(ligand.expression.df))]
  write.csv2(inner_join(res, ligand.expression.df %>% as.data.frame %>% rownames_to_column('test_ligand')), 
             file = paste0(out_Nichenet, '/predicted_ligand_activity_avExpr_for_receivers_senderss.csv'),
             row.names = FALSE)
  
  #p2 = Heatmap(ligand.expression.df, cluster_columns = F, column_title = "Ligand Expression", name="Expression (log2)", row_order = row_order(p1))
  
  # Visualization top ligand per interaction
  p1 = res %>% top_n(1000, pearson) %>% transmute(test_ligand=test_ligand, x=interaction, pearson=pearson) %>% 
    spread(x, pearson)  %>% 
    replace(is.na(.), 0) %>% 
    ungroup() %>% 
    column_to_rownames('test_ligand')
  p2 = ligand.expression.df[match(rownames(p1), rownames(ligand.expression.df)), ]
  
  #ha.col = list(Interaction=setNames(c(brewer.pal(9, 'YlOrRd')), colnames(p1)))
  ha.col = c(brewer.pal(9, 'YlOrRd'))[1:length(colnames(p1))]
  names(ha.col) = colnames(p1)
  ha.col = as.list(ha.col)
  ha = HeatmapAnnotation(Interaction=colnames(p1), annotation_name_side = 'left', col=ha.col)
  
  p3 = p1 %>% Heatmap(
    column_title = "Average Ligand pearson correlation\n(NAs replaced with zeros)", 
    cluster_columns = TRUE, show_row_names = F, show_column_names = , 
    top_annotation = ha, col=brewer.pal(4, 'Oranges'), 
    name="Ligand Score(pearson)")
  
  ha = colnames(p2)
  ha = HeatmapAnnotation(Cluster=ha, col=ha.col)
  p4 = Heatmap(p2, 
               cluster_columns = F, show_column_names = F,
               column_title = "Ligand Expression", name="log10(Expression)", 
               row_order = row_order(p3), 
               top_annotation = ha)
  
  pdfname = paste0(out_Nichenet, '/heatmap_top_ligand_pearson_correlation_individualInteraction_v2.pdf')
  pdf(pdfname, width=16, height = 25)
  
  pp = p3 + p4
  plot(pp)
  dev.off()
  
  
}




