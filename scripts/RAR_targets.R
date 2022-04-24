##########################################################################
##########################################################################
# Project: RA competence with Hannah
# Script purpose: analyze the ChIP-seq, ChIA-Pet and RNA-seq, microarray data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Aug 24 18:47:03 2021
##########################################################################
##########################################################################
rm(list = ls())
resDir = '../results/RAR_targets'
RdataDir = '../results/RAR_targets/Rdata'
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)


########################################################
########################################################
# Section I : Run CID for ChIA-Pet
# 
########################################################
########################################################
Run.CID.ChIAPet = FALSE
if(Run.CID.ChIAPet){
  outDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/publishedDataset/ChIPSeq/'
  
  xx = read.delim(paste0(outDir, 'filereport_read_run_PRJNA123691_tsv.txt'), header = TRUE)
  
  xx = data.frame(xx$sample_title, xx$run_accession, stringsAsFactors = FALSE)
  xx$control = ''
  xx$SRC = paste0('/scratch/jiwang/RA_chipseq/ngs_raw/', xx[,2], '.fastq')
  colnames(xx) = c('CONDITION','SAMPLEID','SAMPLEID_CONTROL','SRC')
  xx$CONDITION = gsub('-','.', xx$CONDITION)
  xx$CONDITION = gsub('[+]', '_', xx$CONDITION)
  
  xx = xx[c(13,14, 1:12, 15:19), ]
  xx$SAMPLEID = as.character(xx$SAMPLEID)
  
  xx$SAMPLEID_CONTROL[-1] = xx$SAMPLEID[1]
  
  write.csv(xx, file = paste0(outDir, 'paramfile.csv'), quote = FALSE, row.names = FALSE)
  xx = xx[, c(2, 1)]
  colnames(xx) = c('sampleID', 'fileName')
  write.csv(xx, file = paste0(outDir, 'sampleInfos.csv'), quote = FALSE, row.names = FALSE)
  
  
  library(MICC)
  bedpe = '/Volumes/clustertmp/jiwang/RA_chiapet/test_cid_v2.bedpe'
  cid.bedpe <- read.table(bedpe, sep = "\t", header = FALSE)
  
  library(tictoc)
  tic()
  MICCoutput(cid.bedpe, "../results/micc_out.txt")
  toc()
  
}

########################################################
########################################################
# Section : Elena's time sereis RNA-seq data first two time points
# before RA and -10h with and without RA
########################################################
########################################################
Reanalyze.RNAseq = FALSE
if(Reanalyze.RNAseq){
  
  
  load(file = 
  paste0('/Users/jiwang/workspace/imp/organoid_patterning/results/Rdata/RNAseq_timeSeries_sortedDay5_count_design_geneSymbol.Rdata'))
  
  # compare the first time points 8h after RA treatment
  #jj = which(design$condition == 'before_RA'| design$condition == 'min10h' | design$condition == 'min10h_RA')
  
  # 8h and 18h after RA treatment
  jj = which(design$condition == 'before_RA'| design$condition == 'min10h' | design$condition == 'min10h_RA'|
               design$condition == '0h'| design$condition == '0h_RA')
  
  design = design[jj, ]
  counts = counts[, jj]
  
  # DESeq test
  require(DESeq2)
  require(ggplot2)
  require(dplyr)
  require(pheatmap)
  library(ggrepel)
  library(tibble)
  
  
  dds <- DESeqDataSetFromMatrix(counts, DataFrame(design), design = ~ condition)
  dds$condition = relevel(dds$condition, ref = 'min10h')
  ss = rowSums(counts(dds))
  
  hist(log10(ss), breaks = 200, main = 'log2(sum of reads for each gene)')
  
  cutoff.peak = 100
  cat(length(which(ss > cutoff.peak)), 'peaks selected \n')
  gg.tokeep = 'Nog|Acvrl1|Acvr1|Notch|Dll|Jagn|Gdf|Bmp|Fgf|Wnt'
  
  sels = unique(c(which(ss > cutoff.peak), grep(gg.tokeep, rownames(dds))))
  
  dds <- dds[sels, ]
  
  
  # normalization and dimensionality reduction
  dds <- estimateSizeFactors(dds)
  #fpm = fpm(dds, robust = TRUE)
  ss = colSums(counts(dds))
  plot(sizeFactors(dds), ss/10^6)
  
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  
  ## there is a batch different for sample S10131
  kk = c(1:nrow(design))
  pca=plotPCA(vsd[,kk], intgroup = c('condition'), returnData = TRUE, ntop = 500)
  pca2save = as.data.frame(pca)
  ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition))  + 
    geom_point(size=3) + 
    geom_text(hjust = 0.3, nudge_y = 0.4, size=4)
  
  plot(ggp)
  
  
  dds <- estimateDispersions(dds, fitType = 'parametric')
  plotDispEsts(dds, ymin = 10^-3, main = 'RA')
  
  
  dds <- nbinomWaldTest(dds)
  resultsNames(dds)  
  
  res = results(dds, contrast=c("condition", 'min10h_RA', 'min10h'), alpha = 0.05)
  res <- lfcShrink(dds, coef="condition_min10h_RA_vs_min10h", type="normal")
  plotMA(res)
  
  res2 = results(dds, contrast=c("condition", '0h_RA', '0h'), alpha = 0.05)
  res2 <- lfcShrink(dds, contrast=c("condition", '0h_RA', '0h'))
  plotMA(res2)
  
  colnames(res2) = paste0(colnames(res), '.18h')
  
  xx = data.frame(as.data.frame(res), as.data.frame(res2))
  
  
  res = xx
  #gg.signif = rownames(res1)[which(res1$padj < 0.05)]
  # save the average expression of Foxa2 positive and negative cells
  cpm = fpm(dds)
  #cpm = cpm[, c(grep('pos', colnames(cpm)), grep('neg', colnames(cpm)))]
  
  res = data.frame(res, cpm, stringsAsFactors = FALSE)
  res = res[order(res$padj), ]
  
  #res = res[which(res$padj < 0.05), ]
  #res1 = data.frame(res1, cpm.foxa2.pos = apply(as.matrix(cpm[, grep('pos', colnames(cpm))]), 1, mean), 
  #                  cpm.foxa2.neg = apply(as.matrix(cpm[, grep('neg', colnames(cpm))]), 1, mean))
  
  res$gene = rownames(res)
  res$padj = -log10(res$padj)
  res$padj.18h = -log10(res$padj.18h)
  
  examples.sel = unique(rownames(res)[which(res$padj > 50 | res$log2FoldChange > 2.)])
  
  ggplot(data=res, aes(x=log2FoldChange, y=padj, label = gene)) +
    geom_point(size = 1.5) + 
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) +
    geom_text_repel(data=subset(res, abs(log2FoldChange) > 2 | padj > 50), size = 3.5) +
    #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel),
    #                 size = 4) + 
    #scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-2.5, 2.5), col="red") +
    geom_hline(yintercept=50, col="red") 
  
  ggsave(paste0(resDir, "/DESeq2_DEgenes_RA.vs.noRA.day2.8h_v2.pdf"), width=16, height = 10)
  
  
  ggplot(data=res, aes(x=log2FoldChange.18h, y=padj.18h, label = gene)) +
    geom_point(size = 1.5) + 
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12)) +
    geom_text_repel(data=subset(res, abs(log2FoldChange.18h) > 2.5 | padj.18h > 50), size = 3.5) +
    #geom_label_repel(data=  as.tibble(res) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel),
    #                 size = 4) + 
    #scale_color_manual(values=c("blue", "black", "red")) +
    geom_vline(xintercept=c(-2.5, 2.5), col="red") +
    geom_hline(yintercept=50, col="red") 
  
  ggsave(paste0(resDir, "/DESeq2_DEgenes_RA.vs.noRA.day2.18h_v2.pdf"), width=16, height = 10)
  
  
  # compare the logFC
  ggplot(res, aes(x=log2FoldChange, y=log2FoldChange.18h, label = gene)) +
    geom_point(size = 0.7) +
    geom_text_repel(data=subset(res, abs(log2FoldChange) > 2.5 | padj > 50), size = 3.5) +
    geom_abline(slope = 1, intercept = 0, colour = "red") + 
    geom_vline(xintercept = c(2.5, -2.5), colour = "darkgray") + 
    geom_hline(yintercept = c(2.5, -2.5), colour = "darkgray") + 
    ggtitle('log2FC 8h vs 18h') 
    
  ggsave(paste0(resDir, "/DESeq2_DEgenes_log2FC_8h.vs.18h.pdf"), width=16, height = 10)
  
  
  ggplot(res, aes(x=padj, y=padj.18h, label = gene)) +
    geom_point(size = 0.7) +
    geom_text_repel(data=subset(res, abs(log2FoldChange) > 2.5 | padj > 30), size = 3.5) +
    geom_abline(slope = 1, intercept = 0, colour = "red") + 
    geom_vline(xintercept = c(30), colour = "darkgray") + 
    geom_hline(yintercept = c(30), colour = "darkgray") + 
    ggtitle('FDR 8h vs 18h') 
  
  ggsave(paste0(resDir, "/DESeq2_DEgenes_FDR_8h.vs.18h.pdf"), width=16, height = 10)
  
  # ggplot(data = yy,  aes(y = pvalue_neg, x = lfc_neg,  )) + 
  #   geom_point(size = 1) + 
  #   labs(title = paste0(cc[n], " - neg  "), x = '', y = '-log10(pval)') + 
  #   theme(axis.text.x = element_text(size = 12), 
  #         axis.text.y = element_text(size = 12)) + 
  #   geom_hline(yintercept = 2, colour = "red") +
  #   #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
  #   geom_text_repel(data=subset(yy, pvalue_neg > 2), size = 4)
  
  write.csv(res, file = paste0(resDir, '/DESeq2_DEgenes_fdr.0.05_normalizedData_8h.18h.csv'), row.names = TRUE, quote = FALSE)
  
}

########################################################
########################################################
# Section : RAR chipseq and ChIApet and main assign peaks to genes
# 
########################################################
########################################################
peak.to.gene.assignment = FALSE
if(peak.to.gene.assignment){
  
  library("ChIPseeker");
  library("rtracklayer")
  require(ChIPpeakAnno)
  
  # read peak file
  peakDir = '/Volumes//groups/tanaka/People/current/jiwang/projects/RA_competence/publishedDataset/chipseq/calledPeaks/macs2/'
  peakfile = paste0(peakDir, 'RAR_Day2_8hrsRA_ChIP.seq_SRR039161.SRR039162.SRR039163.SRR039164.SRR039165.SRR039166._merged_macs2_peaks.xls')
  
  p = readPeakFile(peakfile, as = "GRanges")
  
  p = p[!is.na(match(seqnames(p), c(paste0('chr', c(1:19)), 'chrX', 'chrY')))]
  
  # define the data frame to save
  keep = data.frame(p ,stringsAsFactors = FALSE)
  keep$name = gsub('RAR_Day2_8hrsRA_ChIP.seq_SRR039161.SRR039162.SRR039163.SRR039164.SRR039165.SRR039166._merged_', '',
                   keep$name)
  colnames(keep)[ncol(keep)] = 'peak.name'
  rownames(keep) = paste0(keep$seqnames, ':', keep$start, '-', keep$end)
  names(p) = rownames(keep)
  
  keep$promoter = NA
  keep$geneBody = NA
  keep$downstream = NA
  keep$chiapet = NA
  keep$distalIntergenic_closestTSS = NA
  
  ##########################################
  #  gene annotation to use
  ##########################################
  #library("TxDb.Mmusculus.UCSC.mm10.ensGene")
  library(org.Mm.eg.db)
  data("TSS.mouse.GRCm38") # this data is for the gene coordinates
  
  seqlevelsStyle(TSS.mouse.GRCm38) <- "UCSC"
  
  Select.protein.coding.genes = TRUE
  if(Select.protein.coding.genes){
    annot = read.delim(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/', 
                              'ens_BioMart_GRCm38.p6.txt'), sep = '\t', header = TRUE)
    annot = annot[which(annot$Gene.type == 'protein_coding'), ]
    annot = annot[!is.na(match(annot$Chromosome.scaffold.name, as.character(c(1:19, 'X', 'Y')))), ]
    
    mm = match(names(TSS.mouse.GRCm38), annot$Gene.stable.ID)
    cat(length(which(!is.na(mm))), ' protein coding genes in the TSS.mouse.GRCm38 \n')
    
    TSS.mouse.GRCm38 = TSS.mouse.GRCm38[!is.na(mm)]
  }
  
  # extend gene annotation to promoter 1000bp and 300 downstream for chiapet pair assignment
  genes = TSS.mouse.GRCm38
  genes$feature = names(genes)
  genes = addGeneIDs(genes, orgAnn="org.Mm.eg.db", 
                     IDs2Add="symbol", feature_id_type = 'ensembl_gene_id')
  
  #xx = genes
  start(genes) <- start(genes) - 1000
  end(genes) <- end(genes) + 300
  strand(genes) = '*'
  
  # quick check the feature distribution
  #annoData <- toGRanges(TxDb.Mmusculus.UCSC.mm10.ensGene, feature="gene")
  #annoData[1:2]
  
  #peakAnnots= annotatePeak(p, TxDb=TxDb.Mmusculus.UCSC.mm10.ensGene, tssRegion = c(-2000, 2000)) 
  #print(plotAnnoBar(peakAnnots))
  #print(plotDistToTSS(peakAnnots))
  
  ##########################################
  # first round: around promoters
  ##########################################
  # define promoter regions -5000 and 500bp around TSS
  annotPromoters <- promoters(TSS.mouse.GRCm38, upstream=2000, downstream=500)
  
  # save all promoters overlapping the peaks
  peakAnnots <- annotatePeakInBatch(p, AnnotationData=annotPromoters, output = 'overlapping', multiple = TRUE, 
                                    select = 'all')
  peakAnnots = addGeneIDs(annotatedPeak=peakAnnots, 
                          orgAnn="org.Mm.eg.db", 
                          IDs2Add="symbol")
  
  pa = data.frame(peakAnnots, stringsAsFactors = FALSE)
  pa = pa[which(!is.na(pa$feature)), ]
  
  assignedPeaks = unique(pa$peak)
  for(n in 1:length(assignedPeaks))
  {
    # n = 4
    ii = which(rownames(keep) == assignedPeaks[n]) 
    kk = which(pa$peak == assignedPeaks[n])
    
    if(length(kk) >1 ) cat(ii, '\n')
    gene.symbols = pa$symbol[kk]
    gene.ens = pa$feature[kk]
    
    gene.symbols[which(is.na(gene.symbols))] = gene.ens[which(is.na(gene.symbols))]
    
    keep$promoter[ii] = paste0(gene.symbols, collapse = ';')
  }
  
  cat(length(which(!is.na(keep$promoter))), 'peaks assigned to promoters \n')
  
  
  ##########################################
  # second round of peak-to-gene assignment: gene body 
  ##########################################
  peakAnnots <- annotatePeakInBatch(p, AnnotationData=TSS.mouse.GRCm38, output = 'overlapping', multiple = TRUE, 
                                    select = 'all')
  peakAnnots = addGeneIDs(annotatedPeak=peakAnnots, 
                          orgAnn="org.Mm.eg.db", 
                          IDs2Add="symbol")
  
  pa = data.frame(peakAnnots, stringsAsFactors = FALSE)
  pa = pa[which(!is.na(pa$feature)), ]
  
  assignedPeaks = unique(pa$peak)
  for(n in 1:length(assignedPeaks))
  {
    # n = 4
    ii = which(rownames(keep) == assignedPeaks[n]) 
    kk = which(pa$peak == assignedPeaks[n])
    
    if(length(kk) >1 ) cat(ii, '\n')
    gene.symbols = pa$symbol[kk]
    gene.ens = pa$feature[kk]
    
    gene.symbols[which(is.na(gene.symbols))] = gene.ens[which(is.na(gene.symbols))]
    
    keep$geneBody[ii] = paste0(gene.symbols, collapse = ';')
  }
  
  cat(length(which(!is.na(keep$geneBody))), 'peaks assigned to gene bodies \n')
  cat(length((which(is.na(keep$promoter) & is.na(keep$geneBody) & is.na(keep$chiapet) & is.na(keep$downstream)))), 
      ' peaks without assignments \n')
  
  ##########################################
  # third round of peak-to-gene assignment: downstream
  ##########################################
  peakAnnots <- annotatePeakInBatch(p, AnnotationData=TSS.mouse.GRCm38, output = 'downstream', multiple = TRUE, maxgap = 300,
                                    select = 'all')
  peakAnnots = addGeneIDs(annotatedPeak=peakAnnots, 
                          orgAnn="org.Mm.eg.db", 
                          IDs2Add="symbol")
  
  pa = data.frame(peakAnnots, stringsAsFactors = FALSE)
  pa = pa[which(!is.na(pa$feature)), ]
  
  assignedPeaks = unique(pa$peak)
  
  for(n in 1:length(assignedPeaks))
  {
    # n = 4
    ii = which(rownames(keep) == assignedPeaks[n]) 
    kk = which(pa$peak == assignedPeaks[n])
    
    if(length(kk) >1 ) cat(ii, '\n')
    gene.symbols = pa$symbol[kk]
    gene.ens = pa$feature[kk]
    
    gene.symbols[which(is.na(gene.symbols))] = gene.ens[which(is.na(gene.symbols))]
    
    keep$downstream[ii] = paste0(gene.symbols, collapse = ';')
  }
  
  cat(length(which(!is.na(keep$downstream))), 'peaks assigned to downstream \n')
  
  cat(length((which(is.na(keep$promoter) & is.na(keep$geneBody) & is.na(keep$chiapet) & is.na(keep$downstream)))), 
      ' peaks without assignments \n')
  
  ##########################################
  # distal intergenic peaks assignment with ChIA-Pet 
  ##########################################
  library(MICC)
  chiapetDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/publishedDataset/chiapet/'
  Save.significant.chiapet.interactions = TRUE
  
  fdr.cutoff = 0.2
  
  Rerun.MICC.thresholding = FALSE
  if(Rerun.MICC.thresholding){
    # cid output wiht --micc=2, 3, 4 are the same (https://groups.csail.mit.edu/cgs/gem/cid/)
    cid.output = paste0(chiapetDir, 'out_cid_micc_mim2PET.bedpe')
    cid.bedpe <- read.table(cid.output, sep = "\t", header = FALSE)
    
    
    cid.micc.out.file = paste0(chiapetDir, 'cid_out_min2PET_micc_out_min5PET.txt')
    MICCoutput(cid.bedpe, cid.micc.out.file, MinConfident = 5)
    
    mc = read.table(file = cid.micc.out.file, header = TRUE)
    length(which(mc$fdr<0.05))
    length(which(mc$fdr<0.1))
    length(which(mc$fdr<0.15))
    length(which(mc$fdr<0.2))
    
    
    mc = mc[which(mc$fdr < fdr.cutoff), ]
    
    saveRDS(mc, file = paste0(RdataDir, '/chiapet_interactions_cid_min2PET_micc_min5PET_fdr_', fdr.cutoff, '.rds'))
    
    if(Save.significant.chiapet.interactions){
      write.table(mc, file = paste0(chiapetDir, 'chiapet_interactions_cid_min2PET_micc_min5PET_fdr_', fdr.cutoff, '.bedpe'), 
                  sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
    }
    
  }else{
    mc = readRDS(file = paste0(RdataDir, '/chiapet_interactions_cid_min2PET_micc_min5PET_fdr_', fdr.cutoff, '.rds'))
  }
  
  
  library("GenomicFeatures")
  library("GenomicRanges")
  
  mc_p1 = data.frame(mc[, c(1:3)], stringsAsFactors = FALSE)
  mc_p2 = data.frame(mc[, c(4:6)], stringsAsFactors = FALSE)
  mc_p1$strand = '*'
  mc_p2$strand = '*'
  colnames(mc_p1) = c('seqnames', 'start', 'end', 'strand')
  colnames(mc_p2) = colnames(mc_p1)
  
  mc_p1 = makeGRangesFromDataFrame(mc_p1, seqnames.field=c("seqnames"),
                                 start.field="start", end.field=c("end"), strand.field="strand")  
  mc_p2 = makeGRangesFromDataFrame(mc_p2, seqnames.field=c("seqnames"),
                                   start.field="start", end.field=c("end"), strand.field="strand")  
  
  #gap.size = -1
  hits1 = data.frame(findOverlaps(p, mc_p1, type = 'any', select = 'all'), maxgap=-1L)
  hits2 = data.frame(findOverlaps(p, mc_p2, type = 'any', select = 'all'), maxgap=-1L)
  
  for(ii1 in unique(hits1$queryHits))
  {
    # ii1 = 4430
    index_pet = unique(hits1$subjectHits[which(hits1$queryHits == ii1)])
    distal.genes = findOverlaps(mc_p2[index_pet], genes, type = 'any', select = 'all')
    
    index_found = unique(distal.genes@to)
    
    if(length(index_found) > 0){
      gene.symbols = genes$symbol[index_found]
      gene.ens = genes$feature[index_found]
      gene.symbols[which(is.na(gene.symbols))] = gene.ens[which(is.na(gene.symbols))]
      gene.symbols = unique(gene.symbols)
      cat(ii1, '\n')
      if(!is.na(keep$chiapet[ii1])){
        keep$chiapet[ii1] = paste0(c(keep$chiapet[ii1], gene.symbols), collapse = ';')
      }else{
        keep$chiapet[ii1] = paste0(gene.symbols, collapse = ';')
      }
    }
  }
  
  for(ii2 in unique(hits2$queryHits))
  {
    # ii1 = 4430
    index_pet = unique(hits2$subjectHits[which(hits2$queryHits == ii2)])
    distal.genes = findOverlaps(mc_p1[index_pet], genes, type = 'any', select = 'all')
    index_found = unique(distal.genes@to)
    
    if(length(index_found) > 0){
      gene.symbols = genes$symbol[index_found]
      gene.ens = genes$feature[index_found]
      gene.symbols[which(is.na(gene.symbols))] = gene.ens[which(is.na(gene.symbols))]
      gene.symbols = unique(gene.symbols)
      cat(ii2, '\n')
      if(!is.na(keep$chiapet[ii2])){
        keep$chiapet[ii2] = paste0(c(keep$chiapet[ii2], gene.symbols), collapse = ';')
      }else{
        keep$chiapet[ii2] = paste0(gene.symbols, collapse = ';')
      }
    }
  }
  
  cat(length(which(!is.na(keep$chiapet))), 'peaks assigned with chiapet \n')
  
  cat(length((which(is.na(keep$promoter) & is.na(keep$geneBody) & is.na(keep$chiapet) & is.na(keep$downstream)))), 
      ' peaks without assignments \n')
  
  ##########################################
  # last round: closest TSS if the peaks are still not assigned  
  ##########################################
  peakAnnots <- annotatePeakInBatch(p, AnnotationData=annotPromoters, output = 'nearestLocation', multiple = FALSE)
  peakAnnots = addGeneIDs(annotatedPeak=peakAnnots, 
                          orgAnn="org.Mm.eg.db", 
                          IDs2Add="symbol")
  
  pa = data.frame(peakAnnots, stringsAsFactors = FALSE)
  pa = pa[which(!is.na(pa$feature)), ]
  
  index_leftover = (which(is.na(keep$promoter) & is.na(keep$geneBody) & is.na(keep$chiapet) & is.na(keep$downstream)))
  
  #assignedPeaks = unique(pa$peak)
  for(kk in index_leftover)
  {
    # n = 4
    #ii = which(rownames(keep) == assignedPeaks[n]) 
    #kk = which(pa$peak == assignedPeaks[n])
    gene.symbols = pa$symbol[kk]
    gene.ens = pa$feature[kk]
    gene.symbols[which(is.na(gene.symbols))] = gene.ens[which(is.na(gene.symbols))]
    
    #gene.symbols = pa$symbol[kk]
    #gene.ens = pa$feature[kk]
    keep$distalIntergenic_closestTSS[kk] = paste0(gene.symbols, collapse = ';')
    
  }
  
  cat(length(which(!is.na(keep$distalIntergenic_closestTSS))), 'distal intergenic peaks assigned to closest TSS \n')
  
  saveRDS(keep, file = paste0(RdataDir, '/RAR_chipseq.peak.assignment_promoters.genebody.downstream.chiapet.closestTSS.rds'))
  write.csv(keep, file = paste0(resDir, '/RAR_chipseq.peak.assignment_promoters.genebody.downstream.chiapet.closestTSS.csv'), 
            row.names = TRUE, quote = FALSE)
  
}

########################################################
########################################################
# Section : RAR targets identification by combining peak targets and RNA-seq data
# 
########################################################
########################################################
Combine.RNAseq.peakAssignment = FALSE
if(Combine.RNAseq.peakAssignment){
  
  peaks = readRDS(file = paste0(RdataDir, '/RAR_chipseq.peak.assignment_promoters.genebody.downstream.chiapet.closestTSS.rds'))
  
  ggs = unlist(peaks[, c(13:17)])
  ggs = ggs[!is.na(ggs)]
  ggs = unique(ggs)
  xx = c()
  for(n in 1:length(ggs))
  {
    xx = c(xx, unlist(strsplit(as.character(ggs[n]), ';')))
  }
  ggs = unique(xx)
  
  targets = data.frame(matrix(NA, ncol = 5, nrow = length(ggs)))
  rownames(targets) = ggs
  colnames(targets) = colnames(peaks)[13:17]
  
  index.col = c(13:17)
  
  for(n in 1:nrow(targets))
  {
    cat(n, '\n')
    for(m in 1:length(index.col))
    {
      # n = 1; m = 1
      
      jjj = grep(rownames(targets)[n], peaks[, index.col[m]])
      if(length(jjj)>0){
        for(jj in jjj){
          test = unique(unlist(strsplit(as.character(peaks[jj, index.col[m]]), ';')))
          
          if(length(which(test == rownames(targets)[n])) > 0){
            if(!is.na(targets[n, m])) {
              targets[n, m] = paste0(c(targets[n, m], rownames(peaks)[jj]), collapse = ';')
            }else{
              targets[n, m] = rownames(peaks)[jj]
            }
          }
        }
      }
    }
  }
  
  res = read.csv(file = paste0(resDir, '/DESeq2_DEgenes_fdr.0.05_normalizedData.csv'), row.names = 1, header = TRUE)
  mm = match(rownames(res), rownames(targets))
  xx = data.frame(res[which(!is.na(mm)), ], targets[mm[which(!is.na(mm))], ], stringsAsFactors = FALSE)
  
  write.csv(xx, file = paste0(resDir, '/RAR_targetsLists_DEgeneRNAseq_peakAssignment.csv'), 
            row.names = TRUE, quote = FALSE)
  
  
  
}


