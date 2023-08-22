##########################################################################
##########################################################################
# Project:
# Script purpose: function from 
# https://github.com/juliendelile/MouseSpinalCordAtlas#2-knowledge-based-identification-of-all-cell-populations
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Sep 30 17:15:14 2022
##########################################################################
##########################################################################
#' We use exclusive markers to define partitionning population centers to which 
#' we associate the closest cells (no combinatorial genes)
doCellPartition <- function(literature_markers.df = refs,
                            readcounts=m$getReadcounts(data_status='Raw'),
                            cell_level_min_step1 = 2,
                            cell_level_min_step2 = 1,
                            normalize=TRUE)
{
  
  cat(paste0("Load population definitions \n"))
  
  # Load population definitions
  #literature_markers.df = read.table(file=known_template_file, sep='\t', stringsAsFactors=F, header=TRUE, skipNul=T)
  
  # Make list of genes
  literature_markers.df <- literature_markers.df %>%
    dplyr::mutate(
      gmap_step1=strsplit(gsub(" ", "", Genes_map_step1, fixed = TRUE), split=","),
      gmap_step2=strsplit(gsub(" ", "", Genes_map_step2, fixed = TRUE), split=","),
      Neural_pop_unique=make.unique(Neural_pop), # New for MN
      name=case_when(Type %in% c('Progenitor', 'Neuron') ~ Neural_pop, !Type %in% c('Progenitor', 'Neuron') ~ Type), # No need for rownames (?)
      name_unique=make.unique(name) # No need for rownames (?)
    )
  
  #' ### Step 1 / Progenitor vs Neuron vs extra populations
  cat(paste0("Step 1: Split macro populations (", paste0(unique(literature_markers.df$Type), collapse=", "), ") \n"))
  
  #' Create matrix of target population "centers"
  pop_def_mask_step1 = matrix(0,
                              nrow=length(unique(unlist(literature_markers.df$gmap_step1))), 
                              ncol=length(unique(literature_markers.df$Type)),
                              dimnames=list(
                                unique(unlist(literature_markers.df$gmap_step1)),
                                unique(literature_markers.df$Type)
                              )
  )
  
  for(gl in unique(literature_markers.df$Type)) {
    pop_def_mask_step1[
      literature_markers.df %>% dplyr::filter(Type==gl) %>% dplyr::pull(gmap_step1) %>% unlist
      , gl] <- 1
  }
  
  #' Filter lower UMI count for distance calculation (less than 2)
  
  #' Associate cells to the closest pre-defined center (Euclidean distance)
  cell_levels_step1 = 1*(readcounts[unique(unlist(literature_markers.df$gmap_step1)),] >= cell_level_min_step1)
  
  if(normalize){
    # cell_levels_step1 <- apply(cell_levels_step1, 2, function(x) x/sum(x))
    cell_levels_step1 <- apply(cell_levels_step1, 2, function(x) x/sqrt(sum(x^2)))
    cell_levels_step1[is.na(cell_levels_step1)] <- 0
    
    # pop_def_mask_step1 <- apply(pop_def_mask_step1, 2, function(x) x/sum(x))
    pop_def_mask_step1 <- apply(pop_def_mask_step1, 2, function(x) x/sqrt(sum(x^2)))
    pop_def_mask_step1[is.na(pop_def_mask_step1)] <- 0
  }
  
  cell_centers_alldist_step1 = as.matrix(pdist::pdist(t(cell_levels_step1), t(pop_def_mask_step1))) # euclidean distance, no other option in pdist
  cell_center_id_step1 = apply(cell_centers_alldist_step1, 1, which.min)
  
  #' Store associated type
  Type_step1 = factor(colnames(pop_def_mask_step1)[cell_center_id_step1], levels=unique(literature_markers.df$Type))
  
  #' ### Step 2 / Independent partition of progenitors and neurons
  
  cat("Step 2: Split neural populations \n")
  
  #' Create matrices of target progenitor and neuron population centers
  
  Type_step2_unique = rep("dummy", ncol(readcounts))
  
  pop_def_mask_step2 = sapply(c('Progenitor', 'Neuron'), function(pop){
    
    pop_markers = literature_markers.df %>%
      dplyr::filter(Type==pop) %>%
      {'names<-'(.$gmap_step2, .$Neural_pop_unique)} %>%
      {lapply(., function(x) x[x %in% rownames(readcounts)])}
    
    mask = matrix(0,
                  nrow=length(unique(unlist(pop_markers))), 
                  ncol=length(pop_markers),
                  # ncol=length(pop_markers)+1,
                  dimnames=list(
                    unique(unlist(pop_markers)),
                    names(pop_markers)
                    # c(names(pop_markers), paste0(pop, '_NULL'))
                  )
    )
    for(gl in names(pop_markers)) {
      mask[pop_markers[[gl]], gl] <- 1
    }
    mask
  })
  
  #' Associate cells to the closest pre-defined center (Euclidean distance)
  for(pop in c('Progenitor', 'Neuron')){
    pop_ids = which(Type_step1 == pop)
    cell_levels_step2 = 1*(readcounts[rownames(pop_def_mask_step2[[pop]]), pop_ids] >= cell_level_min_step2)
    
    if(normalize){
      # cell_levels_step2 <- apply(cell_levels_step2, 2, function(x) x/sum(x))
      cell_levels_step2 <- apply(cell_levels_step2, 2, function(x) x/sqrt(sum(x^2)))
      cell_levels_step2[is.na(cell_levels_step2)] <- 0
      
      # pop_def_mask_step2[[pop]] <- apply(pop_def_mask_step2[[pop]], 2, function(x) x/sum(x))
      pop_def_mask_step2[[pop]] <- apply(pop_def_mask_step2[[pop]], 2, function(x) x/sqrt(sum(x^2)))
      pop_def_mask_step2[[pop]][is.na(pop_def_mask_step2[[pop]])] <- 0
    }
    
    cell_centers_alldist_step2 = as.matrix(pdist::pdist(t(cell_levels_step2), t(pop_def_mask_step2[[pop]]))) # euclidean distance, no other option in pdist
    cell_center_id_step2 = apply(cell_centers_alldist_step2, 1, which.min)
    Type_step2_unique[pop_ids] <- colnames(pop_def_mask_step2[[pop]])[cell_center_id_step2]
  }
  
  # remaining "dummy" types are non-neural tissues
  Type_step2_unique[which(Type_step2_unique=="dummy")] <- as.character(Type_step1[which(Type_step2_unique=="dummy")])
  
  # Replace temparary neural pop names (for MN)
  popname_map = literature_markers.df %>% dplyr::filter(Type %in% c('Progenitor', 'Neuron')) %>% {"names<-"(.$Neural_pop, .$Neural_pop_unique)}
  Type_step2 <- popname_map[Type_step2_unique]
  
  Type_step2 <- factor(Type_step2, levels=unique(literature_markers.df$name))
  Type_step2_unique <- factor(Type_step2_unique, levels=unique(literature_markers.df$name_unique))
  
  # DV = (literature_markers.df %>% tibble::rownames_to_column("Pop") %>% {'names<-'(.$DV, .$Pop)})[as.character(Type_step2)]
  DV = (literature_markers.df %>% {'names<-'(.$DV, .$name)})[as.character(Type_step2)]
  DV[is.na(DV)] <- max(DV, na.rm=TRUE) + 1
  
  bothSteps_markers = setNames(
    apply(literature_markers.df, 1, function(x) c(x[['gmap_step1']], x[["gmap_step2"]])) %>% {lapply(., function(x) x[x %in% rownames(readcounts)])},
    literature_markers.df$name_unique)
  
  bothSteps_markers_neural_unique = bothSteps_markers[literature_markers.df %>% dplyr::filter(Type %in% c('Neuron', 'Progenitor')) %>% .$name_unique]
  
  bothSteps_markers_neural = literature_markers.df %>%
    dplyr::mutate(allmarkers=bothSteps_markers) %>%           # reuse merged gene list
    dplyr::filter(Type %in% c('Progenitor', 'Neuron')) %>%    # filter non neural
    dplyr::mutate(Neural_pop=factor(Neural_pop, levels=unique(Neural_pop))) %>% # factor to keep order
    dplyr::group_by(Neural_pop) %>%
    dplyr::summarise(gl=list(unique(unlist(allmarkers))), step2_length=length(unique(unlist(gmap_step2)))) %>%
    dplyr::filter(step2_length!=0) %>% # remove null neural populations
    {"names<-"(.$gl, .$Neural_pop)}
  
  return(list(
    Type_step1=Type_step1,
    Type_step2=Type_step2,
    Type_step2_unique=Type_step2_unique,
    DV=DV,
    step1_markers=literature_markers.df %>% dplyr::select(Type, gmap_step1) %>% dplyr::distinct(Type, .keep_all=T) %>% {"names<-"(.$gmap_step1, .$Type)},
    step2_markers=literature_markers.df %>% dplyr::filter(Type %in% c('Neuron', 'Progenitor')) %>% {"names<-"(.$gmap_step2, .$name_unique)},
    bothSteps_markers=bothSteps_markers,
    bothSteps_markers_neural=bothSteps_markers_neural,
    bothSteps_markers_neural_unique=bothSteps_markers_neural_unique
  ))
  
}