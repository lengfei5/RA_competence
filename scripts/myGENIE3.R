#' @title GENIE3
#' 
#' @description \code{GENIE3} Infers a gene regulatory network (in the form of a weighted adjacency matrix) from expression data, using ensembles of regression trees.
#'
#' @param expr.matrix Expression matrix (genes x samples). Every row is a gene, every column is a sample. 
#' @param tree.method Tree-based method used. Must be either "RF" for Random Forests (default) or "ET" for Extra-Trees.
#' @param K Number of candidate regulators randomly selected at each tree node (for the determination of the best split). Must be either "sqrt" for the square root of the total number of candidate regulators (default), "all" for the total number of candidate regulators, or a stricly positive integer.
#' @param ntrees Number of trees in an ensemble for each target gene. Default: 1000.
#' @param regulators Subset of genes used as candidate regulators. Must be either a vector of indices, e.g. \code{c(1,5,6,7)}, or a vector of gene names, e.g. \code{c("at_12377", "at_10912")}. The default value NULL means that all the genes are used as candidate regulators.
#' @param ncores Number of cores to use for parallel computing. Default: 1.
#' @param verbose If set to TRUE, a feedback on the progress of the calculations is given. Default: FALSE.
#' @param seed Random number generator seed for replication of analyses. The default value NULL means that the seed is not reset.
#'
#' @return Weighted adjacency matrix of inferred network. Element w_ij (row i, column j) gives the importance of the link from regulatory gene i to target gene j. 
#' 
#' @examples
#' ## Generate fake expression matrix
#' expr.matrix <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
#' rownames(expr.matrix) <- paste("Gene", 1:20, sep="")
#' colnames(expr.matrix) <- paste("Sample", 1:5, sep="")
#'
#' ## Run GENIE3
#' weight.matrix <- GENIE3(expr.matrix, regulators=paste("Gene", 1:5, sep=""))
#' 
#' ## Get ranking of edges 
#' link.list <- get.link.list(weight.matrix)
#' head(link.list)
#' @export
GENIE3 <- function(expr.matrix, priorRegulators = NULL, 
                   tree.method="RF", K="sqrt", ntrees=1000, regulators=NULL, ncores=1, 
                   verbose=TRUE, seed=NULL) 
{
  # expr.matrix = E;  priorRegulators = target.tfs; tree.method="RF"; K="sqrt"; ntrees=1000; regulators=NULL; ncores=1; 
  # verbose=TRUE; seed=NULL;
  dyn.load("GENIE3.so")
  
	# check input arguments
	if (!is.matrix(expr.matrix) && !is.array(expr.matrix)) {
		stop("Parameter expr.matrix must be a two-dimensional matrix where each row corresponds to a gene 
		     and each column corresponds to a condition/sample.")
	}
	
	if (length(dim(expr.matrix)) != 2) {
		stop("Parameter expr.matrix must be a two-dimensional matrix 
		     where each row corresponds to a gene and each column corresponds to a condition/sample.")
	}
	
	if (is.null(rownames(expr.matrix))) {
		stop("expr.matrix must specify the names of the genes in rownames(expr.matrix).")
	}
	
	if (tree.method != "RF" && tree.method != "ET") {
		stop("Parameter tree.method must be \"RF\" (Random Forests) or \"ET\" (Extra-Trees).")
	}
	
	if (K != "sqrt" && K != "all" && !is.numeric(K)) {
		stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
	}
	
	if (is.numeric(K) && K<1) {
		stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
	}
	
	if (!is.numeric(ntrees) || ntrees<1) {
		stop("Parameter ntrees should be a stricly positive integer.")
	}
	
	if (!is.null(regulators)) {
		if (!is.vector(regulators)) {
			stop("Parameter regulators must be either a vector of indices or a vector of gene names.")
		}
		
		if (is.character(regulators) && length(intersect(regulators,rownames(expr.matrix))) == 0) {
			stop("The genes must contain at least one candidate regulator.")
		}
		
		if (is.numeric(regulators) && max(regulators) > dim(expr.matrix)[1]) {
			stop("At least one index in regulators exceeds the number of genes.")
		}
	}
	
	if (!is.numeric(ncores) || ncores<1) {
		stop("Parameter ncores should be a stricly positive integer.")
	}
	
	# set random number generator seed if seed is given
  if (!is.null(seed)) {
     set.seed(seed)
  }
    
  # transpose expression matrix to (samples x genes)
  expr.matrix <- t(expr.matrix)
  
  # setup weight matrix
  num.samples <- dim(expr.matrix)[1]
  num.genes <- dim(expr.matrix)[2]
  gene.names <- colnames(expr.matrix)
  weight.matrix <- matrix(0.0, nrow=num.genes, ncol=num.genes)
  rownames(weight.matrix) <- gene.names
  colnames(weight.matrix) <- gene.names
  
  # get names of input genes
  if (is.null(regulators)) {
    input.gene.names <- gene.names
  } else {
    
    # input gene indices given as integers
    if (is.numeric(regulators)) {
      input.gene.names <- gene.names[regulators]
      # input gene indices given as names
    } else {
      input.gene.names <- regulators
      # for security, abort if some input gene name is not in gene names
      missing.gene.names <- setdiff(input.gene.names, gene.names)
      if (length(missing.gene.names) != 0) {
        for (missing.gene.name in missing.gene.names) {
          cat(paste("Gene ", missing.gene.name,
                    " was not in the expression matrix\n", sep=""))
        }
        stop("Aborting computation")
      }
    }
  }
  
  
	# tree method
	if (tree.method == 'RF') {
		RF_randomisation <- 1
		ET_randomisation <- 0
		bootstrap_sampling <- 1
	} else {
		RF_randomisation <- 0
		ET_randomisation <- 1
		bootstrap_sampling <- 0
	} 
	
	if (verbose) {
	  cat(paste("Tree method: ", tree.method, "\nK: ", K,
	            "\nNumber of trees: ", ntrees, "\n\n",
	            sep=""))
	  flush.console()
	  
	}
    
  # compute importances for every target gene
	if (ncores==1) {
		# serial computing
		if (verbose) {
		  cat("Using 1 core.\n\n")
		  flush.console()
		}
	  
	  for (target.gene.idx in seq(from=1, to=num.genes)) 
	  {
	    # target.gene.idx = 244; verbose = TRUE;
	    if (verbose) {	
	      cat(paste("Computing gene ", target.gene.idx, "/", num.genes, "\n", sep=""))
	      flush.console()
	    }
	    
	    target.gene.name <- gene.names[target.gene.idx]
	    
	    # remove target gene from input genes
	    these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
	    if(!is.null(priorRegulators)){
	      preSels = priorRegulators[which(rownames(priorRegulators) == target.gene.name), ]
	      preSels_regulators = colnames(priorRegulators)[which(preSels>0)]
	      these.input.gene.names = intersect(these.input.gene.names, preSels_regulators)
	      #cat('nb of prior regulators : ', length(these.input.gene.names), '\n')
	    }
	    
	    num.input.genes <- length(these.input.gene.names)
	    
	    x <- expr.matrix[ ,these.input.gene.names]
	    y <- expr.matrix[ ,target.gene.name]
	    
	    # normalize output data
	    y <- y / sd(y)
	    
	    # set mtry
	    if (class(K) == "numeric") {
	      mtry <- K
	    } else if (K == "sqrt") {
	      mtry <- round(sqrt(num.input.genes))
	    } else {
	      mtry <- num.input.genes
	    } 
	    
	    # some default parameters 
	    nmin <- 1
	    permutation_importance <- 0
	    
	    im <- .C("BuildTreeEns",
	             as.integer(num.samples),
	             as.integer(num.input.genes),
	             as.single(c(x)),
	             as.single(c(y)),
	             as.integer(nmin),
	             as.integer(ET_randomisation),
	             as.integer(RF_randomisation),
	             as.integer(mtry),
	             as.integer(ntrees),
	             as.integer(bootstrap_sampling),
	             as.integer(permutation_importance),
	             as.double(vector("double",num.input.genes)))[[12]]
	    
	    # some variable importances might be slighly negative due to some rounding error
	    im[im<0] <- 0
	    weight.matrix[these.input.gene.names, target.gene.name] <- im
	    
	  }
	  
	}else{
		# parallel computing
	  library(doRNG); 
	  library(doParallel); 
	  registerDoParallel(); 
	  options(cores=ncores)
		
		if (verbose) {
		  message(paste("\nUsing", getDoParWorkers(), "cores."))
		}
		
	  weight.matrix.reg <- foreach(target.gene.name=gene.names, .combine=cbind) %dorng% 
	  {
	    # remove target gene from input genes
	    these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
	    if(!is.null(priorRegulators)){
	      preSels = priorRegulators[which(rownames(priorRegulators) == target.gene.name), ]
	      preSels_regulators = colnames(priorRegulators)[which(preSels>0)]
	      these.input.gene.names = intersect(these.input.gene.names, preSels_regulators)
	      #cat('nb of prior regulators : ', length(these.input.gene.names), '\n')
	    }
	    
	    num.input.genes <- length(these.input.gene.names)
	    
	    x <- expr.matrix[ ,these.input.gene.names]
	    y <- expr.matrix[ ,target.gene.name]
	    
	    # normalize output data
	    y <- y / sd(y)
	    
	    # set mtry
	    if (class(K) == "numeric") {
	      mtry <- K
	    } else if (K == "sqrt") {
	      mtry <- round(sqrt(num.input.genes))
	    } else {
	      mtry <- num.input.genes
	    } 
	    
	    # some default parameters 
	    nmin <- 1
	    permutation_importance <- 0
	    
	    im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(num.input.genes),
	             as.single(c(x)),as.single(c(y)),as.integer(nmin),
	             as.integer(ET_randomisation),as.integer(RF_randomisation),
	             as.integer(mtry),as.integer(ntrees),
	             as.integer(bootstrap_sampling),as.integer(permutation_importance),
	             as.double(vector("double",num.input.genes)))[[12]]
	    
	    # some variable importances might be slighly negative due to some rounding error
	    im[im<0] <- 0
	    c(setNames(0, target.gene.name), setNames(im, these.input.gene.names))[input.gene.names]
	        
	  }
	  
	  attr(weight.matrix.reg, "rng") <- NULL
	  weight.matrix[input.gene.names,] <- weight.matrix.reg
	    
	}
  
  return(weight.matrix / num.samples)
  
}


GENIE3_withSpecificTargets <- function(expr.matrix, 
                                       priorTargets = c('Foxa2'),
                                       tree.method="RF", 
                                       K="sqrt", 
                                       ntrees=1000, regulators=NULL, ncores=1, verbose=FALSE, 
                                       seed=NULL) {
  
  # expr.matrix = E; priorTargets = c('Foxa2'); tree.method="RF"; K="sqrt";ntrees=1000;regulators=NULL; ncores=1
  dyn.load("GENIE3.so")
  
  # check input arguments
  if (!is.matrix(expr.matrix) && !is.array(expr.matrix)) {
    stop("Parameter expr.matrix must be a two-dimensional matrix where 
         each row corresponds to a gene and each column corresponds to a condition/sample.")
  }
  
  if (length(dim(expr.matrix)) != 2) {
    stop("Parameter expr.matrix must be a two-dimensional 
         matrix where each row corresponds to a gene and each column corresponds to a condition/sample.")
  }
  
  if (is.null(rownames(expr.matrix))) {
    stop("expr.matrix must specify the names of the genes in rownames(expr.matrix).")
  }
  
  if (tree.method != "RF" && tree.method != "ET") {
    stop("Parameter tree.method must be \"RF\" (Random Forests) or \"ET\" (Extra-Trees).")
  }
  
  if (K != "sqrt" && K != "all" && !is.numeric(K)) {
    stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
  }
  
  if (is.numeric(K) && K<1) {
    stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
  }
  
  if (!is.numeric(ntrees) || ntrees<1) {
    stop("Parameter ntrees should be a stricly positive integer.")
  }
  
  if (!is.null(regulators)) {
    if (!is.vector(regulators)) {
      stop("Parameter regulators must be either a vector of indices or a vector of gene names.")
    }
    
    if (is.character(regulators) && length(intersect(regulators,rownames(expr.matrix))) == 0) {
      stop("The genes must contain at least one candidate regulator.")
    }
    
    if (is.numeric(regulators) && max(regulators) > dim(expr.matrix)[1]) {
      stop("At least one index in regulators exceeds the number of genes.")
    }
  }
  
  if (!is.numeric(ncores) || ncores<1) {
    stop("Parameter ncores should be a stricly positive integer.")
  }
  
  
  # set random number generator seed if seed is given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # transpose expression matrix to (samples x genes)
  expr.matrix <- t(expr.matrix)
  
  # setup weight matrix
  num.samples <- dim(expr.matrix)[1]
  num.genes <- dim(expr.matrix)[2]
  gene.names <- colnames(expr.matrix)
  
  weight.matrix <- matrix(0.0, nrow=num.genes, ncol=length(priorTargets))
  rownames(weight.matrix) <- gene.names
  colnames(weight.matrix) <- priorTargets
  
  # get names of input genes
  if (is.null(regulators)) {
    input.gene.names <- gene.names
  } else {
    # input gene indices given as integers
    if (is.numeric(regulators)) {
      input.gene.names <- gene.names[regulators]
      # input gene indices given as names
    } else {
      input.gene.names <- regulators
      # for security, abort if some input gene name is not in gene names
      missing.gene.names <- setdiff(input.gene.names, gene.names)
      if (length(missing.gene.names) != 0) {
        for (missing.gene.name in missing.gene.names) {
          cat(paste("Gene ", missing.gene.name,
                    " was not in the expression matrix\n", sep=""))
        }
        stop("Aborting computation")
      }
    }
  }
  
  # tree method
  if (tree.method == 'RF') {
    RF_randomisation <- 1
    ET_randomisation <- 0
    bootstrap_sampling <- 1
  } else {
    RF_randomisation <- 0
    ET_randomisation <- 1
    bootstrap_sampling <- 0
  } 
  
  if (verbose) {
    cat(paste("Tree method: ", tree.method, "\nK: ", K,
              "\nNumber of trees: ", ntrees, "\n\n",
              sep=""))
    flush.console()
  }
  
  # compute importances for every target gene
  
  if (ncores==1) {
    # serial computing
    if (verbose) {
      cat("Using 1 core.\n\n")
      flush.console()
    }
    
    for (target.gene.idx in seq(from=1, to=length(priorTargets))) {
      
      # target.gene.idx = 1
      if (verbose) {	
        cat(paste("Computing gene ", target.gene.idx, "/", length(priorTargets), "\n", sep=""))
        flush.console()
      }
      
      target.gene.name <- priorTargets[target.gene.idx]
      # remove target gene from input genes
      these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
      num.input.genes <- length(these.input.gene.names)
      
      x <- expr.matrix[,these.input.gene.names]
      y <- expr.matrix[,target.gene.name]
      
      # normalize output data
      y <- y / sd(y)
      
      # set mtry
      if (class(K) == "numeric") {
        mtry <- K
      } else if (K == "sqrt") {
        mtry <- round(sqrt(num.input.genes))
      } else {
        mtry <- num.input.genes
      } 
      
      # some default parameters 
      nmin <- 1
      permutation_importance <- 0
      
      im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(num.input.genes),
               as.single(c(x)),as.single(c(y)),as.integer(nmin),
               as.integer(ET_randomisation),as.integer(RF_randomisation),
               as.integer(mtry),as.integer(ntrees),
               as.integer(bootstrap_sampling),as.integer(permutation_importance),
               as.double(vector("double",num.input.genes)))[[12]]
      
      # some variable importances might be slighly negative due to some rounding error
      im[im<0] <- 0
      
      weight.matrix[these.input.gene.names, target.gene.name] <- im
    }
  } else {
    # parallel computing
    library(doRNG); library(doParallel); registerDoParallel(); options(cores=ncores)
    
    if (verbose) {
      message(paste("\nUsing", getDoParWorkers(), "cores."))
    }
    
    weight.matrix.reg <- foreach(target.gene.name=gene.names, .combine=cbind) %dorng% 
      {
        # remove target gene from input genes
        these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
        num.input.genes <- length(these.input.gene.names)
        
        x <- expr.matrix[,these.input.gene.names]
        y <- expr.matrix[,target.gene.name]
        
        # normalize output data
        y <- y / sd(y)
        
        # set mtry
        if (class(K) == "numeric") {
          mtry <- K
        } else if (K == "sqrt") {
          mtry <- round(sqrt(num.input.genes))
        } else {
          mtry <- num.input.genes
        } 
        
        # some default parameters 
        nmin <- 1
        permutation_importance <- 0
        
        im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(num.input.genes),
                 as.single(c(x)),as.single(c(y)),as.integer(nmin),
                 as.integer(ET_randomisation),as.integer(RF_randomisation),
                 as.integer(mtry),as.integer(ntrees),
                 as.integer(bootstrap_sampling),as.integer(permutation_importance),
                 as.double(vector("double",num.input.genes)))[[12]]
        
        # some variable importances might be slighly negative due to some rounding error
        im[im<0] <- 0
        
        c(setNames(0, target.gene.name), setNames(im, these.input.gene.names))[input.gene.names]
      }
    attr(weight.matrix.reg, "rng") <- NULL
    weight.matrix[input.gene.names,] <- weight.matrix.reg
  }
  return(weight.matrix / num.samples)
}       




#' @title get.link.list
#' 
#' @description \code{get.link.list} Converts the weight matrix, as returned by \code{\link{GENIE3}}, to a sorted list of regulatory links (most likely links first).
#' 
#' @param weight.matrix Weighted adjacency matrix as returned by \code{\link{GENIE3}}.
#' @param report.max Maximum number of links to report. The default value NULL means that all the links are reported.
#' @param threshold Only links with a weight above the threshold are reported. Default: threshold = 0, i.e. all the links are reported.
#' 
#' @return List of regulatory links in a data frame. Each line of the data frame corresponds to a link. The first column is the regulatory gene, the second column is the target gene, and the third column is the weight of the link.
#'
#' @seealso \code{\link{GENIE3}}
#'
#' @examples
#' ## Generate fake expression matrix
#' expr.matrix <- matrix(sample(1:10, 100, replace=TRUE), nrow=20)
#' rownames(expr.matrix) <- paste("Gene", 1:20, sep="")
#' colnames(expr.matrix) <- paste("Sample", 1:5, sep="")
#'
#' ## Run GENIE3
#' weight.matrix <- GENIE3(expr.matrix, regulators=paste("Gene", 1:5, sep=""))
#' 
#' ## Get ranking of edges 
#' link.list <- get.link.list(weight.matrix)
#' head(link.list)
#' @export
get.link.list <- function(weight.matrix, report.max=NULL, threshold=0) {
  if(!is.numeric(threshold)) {
    stop("threshold must be a number.")
  } 
	
	library(reshape2)
	
	# Only process weights off-diagonal
	diag(weight.matrix) <- NA
    link.list <- melt(weight.matrix, na.rm=TRUE)
    colnames(link.list) <- c("regulatory.gene", "target.gene", "weight")
    link.list <- link.list[link.list$weight>=threshold,]
    link.list <- link.list[order(link.list$weight, decreasing=TRUE),]
  
    if(!is.null(report.max)) {
    	link.list <- link.list[1:min(nrow(link.list), report.max),]
    } 
    
    rownames(link.list) <- NULL
    
    return(link.list)
}


read.expr.matrix <- function(filename, form="", sep="", default.gene.label="gene_", default.sample.label="sample_") {
    # Report when form is not correctly set
    if (form != "rows.are.genes" && form != "rows.are.samples") {
        stop("Parameter form must be set to \"rows.are.genes\" or \"rows.are.samples\"")
    }
    # read data
    m <- read.table(filename, sep=sep, as.is=TRUE)
    has.row.names <- FALSE
    has.col.names <- FALSE
    # have row and column names been recognized (cell (1,1) is empty in file) ?
    if (colnames(m)[1] != "V1" && rownames(m)[1] != "1") {
        has.row.names <- TRUE
        has.col.names <- TRUE
    }
    # is first column alphanumeric ?
    if (all(grepl("[[:alpha:]]", m[,1]))) {
        # check duplicated names
        if (any(duplicated(m[,1]))) {
            stop("Duplicated names in first column\n")
        }
        rownames(m) <- m[,1]
        m <- m[,-1]
        has.row.names <- TRUE
    }
    # is first row alphanumeric ?
    if (all(grepl("[[:alpha:]]", m[1,]))) {
        # check duplicated names
        if (any(duplicated(m[1,]))) {
            stop("Duplicated names in first row\n")
        }
        colnames(m) <- m[1,]
        m <- m[-1,]
        has.col.names <- TRUE
    }
    # coerce matrix data to numeric
    col.names <- colnames(m)
    row.names <- rownames(m)
    m <- as.matrix(m)
    m <- apply(m, 2, function(x) { as.numeric(x) })
    colnames(m) <- col.names
    rownames(m) <- row.names
    num.rows <- dim(m)[1]
    num.cols <- dim(m)[2]
    # fill in default gene names in rows if needed
    if (!has.row.names && form=="rows.are.genes") {
        rownames(m) <- paste(default.gene.label, seq(from=1, to=num.rows), sep="")
    }
    # fill in default sample names in rows if needed
    if (!has.row.names && form=="rows.are.samples") {
        rownames(m) <- paste(default.sample.label, seq(from=1, to=num.rows), sep="")
    }
    # fill in default sample names in columns if needed
    if (!has.col.names && form=="rows.are.genes") {
        colnames(m) <- paste(default.sample.label, seq(from=1, to=num.cols), sep="")
    }
    # fill in default gene names in columns if needed
    if (!has.col.names && form=="rows.are.samples") {
        colnames(m) <- paste(default.gene.label, seq(from=1, to=num.cols), sep="")
    }
    # transpose matrix to (genes x samples) if needed
    if (form == "rows.are.samples") m <- t(m)
    return(m)
}


##########################################
# construct the target-regulatory region matrix for target genes 
##########################################
build_target_CRE_matrix = function(targets.ids)
{
  # targets.ids = target.ids_all
  bed = read.table(file = 
                     paste0('/Volumes/groups/tanaka/People/current/jiwang/projects/positional_memory/motif_analysis/', 
                            'FIMO_atacPeak_tss_mediumQ_core.unvalided_64Gmem/peaks_for_fimo_sorted.bed'), header = FALSE)
  
  kk = match(bed$V4, targets.ids)
  mm = which(!is.na(kk))
  targets = data.frame(geneID = bed$V4[mm], CREs = paste0(bed$V1[mm], ':', bed$V2[mm], '-', bed$V3[mm]), stringsAsFactors = FALSE)
  
  enhancers = readRDS(file = paste0(RdataDir, '/enhancers_candidates_55k_atacPeaks_histM_H3K4me1_chipseekerAnnot_manual.rds'))
  mm = match(enhancers$targets, targets.ids)
  mm = which(!is.na(mm))
  
  targets2 = data.frame(geneID = enhancers$targets[mm], CREs = rownames(enhancers)[mm], stringsAsFactors = FALSE)
  
  # shift one bp due to the bed format and fasta formation
  targets2$CREs = sapply(targets2$CREs, function(x) {x = unlist(strsplit(gsub(':', '-', as.character(x)), '-'));
  paste0(x[1], ':', (as.numeric(x[2]) -1), '-', x[3])}) 
  
  targets = rbind(targets, targets2)
  rm(targets2)
  
  targets = targets[match(unique(targets$CREs), targets$CREs), ]
  targets$geneID = as.character(droplevels(targets$geneID))
  #targets$CREs =droplevels(targets$CREs)
  
  ## target-CRE matrix: each row is tf and each column is associated tss and enhancer peaks
  xx = table(targets$geneID, targets$CREs)
  ss = apply(xx, 1, sum)
  xx = xx[which(ss>0), ]
  
  mm = match(rownames(xx), targets.ids)
  # rownames(xx) = targets.ids
  
  return(xx)
  
}


##########################################
# plot the GRN graph 
# some original code from https://mr.schochastics.net/material/netvizr/
# https://github.com/quadbiolab/Pando/blob/main/R/plots.R (Pando from Jonas)
# Some basic network analysis 
# from https://github.com/Arizonagong/vCIES2020_Network-Analysis/blob/master/igraph/vCIES_igraph.Rscript.R
##########################################
plot_tf_network = function(link.list)
{
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(graphlayouts) 
  library(RColorBrewer) # This is the color library
  # require(uwot)
  
  trn <- graph_from_data_frame(link.list, directed = TRUE)
  
  #1. Igraph object summary
  gorder(trn) # nb of node
  gsize(trn) # nb of edges
  
  #2. Nodelist
  V(trn)
  
  #3. Edgelist
  E(trn)
  
  #4. Attributes
  V(trn)$name
  V(trn)$weight
  # compute a clustering for node colors
  V(trn)$module <- as.character(membership(cluster_louvain(graph_from_data_frame(link.list, directed = FALSE))))
  # compute degree as node size
  V(trn)$size <- degree(trn)
  
  #5. Adjacency matrix
  trn[c(1:10),c(1:10)]
  
  ##########################################
  # Measuring Centrality and add them into attributes
  ##########################################
  #1. Degree centrality
  V(trn)$degree<-degree(trn, mode = 'All')
  V(trn)$degreeIn = degree(trn, mode = 'in')
  V(trn)$degreeOut = degree(trn, mode = 'out')
  V(trn)$name[which.max(V(trn)$degree)]
  V(trn)$name[which.max(V(trn)$degreeIn)]
  V(trn)$name[which.max(V(trn)$degreeOut)]
  
  V(trn)$degree[which(V(trn)$name == 'ZNF281')]
    
  #2. Eigenvector centrality
  V(trn)$Eigen<-evcent(trn)$vector
  V(trn)$Eigen
  V(trn)$name[which.max(V(trn)$Eigen)]
  
  #3. Betweenness centrality
  V(trn)$betweenness <- betweenness(trn, directed = FALSE)
  V(trn)$betweenness
  V(trn)$name[which.max(V(trn)$betweenness)]
  
  DF <- as_long_data_frame(trn)
  trn
  
  # define a custom color palette
  V(trn)$cluster <- as.character(membership(cluster_louvain(graph_from_data_frame(link.list, directed = FALSE), 
                                                            resolution = 1)))
  V(trn)$module = V(trn)$cluster
  
  nb_clusters = length(unique(V(trn)$cluster))
  cat(nb_clusters, ' clusters used here \n')
  
  if(length(nb_clusters)<=9) {
    library(khroma)
    muted <- colour("muted")
    got_palette = c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99", 
                    "#999933", "#AA4499")[1:nb_clusters] 
  }else{
    cat('More than 9 colors needed \n')
    
    got_palette = c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99", 
                    "#1A5878", "#C44237", "#AD8941", "#E99093",
                    "#50594B", "#8968CD", "#9ACD32")[1:nb_clusters]
  }
 
  pal<-brewer.pal(nb_clusters, "Set3") # Vertex color assigned per each class number
  
  ###################
  ## test umap layout
  ###################
  kk = match(gnames$gnames[match(V(trn)$name, gnames$node)], rownames(E))
  matE = E[kk, ]
  pcs <- prcomp((matE), scale = TRUE)
  pcs = data.frame(pcs$x)
  
  kk2 = match(V(trn)$name, rownames(wtm))
  matC = wtm[kk2, ]
  matC[is.na(matC)] = 0
  pcs2 <- prcomp((matC), scale = TRUE)
  pcs2 = data.frame(pcs2$x)
  
  nb_pcs = 10; alpha = 0.;
  dist = dist(pcs[, c(1:nb_pcs)])
  dist2 = dist(pcs2[, c(1:nb_pcs)])
  
  dist_weighted = dist + alpha*dist2
  #pcs.use = data.frame(pcs[, c(1:10)], pcs2[, c(1:5)])
  
  set.seed(7)
  #weight_coex_umap <- uwot::umap(dist_weighted, n_neighbors=5)
  nb_pcs = 10;
  weight_coex_umap <- uwot::umap(pcs[, c(1:nb_pcs)], n_neighbors=5, min_dist = 0.01)
  #weight_coex_umap <- uwot::umap(, n_neighbors=5)
  rownames(weight_coex_umap) <- V(trn)$name
  colnames(weight_coex_umap) <- c('UMAP1', 'UMAP2')
  
  V(trn)$umap1 <-weight_coex_umap[,1]
  V(trn)$umap2 <-weight_coex_umap[,2]
  
  saveRDS(trn, file = paste0(RdataDir, '/TRN_umap_layout.rds'))
  
  trn = readRDS(file = paste0(RdataDir, '/TRN_umap_layout.rds'))
  nb_clusters = length(unique(V(trn)$cluster))
  cat(nb_clusters, ' clusters used here \n')
  got_palette = c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99", "#999999")
  
  #set.seed(2022)
  ggraph(trn, x=umap1, y=umap2) +
    geom_edge_link(aes(edge_width = weight), edge_colour = "gray85" ) +
    #geom_edge_diagonal(color='darkgray', width=0.2) +
    geom_node_point(aes(fill = module, size = degreeOut), shape = 21) +
    scale_size_continuous(range = c(1, 10)) +
    #geom_node_text(aes(filter = size >= 20, label = name), family = "serif") +
    geom_node_text(aes(filter = size >= 20, label=name), size=5/ggplot2::.pt, repel=T, family = "serif")+
    scale_edge_width_continuous(range = c(0.05, 0.2)) +
    scale_edge_color_gradientn(colors=rev(brewer.pal(n=11, name="RdBu")), limits=c(0.01, 0.6)) +
    #scale_edge_alpha_continuous(range=c(0.05, 0.4), limits=c(2,20)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) + 
    scale_fill_manual(values = got_palette) +
    #scale_fill_viridis(option='magma') + 
    #scale_fill_manual(values = pal) +
    coord_fixed() +
    theme_graph() +
    #theme(legend.position = "bottom") 
    theme(legend.position = "none") 
  ggsave(paste0(figureDir, "TRN_tfExpr.umap_connectivity.clusters_v7.pdf"), width=10, height = 10)
  
  p1 = ggraph(trn, x=umap1, y=umap2) +
    geom_edge_link(aes(edge_width = weight), edge_colour = "gray80" ) +
    #geom_edge_diagonal(color='darkgray', width=0.2) +
    geom_node_point(aes(fill = module, size = degreeOut), shape = 21) +
    scale_size_continuous(range = c(1, 8)) +
    #geom_node_text(aes(filter = size >= 20, label = name), family = "serif") +
    geom_node_text(aes(filter = size >= 20, label=name), size=6/ggplot2::.pt, repel=T, family = "serif")+
    scale_edge_width_continuous(range = c(0.05, 0.2)) +
    scale_edge_color_gradientn(colors=rev(brewer.pal(n=11, name="RdBu")), limits=c(0.01, 0.6)) +
    #scale_edge_alpha_continuous(range=c(0.05, 0.4), limits=c(2,20)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) + 
    scale_fill_manual(values = got_palette) +
    #scale_fill_viridis(option='magma') + 
    #scale_fill_manual(values = pal) +
    coord_fixed() + 
    theme_void() +
    theme(legend.position = "right") 
    #theme(legend.position = "none") 
  
  #ggsave(paste0(figureDir, "TRN_tfExpr.umap_connectivity.clusters_withLegends_v6.pdf"), width=10, height = 10)
  pdfname = paste0(figureDir,  "TRN_tfExpr.umap_connectivity.clusters_withLegends_v6.pdf")
  pdf(pdfname, width = 10, height = 10)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  p1
  #grid.arrange(p1, ncol = 1, nrow = 1)
  
  dev.off()
  
  
  ##########################################
  # summary of centrality
  ##########################################
  library(viridis)
  xx = data.frame(degree = V(trn)$degree, gene = V(trn)$name)
  xx = xx[order(-xx$degree), ]
  xx = xx[c(1:50), ]
  xx$gene = sapply(xx$gene, function(x) unlist(strsplit(as.character(x), '_'))[1])
  as_tibble(xx) %>% 
    ggplot(aes(y=degree, x=reorder(gene, -degree), fill = reorder(gene, -degree))) + 
    geom_bar(position="dodge", stat="identity") +
    theme_classic() +
    #theme(axis.text.x = element_text(angle = 90, size = 10)) +
    scale_fill_viridis_d(option = 'magma', direction = -1) +
    labs(x = '', y = 'number of connections') +
    theme(axis.text.x = element_text(angle = 60, size = 12, hjust = 1), 
          axis.text.y = element_text(angle = 0, size = 12), 
          axis.title =  element_text(size = 12),
          legend.text = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.position='none',
          #plot.margin = margin()
          #legend.key.size = unit(1, 'cm')
          #legend.key.width= unit(1, 'cm')
    )
  ggsave(paste0(figureDir, "GRN_centrality_all.pdf"),  width = 10, height = 6)
  
  xx = data.frame(degree = V(trn)$degreeOut, gene = V(trn)$name)
  xx = xx[order(-xx$degree), ]
  xx = xx[c(1:50), ]
  xx$gene = sapply(xx$gene, function(x) unlist(strsplit(as.character(x), '_'))[1])
  as_tibble(xx) %>% 
    ggplot(aes(y=degree, x=reorder(gene, -degree), fill = reorder(gene, -degree))) + 
    geom_bar(position="dodge", stat="identity") +
    theme_classic() +
    #theme(axis.text.x = element_text(angle = 90, size = 10)) +
    scale_fill_viridis_d(option = 'magma', direction = -1) +
    labs(x = '', y = 'centrality') +
    theme(axis.text.x = element_text(angle = 90, size = 14), 
          axis.text.y = element_text(angle = 0, size = 14), 
          axis.title =  element_text(size = 14),
          legend.text = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.position='none',
          #plot.margin = margin()
          #legend.key.size = unit(1, 'cm')
          #legend.key.width= unit(1, 'cm')
    )
  ggsave(paste0(figureDir, "GRN_centrality_outDegree.pdf"),  width = 10, height = 4)
  
  
  #plot(p1)
  # pdf(paste0(figureDir, "/GRN_UMPA_Test_v3.pdf"),  width = 14, height = 10) # Open a new pdf file
  # for(n in c(1:25, 1011, 2022))
  # {
  #   cat(n, '\n')
  # }
  # dev.off()
  
  # # basic graph
  # ggraph(trn, layout = "stress") +
  #   geom_edge_link0(aes(edge_width = weight), edge_colour = "grey66") +
  #   geom_node_point(aes(fill = cluster, size = size), shape = 21) +
  #   geom_node_text(aes(filter = size >= 20, label = name), family = "serif") +
  #   scale_fill_manual(values = c(got_palette, 'red', 'blue')) +
  #   scale_edge_width(range = c(0.2, 3)) +
  #   scale_size(range = c(1, 6)) +
  #   theme_graph() +
  #   theme(legend.position = "none")

  # trn.undirected = graph_from_data_frame(link.list, directed = FALSE)
  # bb <- layout_as_backbone(trn.undirected,keep=0.4)
  # E(g)$col <- F
  # E(g)$col[bb$backbone] <- T
  
  # ggraph(g,layout="manual",x=bb$xy[,1],y=bb$xy[,2])+
  #   geom_edge_link0(aes(col=col),width=0.1)+
  #   geom_node_point(aes(col=grp))+
  #   scale_color_brewer(palette = "Set1")+
  #   scale_edge_color_manual(values=c(rgb(0,0,0,0.3),rgb(0,0,0,1)))+
  #   theme_graph()+
  #   theme(legend.position = "none")
  
  # centrality layout
  # https://github.com/schochastics/graphlayouts
  
  # set.seed(2022)
  # ggraph(trn, layout = "centrality", cent = graph.strength(trn)) +
  #   geom_edge_link0(aes(edge_width = weight), edge_colour = "grey66") +
  #   geom_node_point(aes(fill = cluster, size = size), shape = 21) +
  #   geom_node_text(aes(filter = size >= 20, label = name), family = "serif") +
  #   scale_edge_width_continuous(range = c(0.02, 0.1)) +
  #   scale_size_continuous(range = c(1, 10)) +
  #   scale_fill_manual(values = got_palette) +
  #   coord_fixed() +
  #   theme_graph() +
  #   #theme(legend.position = "bottom")
  #   theme(legend.position = "none")
  # 
  # ggsave(paste0(resDir, "/TRN_secondTest.pdf"), width=12, height = 10)
  # 
  ######################################################################
  ## following code from https://github.com/lengfei5/pallium_evo/blob/main/analysis/GRN_analysis/moo_graph_layout.R
  ######################################################################
  #### Net with only TFs ####
  # #### Get avg expression for genes ####
  # region_summary <- Pando::aggregate_matrix(rna_expr[, union(grn_net$tf, grn_net$target)], groups=mome_atac$pred_regions_all)
  # subclass_summary <- Pando::aggregate_matrix(rna_expr[, union(grn_net$tf, grn_net$target)], groups=mome_atac$subclasses)
  # 
  # region_summary_df <- region_summary %>% t() %>% 
  #   as_tibble(rownames='gene') 
  # 
  # subclass_summary_df <- subclass_summary %>% t() %>% 
  #   as_tibble(rownames='gene') 
  # 
  # gene_scores <- inner_join(region_summary_df, subclass_summary_df)
  # 
  # ### Get coex and umap ####
  # gene_cor <- Pando::sparse_cor(rna_expr[, union(grn_net$tf, grn_net$target)])
  # 
  # reg_mat <- grn_net %>% 
  #   filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
  #   filter(target%in%tfs$symbol) %>%
  #   distinct(target, tf, estimate) %>%
  #   pivot_wider(names_from=tf, values_from=estimate, values_fill=0) %>% 
  #   column_to_rownames('target') %>% as.matrix() %>% Matrix::Matrix(sparse=T)
  # reg_factor_mat <- abs(reg_mat) + 1
  # 
  # weighted_coex_mat <- gene_cor[rownames(reg_factor_mat), colnames(reg_factor_mat)] * sqrt(reg_factor_mat)
  # weighted_coex_mat <- as.matrix(gene_cor[rownames(reg_factor_mat), colnames(reg_factor_mat)])
  # weight_coex_umap <- uwot::umap(weighted_coex_mat, n_neighbors=5)
  # rownames(weight_coex_umap) <- rownames(weighted_coex_mat)
  # colnames(weight_coex_umap) <- c('UMAP1', 'UMAP2')
  # 
  # 
  # #### Plot network ####
  # weight_coex_meta <- weight_coex_umap %>% 
  #   as_tibble(rownames='gene') %>% 
  #   left_join(gene_scores)
  # 
  # tf_graph <- as_tbl_graph(grn_net) %>% 
  #   activate(edges) %>% 
  #   mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>% 
  #   activate(nodes) %>% 
  #   mutate(
  #     central_pr=centrality_pagerank(weights = estimate),
  #     central_betw=centrality_betweenness(),
  #     central_eig=centrality_eigen(),
  #     central_deg=centrality_degree(),
  #     outdegree=centrality_degree(mode='out'),
  #     indegree=centrality_degree(mode='in')
  #   ) %>% 
  #   inner_join(weight_coex_meta, by=c('name'='gene')) %>% 
  #   activate(edges) %>%
  #   filter(padj<0.05) %N>% 
  #   filter(!node_is_isolated())
  # 
  # ggraph(tf_graph, x=UMAP1, y=UMAP2) + 
  #   geom_edge_diagonal(aes(alpha=-log10(padj), color=factor(sign(estimate))), width=0.5) + 
  #   geom_node_point(aes(size=outdegree), shape=21, color='black', fill='grey') +
  #   geom_node_text(aes(label=name), size=5/ggplot2::.pt, repel=T) +
  #   scale_edge_color_manual(values=c('#f5b7b1', '#7dcea0')) +
  #   scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
  #   theme_void() 
  # ggsave('plots/tf_grn_umap.png', width=6, height=6)
  # ggsave('plots/tf_grn_umap.pdf', width=6, height=6)
  # 
  # 
  # ggraph(tf_graph, x=UMAP1, y=UMAP2) + 
  #   geom_edge_diagonal(aes(alpha=-log10(padj), color=factor(sign(estimate))), width=0.5) + 
  #   geom_node_point(aes(size=outdegree), shape=21, color='black', fill='grey') +
  #   # geom_node_text(aes(label=name), size=5/ggplot2::.pt, repel=T) +
  #   scale_edge_color_manual(values=c('#f5b7b1', '#7dcea0')) +
  #   scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
  #   theme_void() 
  # ggsave('plots/tf_grn_unlabelled_umap.png', width=6, height=6)
  # ggsave('plots/tf_grn_unlabelled_umap.pdf', width=6, height=6)
  
}

build_subgraph_GRN = function(pp, subgroup = 'BL.day5')
{
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(graphlayouts) 
  library(RColorBrewer) 
  
  trn = readRDS(file = paste0(RdataDir, '/TRN_umap_layout.rds'))
  
  # saved matrix: E, mat_gcre, mocs, atm, target.tfs,
  load(file = paste0(RdataDir, '/target_CREs_motifs_tfs_matrix_for_targetRegulator.Rdata')) 
  target.tfs =readRDS(file = paste0(RdataDir, '/GRN_priorNetwork_target_TFs.rds'))
  #target.tfs =readRDS(file = paste0(RdataDir, '/GRN_priorNetwork_target_TFs_final.rds'))
  geneName_mapping = readRDS(file = paste0(RdataDir, '/geneName_mapping.rds'))
  
  mat_gcre_sub = mat_gcre
  mm = match(colnames(mat_gcre), pp)
  jj = which(is.na(mm))
  for(j in jj) mat_gcre_sub[ ,j] = 0
  
  target.tfs_sub = mat_gcre_sub %*% mocs %*% atm
  
  target.tfs_sub = target.tfs_sub[match(rownames(target.tfs), rownames(target.tfs_sub)),
                                  match(colnames(target.tfs), colnames(target.tfs_sub))]
  
  colnames(target.tfs) = geneName_mapping$gene[match(colnames(target.tfs), geneName_mapping$original)]
  rownames(target.tfs) = geneName_mapping$gene[match(rownames(target.tfs), geneName_mapping$original)]
  colnames(target.tfs_sub) = colnames(target.tfs)
  rownames(target.tfs_sub) = rownames(target.tfs)
  
  # compute the group specificity 
  ggs = unique(c(colnames(target.tfs), rownames(target.tfs)))
  specificity = c()
  for(n in 1:length(ggs))
  {
    jj1 = which(colnames(target.tfs) == ggs[n])
    jj2 = which(rownames(target.tfs) == ggs[n])
    motif.total = sum(c(target.tfs[, jj1], target.tfs[jj2, ]), na.rm = TRUE)
    motif.group = sum(c(target.tfs_sub[, jj1], target.tfs_sub[jj2, ]), na.rm = TRUE)
    specificity = c(specificity, motif.group/motif.total)

  }
  names(specificity) = ggs
    
  ### prune edges based on the sub target-tf matrix
  edgeList.full = as_edgelist(trn)
  sub_trn = trn
  for(n in 1:nrow(edgeList.full))
  {
    # n = 1
    kk1 = which(colnames(target.tfs_sub) == edgeList.full[n, 1])
    kk2 = which(rownames(target.tfs_sub) == edgeList.full[n, 2])
    
    if(target.tfs[kk2, kk1] > 0){
      ratio = target.tfs_sub[kk2, kk1] / target.tfs[kk2, kk1]
      if(ratio < 0.1) {
        # cat(n, ' : delete edge -- ', paste0(edgeList.full[n, 1], '|', edgeList.full[n, 2]), '\n')
        sub_trn = delete_edges(sub_trn, paste0(edgeList.full[n, 1], '|', edgeList.full[n, 2]))
      }
    }
  }
  
  sub_trn = delete.vertices(simplify(sub_trn), degree(sub_trn)<2)
  
  nb_clusters = length(unique(V(trn)$cluster))
  cat(nb_clusters, ' clusters used here \n')
  
  if(length(nb_clusters)<=9) {
    library(khroma)
    muted <- colour("muted")
    got_palette = c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99")[1:nb_clusters] 
  }else{
    cat('More than 9 colors needed \n')
    
    got_palette = c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255", "#44AA99", 
                    "#1A5878", "#C44237", "#AD8941", "#E99093",
                    "#50594B", "#8968CD", "#9ACD32")[1:nb_clusters]
  }
  
  pal<-brewer.pal(nb_clusters, "Set3") # Vertex color assigned per each class number
  
  # manuallly control the module color
  modules_all = as.factor(unique(V(trn)$module))
  modules_sub = as.factor(unique(V(sub_trn)$module))
  
  got_palette = c("#CC6677", "#332288", "#DDCC77", "#117733", 
                  "#88CCEE", "#882255", "#44AA99", "#999999")[match(levels(modules_sub), levels(modules_all))]
  
  V(sub_trn)$degree<-degree(sub_trn, mode = 'All')
  V(sub_trn)$degreeIn = degree(sub_trn, mode = 'in')
  V(sub_trn)$degreeOut = degree(sub_trn, mode = 'out')
  V(sub_trn)$name[which.max(V(sub_trn)$degree)]
  V(sub_trn)$name[which.max(V(sub_trn)$degreeIn)]
  V(sub_trn)$name[which.max(V(sub_trn)$degreeOut)]
  
  ggraph(sub_trn, x=umap1, y=umap2) +
    geom_edge_link(aes(edge_width = weight), edge_colour = "gray80" ) +
    #geom_edge_diagonal(color='darkgray', width=0.2) +
    geom_node_point(aes(fill = module, size = degree), shape = 21) +
    scale_size_continuous(range = c(1, 7)) +
    #geom_node_text(aes(filter = size >= 20, label = name), family = "serif") +
    geom_node_text(aes(filter = size >= 10, label=name), size=6/ggplot2::.pt, repel=T, family = "serif")+
    scale_edge_width_continuous(range = c(0.05, 0.2)) +
    scale_edge_color_gradientn(colors=rev(brewer.pal(n=11, name="RdBu")), limits=c(0.01, 0.6)) +
    #scale_edge_alpha_continuous(range=c(0.05, 0.4), limits=c(2,20)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) + 
    scale_fill_manual(values = got_palette) +
    #scale_fill_viridis(option='magma') + 
    #scale_fill_manual(values = pal) +
    coord_fixed() +
    theme_graph() +
    #theme(legend.position = "bottom") 
    theme(legend.position = "none") 
  
  ggsave(paste0(figureDir, "GRN_tfExpr.umap_connectivit.module_subgraph_", subgroup, "_withGeneLabels.pdf"),
         width=8, height = 6)
  
  library(viridis)
  ntop = 20
  xx = data.frame(degree = V(sub_trn)$degree, gene = V(sub_trn)$name)
  xx$specificity = as.numeric(specificity[match(xx$gene, names(specificity))])
  xx = xx[order(-xx$degree), ]
  xx = xx[c(1:ntop), ]
  
  xx$gene = sapply(xx$gene, function(x) unlist(strsplit(as.character(x), '_'))[1])
  
  as_tibble(xx) %>% 
    ggplot(aes(y=degree, x=reorder(gene, -degree), fill = specificity)) + 
    geom_bar(position="dodge", stat="identity") +
    theme_classic() +
    #theme(axis.text.x = element_text(angle = 90, size = 10)) +
    scale_fill_viridis_c(option = 'magma', direction = -1) +
    labs(x = '', y = 'number of connections') +
    theme(axis.text.x = element_text(angle = 60, size = 12, hjust = 1), 
          axis.text.y = element_text(angle = 0, size = 12), 
          axis.title =  element_text(size = 14),
          legend.text = element_text(size=12),
          legend.title = element_text(size = 14),
          legend.position = 'top',
          #plot.margin = margin()
          #legend.key.size = unit(1, 'cm')
          #legend.key.width= unit(1, 'cm')
    )
  
  ggsave(paste0(figureDir, "GRN_subgraph_centrality_all_", subgroup, ".pdf"),  width = 6, height = 4)
  
}


