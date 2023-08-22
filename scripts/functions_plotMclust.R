##########################################################################
##########################################################################
# Project: RA competence 
# Script purpose: manually modify the mclust plotting function
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Nov 15 11:50:27 2022
##########################################################################
##########################################################################
plot.Mclust_cutomized <- function(x, 
                        what = c("BIC", "classification", "uncertainty", "density"), 
                        dimens = NULL, xlab = NULL, ylab = NULL,  
                        addEllipses = TRUE, main = FALSE, cex_clusterlabels = 1.2,
                        ...) 
{
  # x = mb;dimens = NULL; addEllipses = TRUE;
  object <- x # Argh.  Really want to use object anyway
  if(!inherits(object, "Mclust")) 
    stop("object not of class 'Mclust'")
  
  data <- object$data
  p <- ncol(data)
  
  if(p == 1) 
    colnames(data) <- deparse(object$call$data)
  dimens <- if(is.null(dimens)) seq(p) else dimens[dimens <= p]
  d <- length(dimens)
  
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
  
  what <- match.arg(what, several.ok = TRUE)
  oldpar <- par(no.readonly = TRUE)
  
  plot.Mclust.bic <- function(...)
    plot.mclustBIC(object$BIC, xlab = xlab, ...)
  
  plot.Mclust.classification <- function(...)
  {  
    if(d == 1)
    { 
      mclust1Dplot(data = data[,dimens,drop=FALSE], 
                   what = "classification",
                   classification = object$classification,
                   z = object$z, 
                   xlab = if(is.null(xlab)) colnames(data)[dimens] else xlab, 
                   main = main, ...) 
    }
    if(d == 2) 
    { 
      pars <- object$parameters
      pars$mean <- pars$mean[dimens,,drop=FALSE]
      pars$variance$d <- length(dimens)
      pars$variance$sigma <- pars$variance$sigma[dimens,dimens,,drop=FALSE]
      mclust2Dplot(data = data[,dimens,drop=FALSE], 
                   what = "classification", 
                   classification = object$classification, 
                   parameters = if(addEllipses) pars else NULL,
                   xlab = if(is.null(xlab)) colnames(data)[dimens][1] else xlab, 
                   ylab = if(is.null(ylab)) colnames(data)[dimens][2] else ylab,
                   main = main, ...) 
    }
    if(d > 2)
    { 
      pars <- object$parameters
      pars$mean <- pars$mean[dimens,,drop=FALSE]
      pars$variance$d <- length(dimens)
      pars$variance$sigma <- pars$variance$sigma[dimens,dimens,,drop=FALSE]
      on.exit(par(oldpar))
      
      par(mfrow = c(d, d), 
          mar = rep(0.2/2,4), 
          oma = rep(3,4))
      for(i in seq(d))
      {
        for(j in seq(d))
        {
          if(i == j)
          {
            plot(data[, dimens[c(j, i)]],
                 type = "n", xlab = "", ylab = "", axes = FALSE)
            text(mean(par("usr")[1:2]), mean(par("usr")[3:4]),
                 labels = colnames(data[, dimens])[i],
                 cex = 1.5, adj = 0.5)
            box()
          }else{
            # i = 1; j = 2
            centers = object$parameters$mean[c(j,i),]
            coordProj(data = data, 
                      dimens = dimens[c(j,i)], 
                      what = "classification", 
                      classification = object$classification,
                      parameters = object$parameters,
                      addEllipses = addEllipses,
                      main = FALSE, xaxt = "n", yaxt = "n", ...)
            text(x = centers[1, ], y = centers[2, ], labels = c(1:ncol(centers)), 
                 cex = cex_clusterlabels)
            
          }
          if(i == 1 && (!(j%%2))) axis(3)
          if(i == d && (j%%2))    axis(1)
          if(j == 1 && (!(i%%2))) axis(2)
          if(j == d && (i%%2))    axis(4)
        }
      }
    }
  }
  
  plot.Mclust.uncertainty <- function(...) 
  {
    pars <- object$parameters
    if(d > 1)
    {
      pars$mean <- pars$mean[dimens,,drop=FALSE]
      pars$variance$d <- length(dimens)
      pars$variance$sigma <- pars$variance$sigma[dimens,dimens,,drop=FALSE]
    }
    #
    if(p == 1 || d == 1)
    { 
      mclust1Dplot(data = data[,dimens,drop=FALSE], 
                   what = "uncertainty", 
                   parameters = pars, z = object$z, 
                   xlab = if(is.null(xlab)) colnames(data)[dimens] else xlab, 
                   main = main, ...) 
    }
    if(p == 2 || d == 2)
    { 
      mclust2Dplot(data = data[,dimens,drop=FALSE], 
                   what = "uncertainty", 
                   parameters = pars,
                   # uncertainty = object$uncertainty,
                   z = object$z,
                   classification = object$classification,
                   xlab = if(is.null(xlab)) colnames(data)[dimens][1] else xlab, 
                   ylab = if(is.null(ylab)) colnames(data)[dimens][2] else ylab,
                   addEllipses = addEllipses, main = main, ...)
    }
    if(p > 2 && d > 2)
    { 
      on.exit(par(oldpar))
      par(mfrow = c(d, d), 
          mar = rep(0,4),
          mar = rep(0.2/2,4), 
          oma = rep(3,4))
      for(i in seq(d))
      { 
        for(j in seq(d)) 
        { 
          if(i == j) 
          { 
            plot(data[, dimens[c(j, i)]], type="n",
                 xlab = "", ylab = "", axes = FALSE)
            text(mean(par("usr")[1:2]), mean(par("usr")[3:4]),
                 labels = colnames(data[,dimens])[i], 
                 cex = 1.5, adj = 0.5)
            box()
          } else 
          { 
            coordProj(data = data, 
                      what = "uncertainty", 
                      parameters = object$parameters,
                      # uncertainty = object$uncertainty,
                      z = object$z,
                      classification = object$classification,
                      dimens = dimens[c(j,i)], 
                      main = FALSE, 
                      addEllipses = addEllipses,
                      xaxt = "n", yaxt = "n", ...)
          }
          if(i == 1 && (!(j%%2))) axis(3)
          if(i == d && (j%%2))    axis(1)
          if(j == 1 && (!(i%%2))) axis(2)
          if(j == d && (i%%2))    axis(4)
        }
      }
    }
  }
  
  plot.Mclust.density <- function(...)
  {
    if(p == 1)
    { 
      objdens <- as.densityMclust(object)
      plotDensityMclust1(objdens, 
                         xlab = if(is.null(xlab)) colnames(data)[dimens] else xlab, 
                         main = if(main) main else NULL, ...) 
      # mclust1Dplot(data = data,
      #              parameters = object$parameters,
      #              # z = object$z, 
      #              what = "density", 
      #              xlab = if(is.null(xlab)) colnames(data)[dimens] else xlab, 
      #              main = main, ...) 
    }
    if(p == 2) 
    { surfacePlot(data = data, 
                  parameters = object$parameters,
                  what = "density", 
                  xlab = if(is.null(xlab)) colnames(data)[1] else xlab, 
                  ylab = if(is.null(ylab)) colnames(data)[2] else ylab,
                  main = main, ...) 
    }
    if(p > 2) 
    { 
      objdens <- as.densityMclust(object)
      objdens$data <- objdens$data[,dimens,drop=FALSE]
      objdens$varname <- colnames(data)[dimens]
      objdens$range <- apply(data, 2, range)
      objdens$d <- d
      objdens$parameters$mean <- objdens$parameters$mean[dimens,,drop=FALSE]
      objdens$parameters$variance$d <- d
      objdens$parameters$variance$sigma <- 
        objdens$parameters$variance$sigma[dimens,dimens,,drop=FALSE]
      # 
      if (d == 1)
        plotDensityMclust1(objdens, ...)
      else if (d == 2)
        plotDensityMclust2(objdens, ...)
      else
        plotDensityMclustd(objdens, ...)
    }
  }
  
  if(interactive() & length(what) > 1)
  { title <- "Model-based clustering plots:"
  # present menu waiting user choice
  choice <- menu(what, graphics = FALSE, title = title)
  while(choice != 0)
  { if(what[choice] == "BIC")            plot.Mclust.bic(...)
    if(what[choice] == "classification") plot.Mclust.classification(...)
    if(what[choice] == "uncertainty")    plot.Mclust.uncertainty(...)
    if(what[choice] == "density")        plot.Mclust.density(...)
    # re-present menu waiting user choice
    choice <- menu(what, graphics = FALSE, title = title)
  }
  } 
  else 
  { if(any(what == "BIC"))            plot.Mclust.bic(...)
    if(any(what == "classification")) plot.Mclust.classification(...) 
    if(any(what == "uncertainty"))    plot.Mclust.uncertainty(...) 
    if(any(what == "density"))        plot.Mclust.density(...) 
  }
  
  invisible()
}


plot.surface_customized = function(x)
{s
  cat('to update')
}

