
#' Scatterplot of fold changes.
#'
#' A simple function that plots the fold change in mRNA and RPF.
#' Different classes of genes are labeled with different of color.
#' Each dot, representing a particular gene, are color coded by P-value (-log10)
#' of its differential translation across two conditions, which was estimated with xtail.
#'
#' @docType methods
#' @name plotFCs
#' @rdname plotFCs
#' @aliases plotFCs
#'
#' @importFrom scales
#'
#' @param object a xtailResults
#' @ explort

plotFCs <- function(object, log2FC.cutoff = 1, cex=1, xlim, ylim, ..., cex.lab, cex.axis, cex.main, cex.sub, pch)
{
  list.of.packages <- "scales"
  if(list.of.packages %in% installed.packages()[,"Package"]){
    library("scales")
  }else{
    message('please install "scales" before run this function')
  }

  if(is.null(object$resultsTable)) stop("error, the object must be created by using the xtail function.")
  if(!is(object, "xtailResults")) stop("error, the object must be created by using the xtail function.")

	if (log2FC.cutoff <= 0 ) stop("log2FC.cutoff must be larger than 0")
	if (!missing(xlim)){
		if (length(xlim) != 2) stop("invalid xlim value, should have two elemtnts")
	}
	if (!missing(ylim)){
		if (length(ylim) != 2) stop("invalid ylim value, should have two elemtnts")
	}
	resultsTable <- object$resultsTable
	resultsTable <- resultsTable[complete.cases(resultsTable),]
	resultsTable <- resultsTable[order(resultsTable$pvalue_final, decreasing = TRUE),]
	resultsTable$colorscale <- rescale(-log10(resultsTable$pvalue_final))

	if (missing(xlim)){
		xmin <- min(resultsTable$mRNA_log2FC,na.rm=T) - 0.1
		xmax <- max(resultsTable$mRNA_log2FC,na.rm=T) + 0.1
		xlim <- c(xmin,xmax)
	}
	if (missing(ylim)){
		ymin <- min(resultsTable$RPF_log2FC,na.rm=T) - 0.1
		ymax <- max(resultsTable$RPF_log2FC,na.rm=T) + 0.1
		ylim <- c(ymin,ymax)
	}
	if (missing(cex.lab)) cex.lab <- 1.2
	if (missing(cex.axis)) cex.axis <- 1
	if (missing(cex.main)) cex.main <- 1.2
	if (missing(cex.sub)) cex.sub <- 1
	if (missing(pch)) pch <- 20

	#mRNA stable, RPF stable.
	variable <- which(abs(resultsTable$mRNA_log2FC - resultsTable$RPF_log2FC) <= log2FC.cutoff |
				(abs(resultsTable$mRNA_log2FC) < log2FC.cutoff & abs(resultsTable$RPF_log2FC) < log2FC.cutoff))
	plot(resultsTable$mRNA_log2FC[variable],resultsTable$RPF_log2FC[variable],pch=pch,col="gray90",
		xlim=xlim,ylim=ylim,xlab = "mRNA log2FC", ylab="RPF log2FC",frame.plot=TRUE,cex=cex,
		cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, cex.sub=cex.sub, ...)

	#mRNA change, RPF stable.
	variable <- which( abs(resultsTable$mRNA_log2FC) >= log2FC.cutoff &
	 			abs(resultsTable$RPF_log2FC) < log2FC.cutoff &
	 			abs(resultsTable$mRNA_log2FC - resultsTable$RPF_log2FC) >= log2FC.cutoff)
	colornames <- seq_gradient_pal("#BFBBFF","#0D5FFF")(resultsTable$colorscale[variable])
	points(resultsTable$mRNA_log2FC[variable],resultsTable$RPF_log2FC[variable],
		pch=pch,col=alpha(colornames,0.8),cex=cex)
	leg <- "transcription only"
	leg.col <- "#0D5FFF"

	#mRNA stable, RPF change.
	variable <- which( abs(resultsTable$mRNA_log2FC) < log2FC.cutoff &
	 			abs(resultsTable$RPF_log2FC) >= log2FC.cutoff &
	 			abs(resultsTable$mRNA_log2FC - resultsTable$RPF_log2FC) >= log2FC.cutoff )
	colornames <- seq_gradient_pal("#FFB69C","#FF2E06")(resultsTable$colorscale[variable])
	points(resultsTable$mRNA_log2FC[variable],resultsTable$RPF_log2FC[variable],pch=pch,col=alpha(colornames,0.8),cex=cex)
	leg <- c(leg, "translation only")
	leg.col <- c(leg.col, "#FF2E06")

	#mRNA change, RPF change, homodirectional.
	variable <- which( sign(resultsTable$mRNA_log2FC) * sign(resultsTable$RPF_log2FC) > 0 &
				abs(resultsTable$mRNA_log2FC) >= log2FC.cutoff &
				abs(resultsTable$RPF_log2FC) >= log2FC.cutoff &
	 			abs(resultsTable$mRNA_log2FC - resultsTable$RPF_log2FC) >= log2FC.cutoff )
	colornames <- seq_gradient_pal("#D1F7AD","#75E805")(resultsTable$colorscale[variable])
	points(resultsTable$mRNA_log2FC[variable],resultsTable$RPF_log2FC[variable],pch=pch,col=alpha(colornames,0.8),cex=cex)
	leg <- c(leg, "homodirectional")
	leg.col <- c(leg.col, "#75E805")

	#mRNA change, RPF change, opposite change.
	variable <- which( sign(resultsTable$mRNA_log2FC) * sign(resultsTable$RPF_log2FC) < 0 &
				abs(resultsTable$mRNA_log2FC) >= log2FC.cutoff &
				abs(resultsTable$RPF_log2FC) >= log2FC.cutoff &
	 			abs(resultsTable$mRNA_log2FC - resultsTable$RPF_log2FC) >= log2FC.cutoff )
	colornames <- seq_gradient_pal("#FFF1AF","#FFDE13")(resultsTable$colorscale[variable])
	points(resultsTable$mRNA_log2FC[variable],resultsTable$RPF_log2FC[variable],pch=pch,col=alpha(colornames,0.8),cex=cex)
	leg <- c(leg, "opposite change")
	leg.col <- c(leg.col, "#FFDE13")

	abline(h= log2FC.cutoff, lty=2, col="gray")
	abline(v= log2FC.cutoff, lty=2, col="gray")
	abline(h= -log2FC.cutoff, lty=2, col="gray")
	abline(v= -log2FC.cutoff, lty=2, col="gray")

	legend("bottomright", legend=leg,pch=pch, col=leg.col, bty="n",cex=cex)
}


#' Scatterplot of log2 RPF-to-mRNA ratios.
#'
#' A simple function that plots the log2 RPF-to-mRNA ratios in two conditions.
#' Different classes of genes are labeled with different of color.
#' Each dot, representing a particular gene, are color coded by P-value (-log10)
#'
#' @docType methods
#' @name plotRs
#' @rdname plotRs
#' @aliases plotRs
#'
#' @importFrom scales
#'
#' @param object a xtailResults
#' @ explort

plotRs <- function(object, log2R.cutoff = 1, cex=1, xlim, ylim, ..., cex.lab, cex.axis, cex.main, cex.sub, pch)
{
  list.of.packages <- "scales"
  if(list.of.packages %in% installed.packages()[,"Package"]){
    library("scales")
  }else{
    message('please install "scales" before run this function')
  }

  if(is.null(object$resultsTable)) stop("error, the object must be created by using the xtail function.")
  if(!is(object, "xtailResults")) stop("error, the object must be created by using the xtail function.")

  if (log2R.cutoff <= 0 ) stop("log2R.cutoff must be larger than 0")
  if (!missing(xlim)){
    if (length(xlim) != 2) stop("invalid xlim value, should have two elemtnts")
  }
  if (!missing(ylim)){
    if (length(ylim) != 2) stop("invalid ylim value, should have two elemtnts")
  }
  resultsTable <- object$resultsTable
  resultsTable <- resultsTable[complete.cases(resultsTable),]
  resultsTable <- resultsTable[order(resultsTable$pvalue_final, decreasing = TRUE),]
  resultsTable$colorscale <- rescale(-log10(resultsTable$pvalue_final))

  condition1_TE <- paste0(object$condition1,"_log2TE")
  condition2_TE <- paste0(object$condition2, "_log2TE")

  if (missing(xlim)){
    xmin <- min(resultsTable[[condition1_TE]],na.rm=T) - 0.1
    xmax <- max(resultsTable[[condition1_TE]],na.rm=T) + 0.1
    xlim <- c(xmin,xmax)
  }
  if (missing(ylim)){
    ymin <- min(resultsTable[[condition2_TE]],na.rm=T) - 0.1
    ymax <- max(resultsTable[[condition2_TE]],na.rm=T) + 0.1
    ylim <- c(ymin,ymax)
  }
  if (missing(cex.lab)) cex.lab <- 1.2
  if (missing(cex.axis)) cex.axis <- 1
  if (missing(cex.main)) cex.main <- 1.2
  if (missing(cex.sub)) cex.sub <- 1
  if (missing(pch)) pch <- 20

  #mRNA stable, RPF stable.
  variable <- which(abs(resultsTable[[condition1_TE]] - resultsTable[[condition2_TE]]) <= log2R.cutoff |
                      (abs(resultsTable[[condition1_TE]]) < log2R.cutoff & abs(resultsTable[[condition2_TE]]) < log2R.cutoff))
  plot(resultsTable[[condition1_TE]][variable],resultsTable[[condition2_TE]][variable],pch=pch,col="gray90",
       xlim=xlim,ylim=ylim,xlab = paste0(object$condition1," log2Rs"), ylab=paste0(object$condition2," log2Rs"),frame.plot=TRUE,cex=cex,
       cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, cex.sub=cex.sub, ...)

  #change, stable.
  variable <- which( abs(resultsTable[[condition1_TE]]) >= log2R.cutoff &
                       abs(resultsTable[[condition2_TE]]) < log2R.cutoff &
                       abs(resultsTable[[condition1_TE]] - resultsTable[[condition2_TE]]) >= log2R.cutoff)
  colornames <- seq_gradient_pal("#BFBBFF","#0D5FFF")(resultsTable$colorscale[variable])
  points(resultsTable[[condition1_TE]][variable],resultsTable[[condition2_TE]][variable],
         pch=pch,col=alpha(colornames,0.8),cex=cex)
  leg <- paste0(object$condition1," only")
  leg.col <- "#0D5FFF"

  # stable, change.
  variable <- which( abs(resultsTable[[condition1_TE]]) < log2R.cutoff &
                       abs(resultsTable[[condition2_TE]]) >= log2R.cutoff &
                       abs(resultsTable[[condition1_TE]] - resultsTable[[condition2_TE]]) >= log2R.cutoff )
  colornames <- seq_gradient_pal("#FFB69C","#FF2E06")(resultsTable$colorscale[variable])
  points(resultsTable[[condition1_TE]][variable],resultsTable[[condition2_TE]][variable],pch=pch,col=alpha(colornames,0.8),cex=cex)
  leg <- c(leg, paste0(object$condition2," only"))
  leg.col <- c(leg.col, "#FF2E06")

  # both change, homodirectional.
  variable <- which( sign(resultsTable[[condition1_TE]]) * sign(resultsTable[[condition2_TE]]) > 0 &
                       abs(resultsTable[[condition1_TE]]) >= log2R.cutoff &
                       abs(resultsTable[[condition2_TE]]) >= log2R.cutoff &
                       abs(resultsTable[[condition1_TE]] - resultsTable[[condition2_TE]]) >= log2R.cutoff )
  colornames <- seq_gradient_pal("#D1F7AD","#75E805")(resultsTable$colorscale[variable])
  points(resultsTable[[condition1_TE]][variable],resultsTable[[condition2_TE]][variable],pch=pch,col=alpha(colornames,0.8),cex=cex)
  leg <- c(leg, "homodirectional")
  leg.col <- c(leg.col, "#75E805")

  # opposite change.
  variable <- which( sign(resultsTable[[condition1_TE]]) * sign(resultsTable[[condition2_TE]]) < 0 &
                       abs(resultsTable[[condition1_TE]]) >= log2R.cutoff &
                       abs(resultsTable[[condition2_TE]]) >= log2R.cutoff &
                       abs(resultsTable[[condition1_TE]] - resultsTable[[condition2_TE]]) >= log2R.cutoff )
  colornames <- seq_gradient_pal("#FFF1AF","#FFDE13")(resultsTable$colorscale[variable])
  points(resultsTable[[condition1_TE]][variable],resultsTable[[condition2_TE]][variable],pch=pch,col=alpha(colornames,0.8),cex=cex)
  leg <- c(leg, "opposite change")
  leg.col <- c(leg.col, "#FFDE13")

  abline(h= log2R.cutoff, lty=2, col="gray")
  abline(v= log2R.cutoff, lty=2, col="gray")
  abline(h= -log2R.cutoff, lty=2, col="gray")
  abline(v= -log2R.cutoff, lty=2, col="gray")

  legend("bottomright", legend=leg,pch=pch, col=leg.col, bty="n",cex=cex)
}

#' volcano plot
#'
#' A simple function that plots log2 fold change of TE and -log10 of pvalues
#' The genes are color coded by P-value (-log10)
#'
#'
#' @docType methods
#' @name volcanoPlot
#' @rdname volcanoPlot
#' @aliases volcanoPlot
#'
#' @importFrom ggplot2
#'
#' @param object a xtailResults
#' @ explort

volcanoPlot <- function(object){
  list.of.packages <- "ggplot2"
  if(list.of.packages %in% installed.packages()[,"Package"]){
    library("ggplot2")
  }else{
    message('please install "ggplot2" before run this function')
  }

  p <- ggplot(data = object$resultsTable, aes(x=log2FC_TE_final, y=-log10(pvalue_final)))
  p <- p + geom_point(size=2.5, colour="#836FFF")
  p <- p+ theme_bw()
  p
}
