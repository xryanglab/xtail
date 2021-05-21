#' Scatterplot of fold changes.
#'
#' A simple function that plots the fold change in mRNA and RPF.
#' Different classes of genes are labeled with different of color.
#' Each dot, representing a particular gene, are color coded by the category of
#' expression change.
#'
#' @param object a \code{xtail} object
#'
#' @param log2FC.cutoff cutoff for categorizing results
#'
#' @param xlim argument for limiting the range plotted on the x axis. See
#'   \code{\link[ggplot2:scale_continuous]{scale_continuous}} for more
#'   details.
#'
#' @param ylim argument for limiting the range plotted on the y axis. See
#'   \code{\link[ggplot2:scale_continuous]{scale_continuous}} for more
#'   details.
#'
#' @param ... optional arguments. Currently not used
#'
#' @name plotFCs
NULL

#' @rdname plotFCs
#' @export
setGeneric("plotFCs",
           function(object,...)
             standardGeneric("plotFCs"))

#' @rdname plotFCs
#'
#' @import ggplot2
#' @importFrom stats complete.cases
#'
#' @export
setMethod("plotFCs", signature = "xtail",
    function(object,
             log2FC.cutoff = 1,
             xlim = NULL,
             ylim = NULL){
        # input check
        if(!is.numeric(log2FC.cutoff) ||
           length(log2FC.cutoff) != 1L ||
           log2FC.cutoff <= 0 ) {
            stop("log2FC.cutoff must be a single numeric value and be larger ",
                 "than 0.", call. = FALSE)
        }
        if(!is.null(xlim)){
            if(!is.numeric(xlim) || length(xlim) != 2L) {
                stop("invalid xlim value, should have two numeric elements.",
                     call. = FALSE)
            }
        }
        if(!is.null(ylim)){
            if(!is.numeric(ylim) || length(ylim) != 2L) {
                stop("invalid ylim value, should have two numeric elements.",
                     call. = FALSE)
            }
        }
        #
        table <- resultsTable(object, sort.by = "pvalue_final", log2FCs = TRUE)
        table <- table[complete.cases(as.data.frame(table)),]
        categories <- c("transcription only" = "#0D5FFF",
                        "translation only" = "#FF2E06",
                        "homodirectional" = "#75E805",
                        "opposite change" = "#FFDE13",
                        "stable" = "gray90")
        x <- "mRNA_log2FC"
        y <- "RPF_log2FC"
        .plot_scatter(table, x, y, categories, log2FC.cutoff, xlim, ylim)
    }
)
