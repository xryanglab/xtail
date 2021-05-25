#' Scatterplot of log2 RPF-to-mRNA ratios.
#'
#' A simple function that plots the log2 RPF-to-mRNA ratios in two conditions.
#' Different classes of genes are labeled with different of color.
#' Each dot, representing a particular gene, are color coded by the category of
#' expression change.
#'
#' @param object a \code{xtail} object
#'
#' @param log2R.cutoff cutoff for categorizing results
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
#' @return a \code{ggplot} object
#'
#' @name plotRs
#'
#' @examples
#' data(xtailres)
#' plotRs(xtailres)
NULL

#' @rdname plotRs
#' @export
setGeneric("plotRs",
           function(object,...)
             standardGeneric("plotRs"))

#' @rdname plotRs
#'
#' @import ggplot2
#' @importFrom stats complete.cases
#'
#' @export
setMethod("plotRs", signature = "xtail",
    function(object,
             log2R.cutoff = 1,
             xlim = NULL,
             ylim = NULL){
        # input check
        if(!is.numeric(log2R.cutoff) ||
           length(log2R.cutoff) != 1L ||
           log2R.cutoff <= 0 ) {
            stop("log2R.cutoff must be a single numeric value and be larger ",
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
        table <- resultsTable(object, sort.by = "pvalue_final", log2Rs = TRUE)
        table <- table[complete.cases(as.data.frame(table)),]
        cond <- paste0(conditions(object),"_log2TE")
        categories <- c("condition1" = "#0D5FFF",
                        "condition2" = "#FF2E06",
                        "homodirectional" = "#75E805",
                        "opposite change" = "#FFDE13",
                        "stable" = "gray90")
        names(categories)[seq_len(2L)] <- paste0(conditions(object), " only")
        x <- cond[1L]
        y <- cond[2L]
        .plot_scatter(table, x, y, categories, log2R.cutoff, xlim, ylim)
    }
)
