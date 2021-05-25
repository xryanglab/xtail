#' Volcano plot
#'
#' A simple function that plots log2 fold change of TE and -log10 of pvalues
#' The genes are color coded by P-value (-log10)
#'
#' @param object a \code{xtail} object
#'
#' @param ... optional arguments. Currently not used
#'
#' @return a \code{ggplot} object
#'
#' @name volcanoPlot
#'
#' @examples
#' data(xtailres)
#' volcanoPlot(xtailres)
NULL

#' @rdname volcanoPlot
#' @export
setGeneric("volcanoPlot",
           function(object, ...)
             standardGeneric("volcanoPlot"))

#' @rdname volcanoPlot
#'
#' @import ggplot2
#'
#' @export
setMethod("volcanoPlot", signature = "xtail",
    function(object){
        ggplot(data = as.data.frame(resultsTable(object)),
               aes_string(x = "log2FC_TE_final",
                          y = "-log10(pvalue_final)"))+
            geom_point(colour="#836FFF", alpha = 0.65, shape = 19L) +
            theme_bw()
    }
)
