
#' @rdname dispersions
#' @export
setGeneric("dispersionMatrix", function(object,...) standardGeneric("dispersionMatrix"))

#' @rdname dispersions
#' @export
setGeneric("dispersionMatrix<-", function(object,...,value) standardGeneric("dispersionMatrix<-"))

#' @rdname resultsTable
#' @export
setGeneric("resultsTable", function(object,...) standardGeneric("resultsTable"))

#' @rdname plotFCs
#' @export
setGeneric("plotFCs", function(object,...) standardGeneric("plotFCs"))

#' @rdname plotFCs
#' @export
setGeneric("plotRs", function(object,...) standardGeneric("plotRs"))

#' @rdname vocanoPlot
#' @export
setGeneric("volcanoPlot", function(object,...) standardGeneric("volcanoPlot"))

#' @rdname xtailResults
#' @export
setClass("xtailResults")
