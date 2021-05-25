#' xtail analysis wrapper for \code{SummarizedExperiment} objects
#'
#' The \code{\link[=xtail]{xtail}} function can be used directly with
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' objects. The \code{mrna} and \code{rpf} data must be stored as two separate
#' assays.
#'
#' See \code{\link[=xtail]{xtail}} for more details on the analysis
#' function.
#'
#' @param x a \code{RpfSummarizedExperiment} object
#'
#' @param mrna_assay a scalar character. The name of the \code{assay} containing
#'   the \code{mRNA} data.
#'
#' @param rpf_assay a scalar character. The name of the \code{assay} containing
#'   the \code{rpf} data.
#'
#' @param ...
#' \itemize{
#'   \item{For \code{runXtail}: additional parameters passed on to
#'     \code{\link[=xtail]{xtail}}.}
#'   \item{For \code{addXtail}: additional parameters passed on to
#'     \code{runXtail} and \code{\link[=xtail]{xtail}}.}
#' }
#'
#' @return A \code{\link[=xtail]{xtail}} object for \code{runXtail} or
#'   an object of \code{class(x)}.
#'
#' @name runXtail
#'
#' @examples
#' data(xtaildata)
#' test.mrna <- xtaildata$mrna[1:100,]
#' test.rpf <- xtaildata$rpf[1:100,]
#' condition <- c("control","control","treat","treat")
#'
#' se <- SummarizedExoeriment(assays = list(mrna = test.mrna, rpf = test.rpf),
#'                            colData = DataFrame(condition = condition))
#' xtail <- runXtail(se, "mrna", "rpf", bins = 1000, threads = 2)
#' xtail
#'
#' se <- addXtail(se, "mrna", "rpf", bins = 1000, threads = 2)
#' rowData(se)
NULL

#' @rdname runXtail
#' @export
setGeneric("runXtail", signature = c("x"),
    function(x, mrna_assay, rpf_assay, ...) standardGeneric("runXtail")
)

#' @rdname runXtail
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("runXtail", signature = c(x = "SummarizedExperiment"),
    function(x, mrna_assay, rpf_assay, ...){
        mrna <- assay(x,mrna_assay)
        rpf <- assay(x,rpf_assay)
        mrna <- .remove_zero_and_NA_rows(mrna)
        rpf <- .remove_zero_and_NA_rows(rpf)
        xtail(mrna = mrna, rpf = rpf, ...)
    }
)

.remove_zero_and_NA_rows <- function(mat){
  mat <- mat[!apply(apply(mat,1,is.na),2,all),]
  mat <- mat[!apply(apply(mat,1,"==",0),2,all),]
  mat
}

#' @rdname runXtail
#' @export
setGeneric("addXtail", signature = c("x"),
           function(x, ...) standardGeneric("addXtail")
)

#' @rdname runXtail
#' @export
setMethod("addXtail", signature = c(x = "SummarizedExperiment"),
    function(x,  ...){
        xtail <- runXtail(x, ...)
        .add_xtail_results(x, xtail)
    }
)

#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom S4Vectors metadata metadata<-
.add_xtail_results <- function(x, xtail){
    # add xtail results
    rd <- rowData(x)
    table <- resultsTable(xtail)
    rd[,colnames(table)] <- NA
    rd[rownames(table),colnames(table)] <- table
    rowData(x) <- rd
    # add xtail metadata
    metadata(x) <- c(metadata(x),
                     list(num = resultsNum(xtail),
                          conditions = conditions(xtail)))
    #
    x
}
