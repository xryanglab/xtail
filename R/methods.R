print.xtailResults <- function(){print ("xtailResults")}

dispersionMatrix.DESeqDataSet <- function(object){
  if (!"dispersionMatrix" %in% names(assays(object))) return (NULL)
  disp <- assays(object)[["dispersionMatrix"]]
  colnames(disp) <- colnames(object)
  disp
}

#' @export
setMethod("dispersionMatrix", signature(object="DESeqDataSet"),
          dispersionMatrix.DESeqDataSet)

#' @name dispersions
#' @rdname dispersions
#' @exportMethod "dispersions<-"
setReplaceMethod("dispersionMatrix", signature(object="DESeqDataSet", value="matrix"),
                 function(object, value) {
                  assays(object)[["dispersionMatrix"]] <- value
                  validObject( object )
                  object
                 })

#' Summarize xtail results
#' Print a summary of the results of xtail
#'
#' @usage
#' \method{summary}{xtailResults}(object,alpha)
#' @param object a \code{xtailResults} obect
#' @param alpha the adjusted p-value cutoff. If not set,
#' the default is 0.1.
#'
#' @docType methods
#' @name summary
#' @rdname summary
#' @aliases summary summary.xtailResults
#'
#' @export
summary.xtailResults <- function(object,alpha){
  if(is.null(object$resultsTable)) stop("error, the object must be created by using the xtail function.")
  if(!is(object, "xtailResults")) stop("error, the object must be created by using the xtail function.")
  if (missing(alpha)){
    alpha <- 0.1
  }
  resultsTable <- object$resultsTable
  cat("\n")
  all_num <- nrow(resultsTable)
  up_num <- sum(resultsTable$log2FC_TE_final > 0 & resultsTable$pvalue.adjust < alpha, na.rm=TRUE)
  down_num <- sum(resultsTable$log2FC_TE_final < 0 & resultsTable$pvalue.adjust < alpha, na.rm=TRUE)

  cat("The total number of gene is:", all_num,"\n")
  cat("Number of the log2FC and log2R used in determining the final p-value:\n")
  cat("     log2FC:", object$log2FC_determine_num,"\n")
  cat("     log2R :", object$log2R_determine_num,"\n")
  cat("\n")
  cat("adjusted pvalue < ", alpha, "\n")
  cat("     log2FC_TE > 0 (up)  :", up_num, "\n")
  cat("     log2FC_TE < 0 (down):", down_num,"\n")
  cat("\n")
}

resultsTable.xtailResults <- function(object,sort.by="pvalue.adjust", log2FCs = FALSE, log2Rs = FALSE){
  #check object
  if(!is(object, "xtailResults")) stop("error, the object must be created by using the xtail function.")
  if(is.null(object$resultsTable)) stop("error, the object must be created by using the xtail function.")

  if (missing(sort.by)){
    tab <-object$resultsTable
  }else{
    #check sort.by
    sort.by <- match.arg(sort.by, c("log2FC_TE_v1", "log2FC_TE_v2", "log2FC_TE_final" ,"pvalue_v1" , "pvalue_v2", "pvalue_final","pvalue.adjust"))

    #absolute TE fold change.
    if(is.element(sort.by, c("log2FC_TE_v1","log2FC_TE_v2","log2FC_TE_final"))){
      tefc <- abs(object$resultsTable[[sort.by]])
    }else{
      tefc <- abs(object$resultsTable[["log2FC_TE_final"]])
    }

    #pvalue order
    if(is.element(sort.by, c("pvalue_v1","pvalue_v2","pvalue_final", "pvalue.adjust"))) pvfc <- object$resultsTable[[sort.by]]

    # out
    o <- switch(sort.by,
                "log2FC_TE_v1" = order(tefc, decreasing=TRUE),
                "log2FC_TE_v2" = order(tefc, decreasing=TRUE),
                "log2FC_TE_final" = order(tefc, decreasing=TRUE),
                "pvalue_v1" = order(pvfc, -tefc),
                "pvalue_v2" = order(pvfc, -tefc),
                "pvalue_final" = order(pvfc, -tefc),
                "pvalue.adjust" = order(pvfc, -tefc)
                )
    tab <- object$resultsTable[o,]
  }
  #remove redundant columns
  if (!log2Rs){
    condition1 <- paste0(object$condition1, "_log2TE")
    condition2 <- paste0(object$condition2, "_log2TE")
    tab[[condition1]] <- NULL
    tab[[condition2]] <- NULL
  }
  if (!log2FCs){
    tab$mRNA_log2FC <- NULL
    tab$RPF_log2FC <- NULL
  }
  tab

}

#' Results table of xtail results
#' Print the results of xtail
#'
#' @usage
#' \method{resultsTable}{xtailResults}(object, sort.by)
#' @param object a \code{xtailResults} obect
#' @param sort.by the sorted column. By default, the pvalue.adjust is used.
#' @docType methods
#' @name resultsTable
#' @rdname resultsTable
#' @aliases resultsTable resultsTable.xtailResults
#'
#' @export
setMethod("resultsTable", signature(object="xtailResults"),resultsTable.xtailResults)


#' write xtail results as table
#'
#' @usage
#' \method{write.xtail}{xtailResults}(object)
#' @param object a \code{xtailResults} obect
#'
#' @docType methods
#' @name write.xtail
#' @rdname write.xtail
#' @aliases write.xtail write.table
#'
#' @export
#'
write.xtail <- function(object, ..., sort.by="pvalue.adjust", log2FCs = FALSE, log2Rs = FALSE){
  #check object
  if(!is(object, "xtailResults")) stop("error, the object must be created by using the xtail function.")
  if(is.null(object$resultsTable)) stop("error, the object must be created by using the xtail function.")

  if (missing(sort.by)){
    tab <- object$resultsTable
  }else{
    #check sort.by
    sort.by <- match.arg(sort.by, c("log2FC_TE_v1", "log2FC_TE_v2" ,"log2FC_TE_final", "pvalue_v1" , "pvalue_v2", "pvalue_final","pvalue.adjust"))

    #absolute TE fold change.
    if(is.element(sort.by, c("log2FC_TE_v1","log2FC_TE_v2","log2FC_TE_final"))){
      tefc <- abs(object$resultsTable[[sort.by]])
    }else{
      tefc <- abs(object$resultsTable[["log2FC_TE_final"]])
    }

    #pvalue order
    if(is.element(sort.by, c("pvalue_v1","pvalue_v2","pvalue_final", "pvalue.adjust"))) pvfc <- object$resultsTable[[sort.by]]

    # out
    o <- switch(sort.by,
                "log2FC_TE_v1" = order(tefc, decreasing=TRUE),
                "log2FC_TE_v2" = order(tefc, decreasing=TRUE),
                "log2FC_TE_final" = order(tefc, decreasing=TRUE),
                "pvalue_v1" = order(pvfc, -tefc),
                "pvalue_v2" = order(pvfc, -tefc),
                "pvalue_final" = order(pvfc, -tefc),
                "pvalue.adjust" = order(pvfc, -tefc)
    )
    tab <- object$resultsTable[o,]
  }

  #remove redundant columns
  if (!log2Rs){
    condition1 <- paste0(object$condition1, "_log2TE")
    condition2 <- paste0(object$condition2, "_log2TE")
    tab[[condition1]] <- NULL
    tab[[condition2]] <- NULL
  }
  if (!log2FCs){
    tab$mRNA_log2FC <- NULL
    tab$RPF_log2FC <- NULL
  }
  write.table(tab, ...)
}

