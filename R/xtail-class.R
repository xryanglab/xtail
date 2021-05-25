#' @rdname xtail
#' @export
setClass("xtail",
         contains = "SimpleList")

################################################################################
# constructor

# only to be used internally

#' @importFrom S4Vectors DataFrame
.combine_results <- function(result_log2FC, result_log2R, condition, baseLevel){
    if(is.null(rownames(result_log2FC)) || is.null(rownames(result_log2R))){
        stop("Something went wrong")
    }
    # get the intersection of results and get them in the same order
    intersect.genes <- intersect(rownames(result_log2FC),
                                 rownames(result_log2R))
    result_log2R <- result_log2R[intersect.genes,]
    result_log2FC <- result_log2FC[intersect.genes,]

    #result data frame
    condition1_TE <- paste0(baseLevel,"_log2TE")
    condition2_TE <- paste0(unique(condition)[2], "_log2TE")
    res <- cbind(result_log2FC[,c("log2Ratio1","log2Ratio2","deltaTE","Pval")],
                 result_log2R[,c("log2Ratio1","log2Ratio2","deltaTE","Pval")])
    colnames(res) <- c("mRNA_log2FC","RPF_log2FC","log2FC_TE_v1","pvalue_v1",
                       condition1_TE,condition2_TE,"log2FC_TE_v2","pvalue_v2")
    res <- DataFrame(res)
    res
}

#' @importFrom stats p.adjust
.adjust_p <- function(res, method){
    res$pvalue.adjust <- p.adjust(res$pvalue_final, method = method)
    res
}

#' @importFrom S4Vectors SimpleList
.new_xtail <- function(result_log2FC, result_log2R, condition, baseLevel, ci,
                       p.adjust.method){
    # input check
    if (length(ci) != 1L || ci<0 || ci > 1 ) {
        stop("ci must be a single numeric value between 0 and 1.",
             call. = FALSE)
    }
    if (length(unique(condition))!=2) {
        stop("There must be exactly two different conditions.",
             call. = FALSE)
    }
    if (length(baseLevel) != 1L || !(baseLevel %in% condition)){
        stop("baseLevel is not in condition", call. = FALSE)
    }
    #
    res <- .combine_results(result_log2FC, result_log2R, condition, baseLevel)
    res$log2FC_TE_final <- 0
    res$pvalue_final <- 0
    res$pvalue.adjust <- 0
    log2FC_determine_num <- 0
    log2R_determine_num <- 0
    # empty results
    na_result <- is.na(res$pvalue_v1) | is.na(res$pvalue_v2)
    res[na_result,"log2FC_TE_final"] <- NA
    res[na_result,"pvalue_final"] <- NA
    # fold change better than ratio
    fc_result <- res$pvalue_v1 > res$pvalue_v2
    log2FC_determine_num <- sum(fc_result)
    res[fc_result,"log2FC_TE_final"] <- res[fc_result,"log2FC_TE_v1"]
    res[fc_result,"pvalue_final"] <- res[fc_result,"pvalue_v1"]
    # ratio better than fold change
    ratio_result <- res$pvalue_v1 <= res$pvalue_v2
    log2R_determine_num <- sum(ratio_result)
    res[ratio_result,"log2FC_TE_final"] <- res[ratio_result,"log2FC_TE_v2"]
    res[ratio_result,"pvalue_final"] <- res[ratio_result,"pvalue_v2"]
    #
    if(ci > 0){
        CI_string = paste0("CI(",100*ci,"%)")
        res[[CI_string]] <- NA
        res[fc_result,CI_string] <- paste0("[",
                                           round(result_log2FC[fc_result,"lowci"],2L),
                                           ",",
                                           round(result_log2FC[fc_result,"highci"],2L),
                                           "]")
        res[ratio_result,CI_string] <- paste0("[",
                                           round(result_log2R[ratio_result,"lowci"],2L),
                                           ",",
                                           round(result_log2R[ratio_result,"highci"],2L),
                                           "]")
    }
    res <- .adjust_p(res, method = p.adjust.method)
    # construct result object
    xtail <- SimpleList(resultsTable = res,
                        log2FC_determine_num = log2FC_determine_num,
                        log2R_determine_num = log2R_determine_num,
                        condition1 = baseLevel,
                        condition2 = unique(condition[condition != baseLevel]) )
    class(xtail) <- "xtail"
    validObject(xtail)
    xtail
}


################################################################################
# validity

XTAIL_OBJ_NAMES <- c("resultsTable","log2FC_determine_num",
                     "log2R_determine_num","condition1","condition2")

.valid_xtail <- function(object){
    f <- XTAIL_OBJ_NAMES %in% names(object)
    if(!all(f)){
        return(paste0("Missing elements in xtail object: '",
                      paste(XTAIL_OBJ_NAMES[f],"', '"),
                      "'"))
    }
    if(!is(object$resultsTable,"DataFrame")){
        return("results table is not a DataFrame")
    }
    if(!is.integer(object$log2FC_determine_num) ||
       length(object$log2FC_determine_num) != 1L){
        return("log2FC_determine_num is not an integer value")
    }
    if(!is.integer(object$log2R_determine_num)||
       length(object$log2R_determine_num) != 1L){
        return("log2R_determine_num is not an integer value")
    }
    if(!is.character(object$condition1) ||
       length(object$condition1) != 1L){
        return("condition1 is not a character value")
    }
    if(!is.character(object$condition2)||
       length(object$condition2) != 1L){
        return("condition2 is not a character value")
    }
    return(NULL)
}

setValidity("xtail", .valid_xtail)

################################################################################
# accessors

#' Results table of xtail results
#'
#' To retrieve the results from the xtail run use one of the accessor functions.
#'
#'
#' @param object a \code{xtail} object
#'
#' @param sort.by the column to sort with. Defaults to \code{NULL} to disable
#'   sorting.
#'
#' @param log2FCs \code{TRUE} or \code{FALSE}: Should log2 fold change values be
#'   returned? (defaults to \code{TRUE})
#'
#' @param log2Rs \code{TRUE} or \code{FALSE}: Should log2 ratio values be
#'   returned? (defaults to \code{TRUE})
#'
#' @param ... optional arguments. Currently not used
#'
#' @return a \code{DataFrame} with the results or numeric vectors
#'
#' @name xtail-accessors
#'
#' @aliases resultsTable
#'
#' @examples
#' data(xtailres)
#' resultsTable(xtailres)
#' conditions(xtailres)
#' resultsNum(xtailres)
#'
#' # sorting or results
#' resultsTable(xtailres, sort.by = "pvalue.adjust")
NULL

#' @rdname xtail-accessors
#' @importFrom BiocGenerics conditions
#' @export
setMethod("conditions", signature = c(object="xtail"),
    function(object){
        c(object$condition1,object$condition2)
    }
)

#' @rdname xtail-accessors
#' @export
setGeneric("resultsNum", signature = c("object"),
           function(object, ...)
               standardGeneric("resultsNum"))

#' @rdname xtail-accessors
#' @export
setMethod("resultsNum", signature = c("xtail"),
    function(object, ...){
        c(numFoldChange = object$log2FC_determine_num,
          numRatio = object$log2R_determine_num)
    }
)

#' @rdname xtail-accessors
#' @export
setGeneric("resultsTable",
           function(object,...)
             standardGeneric("resultsTable"))

#' @rdname xtail-accessors
#' @export
setMethod("resultsTable", signature = c(object="xtail"),
    function(object,
             sort.by = NULL,
             log2FCs = FALSE,
             log2Rs = FALSE){
        # input check
        if(!is.logical(log2FCs) || length(log2FCs) != 1L || is.na(log2FCs)){
            stop("'log2FCs' must be TRUE or FALSE",call. = FALSE)
        }
        if(!is.logical(log2Rs) || length(log2Rs) != 1L || is.na(log2Rs)){
            stop("'log2Rs' must be TRUE or FALSE",call. = FALSE)
        }
        #
        x <- object$resultsTable
        if(!is.null(sort.by)){
            sort.by <- match.arg(sort.by, c("pvalue.adjust", "log2FC_TE_v1",
                                            "log2FC_TE_v2", "log2FC_TE_final",
                                            "pvalue_v1" , "pvalue_v2",
                                            "pvalue_final"))
            # sort
            #absolute TE fold change.
            if(sort.by %in% c("log2FC_TE_v1","log2FC_TE_v2","log2FC_TE_final")){
                tefc <- abs(x[[sort.by]])
            } else {
                tefc <- abs(x[["log2FC_TE_final"]])
            }

            #pvalue order
            if(sort.by %in% c("pvalue_v1","pvalue_v2","pvalue_final", "pvalue.adjust")) {
                pvfc <- object$resultsTable[[sort.by]]
            }

            # out
            o <- switch(sort.by,
                        "log2FC_TE_v1" = order(tefc, decreasing = TRUE),
                        "log2FC_TE_v2" = order(tefc, decreasing = TRUE),
                        "log2FC_TE_final" = order(tefc, decreasing = TRUE),
                        "pvalue_v1" = order(pvfc, -tefc),
                        "pvalue_v2" = order(pvfc, -tefc),
                        "pvalue_final" = order(pvfc, -tefc),
                        "pvalue.adjust" = order(pvfc, -tefc)
            )
            x <- x[o,]
        }
        # remove deselected columns
        if (!log2Rs){
            cond <- paste0(conditions(object), "_log2TE")
            x[[cond[1L]]] <- NULL
            x[[cond[1L]]] <- NULL
        }
        if (!log2FCs){
            x$mRNA_log2FC <- NULL
            x$RPF_log2FC <- NULL
        }
        #
        x
    }
)

################################################################################
# show & summary

.show_xtail <- function(object, alpha = 0.1){
    table <- resultsTable(object)
    nums <- resultsNum(object)
    all_num <- nrow(table)
    up_num <- sum(table$log2FC_TE_final > 0 & table$pvalue.adjust < alpha,
                  na.rm = TRUE)
    down_num <- sum(table$log2FC_TE_final < 0 & table$pvalue.adjust < alpha,
                    na.rm = TRUE)
    cat("A xtail object:\n")
    cat("Number of genes tested:", all_num,"\n")
    cat("Number of the log2FC and log2R used in determining the final p-value:\n")
    cat("     log2FC:", nums["numFoldChange"],"\n")
    cat("     log2R :", nums["numRatio"],"\n")
    cat("\n")
    cat("Number of result with adjusted pvalue < ", alpha, "\n")
    cat("     log2FC_TE > 0 (up)  :", up_num, "\n")
    cat("     log2FC_TE < 0 (down):", down_num,"\n")
}

setMethod("show", signature = "xtail",
    function(object){
        .show_xtail(object)
    }
)

#' @rdname xtail-accessors
#'
#' @param alpha cut off for summarizing results. Only results with a adjusted
#'   p-value lower than \code{alpha} will be reported.
#'
#' @export
summary.xtail <- function(object, alpha = 0.1, ...){
    .show_xtail(object, alpha = alpha)
}

################################################################################

#' Write xtail results as table
#'
#' xtail results can be directly written to file using \code{write.xtail}
#'
#' @param object a \code{xtail} object
#'
#' @param ... arguments passed onto \code{\link[=resultsTable]{resultsTable}}
#'
#' @param write.args a list of arguments passed onto \code{write.table}
#'
#' @return invisible result from \code{write.table}
#'
#' @name write.xtail
#'
#' @seealso
#' \code{\link[utils:write.table]{write.table}}
#'
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' data(xtestres)
#' write.xtail(xtestres, file = tempfile())
write.xtail <- function(object, ..., write.args = list()){
    # input check
    if(!is.list(write.args)){
        stop("'write.args' must be a list.", call. = FALSE)
    }
    #
    x <- resultsTable(object, ...)
    do.call(write.table, c(list(x = x),write.args))
}
