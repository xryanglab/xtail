#' A tool to quantitatively assess differential translations with ribosome profiling data.
#'
#' By pairwise comparisons of ribosome profiling data, Xtail
#' identifies differentially translated genes across two experimental or
#' physiological conditions.
#'
#'
#' @references Zhengtao Xiao, Qin Zou, Yu Liu, and Xuerui Yang: Genome-wide
#'   assessment of differential translations with ribosome profiling data.
#'
#' @docType package
#' @name xtail-package
#' @author Zhengtao xiao
NULL

#' @import methods
#' @import Rcpp
#' @useDynLib xtail
NULL

#' xtail example data
#'
#' The \code{xtail} package includes some example data
#'
#' @name xtaildata
#' @format a list of matrices
#' @usage data(xtaildata)
"xtaildata"
#' @rdname xtaildata
#' @format an example results
#' @usage data(xtailres)
"xtailres"
