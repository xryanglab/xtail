#' A tool to quantitatively assess differential translations with ribosome profiling data.
#'
#' By pairwise comparisons of ribosome profiling data, Xtail
#' identifies differentially translated genes across two experimental or
#' physiological conditions.
#'
#' @references Zhengtao Xiao, Qin Zou, Yu Liu, and Xuerui Yang: Genome-wide
#'   assessment of differential translations with ribosome profiling data.
#'
#' @param mrna a matrix or data frame of raw mRNA count data whose rows
#'   correspond to genes and columns correspond to samples. The column names
#'   should be non-empty, and in same order with condition.
#' @param rpf a matrix or data frame of raw RPF count data whose rows correspond
#'   to genes and columns correspond to samples.The column names should be
#'   non-empty, and in same order with condition.
#' @param condition condition labels corresponding to the order of samples in
#'   mrna and rpf. There must be exactly two unique values.
#' @param baseLevel The baseLevel indicates which one of the two conditions will
#'   be compared against by the other one. If not specified, \code{Xtail} will
#'   return results of comparing the second condition over the first one.
#' @param minMeanCount \code{Xtail} uses the average expression level of each
#'   gene, across all samples as filter criterion and it omits all genes with
#'   mean counts below minMeanCount.
#' @param ci The level of confindence to get credible intervals of log2 fold
#'   change of translational efficiency (TE), for example 0.95.
#' @param normalize Whether normalization should be done (TRUE \\ FALSE). If
#'   missing, \code{Xtail} will perform median-of-ratios normlazation by
#'   default.
#' @param p.adjust.method The method to use for adjusting multiple comparisons, by
#'   default "BH",  see \code{?p.adjust}
#' @param threads The number of CPU cores used. By default, all available  cores
#'   are used.
#' @param bins The number of bins used for calculating the probability density
#'   of log2FC or log2R (default is 10000). This paramater will determine
#'   accuracy of pvalue. Set it small for a very quick test run.
#'
#' @details No missing values are allowed in input data mrna and rpf.
#'
#' Duplicate row names (gene names or gene ids) are not allowed.
#'
#' \code{Xtail} takes in raw read counts of RPF and mRNA, and performs
#' median-of-ratios normalization. Alternatively, users can provide normalized
#' read counts and skip the built-in normal by setting "normalize" to FALSE.
#'
#' The step of estimation of the probability distributions, for log2FC or
#' log2R, will execute slowly in the current implementation, but can be
#' speeded up by running on multiple cores using the parallel library. By
#' default, the "detectCores" function in parallel library is used to
#' determine the number of CPU cores in the machine on which R is running. To
#' adjust the number of cores used, use "threads" argument to assign.
#'
#' @name xtail
#'
#' @author Zhengtao xiao
#'
#' @examples
#' #load the data
#' data(xtaildata)
#' # Get the mrna count data and rpf count data. For the example only the first
#' # 100 are used
#' test.mrna <- xtaildata$mrna[1:100,]
#' test.rpf <- xtaildata$rpf[1:100,]
#'
#' #Assign condition labels to samples.
#' condition <- c("control","control","treat","treat")
#'
#' #run xtail
#' test.results <- xtail(test.mrna,test.rpf,condition, threads = 2)
#' test.results
NULL

XTAIL_DISPERSION_ASSAY <- "dispersion"

.norm_mrna_rpf <- function(mrna, rpf, minMeanCount){
    #
    keep.genes <- rowMeans(mrna) >= minMeanCount
    mrna.keep <- rownames(mrna)[keep.genes]
    keep.genes <- rowMeans(rpf) >= minMeanCount
    rpf.keep <- rownames(rpf)[keep.genes]
    ## merge genes in mrna and rpf
    keep.genes <- intersect(mrna.keep,rpf.keep)
    mrna <- mrna[keep.genes,,drop=FALSE]
    rpf <- rpf[keep.genes,,drop=FALSE]

    ## if the colnames of mrna and rpf are same, add characters to distinguish them.
    if (sum(colnames(mrna) == colnames(rpf)) >= 1){
        colnames(mrna) = paste0("mrna_",colnames(mrna))
        colnames(rpf) = paste0("rpf_",colnames(rpf))
    }
    list(mrna = mrna, rpf = rpf)
}

#' @importFrom stats median
.estimate_size_factors_for_matrix <- function(counts,
                                              locfunc = stats::median){
    loggeomeans <- rowMeans(log(counts))
    if (all(is.infinite(loggeomeans))) {
        stop("every gene contains at least one zero, cannot compute log ",
             "geometric means",
             call. = FALSE)
    }
    sf <- apply(counts, 2, function(cnts) {
        exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
    })
    sf
}

.get_size_factors <- function(mrna, rpf, normalize){
    #normalize together
    if (normalize){
        message("Calculating the library size factors")
        pool_sizeFactor <- .estimate_size_factors_for_matrix(cbind(mrna,rpf))
        mrna_sizeFactor <- pool_sizeFactor[1:ncol(mrna)]
        rpf_sizeFactor <- pool_sizeFactor[(ncol(mrna)+1):(ncol(mrna)+ncol(rpf))]
    } else {
        mrna_sizeFactor <- rep(1, ncol(mrna))
        rpf_sizeFactor <- rep(1,ncol(rpf))
    }
    list(mrna_sizeFactor = mrna_sizeFactor, rpf_sizeFactor = rpf_sizeFactor)
}


#' @rdname xtail
#' @importFrom SummarizedExperiment assay
#' @export
xtail <- function(mrna,
                  rpf,
                  condition,
                  baseLevel = NA,
                  minMeanCount = 1,
                  normalize = TRUE,
                  p.adjust.method ="BH",
                  threads = NA,
                  bins = 10000L,
                  ci = 0){
    # input checks
    ## there must be exactly two different conditon
    if (length(unique(condition))!=2) {
        stop("There must be exactly two different conditions.",
             call. = FALSE)
    }
    # normalize must be TRUE or FALSE
    if(!is.logical(normalize) || length(normalize) != 1L || is.na(normalize)){
        stop("'normalize' must be TRUE or FALSE.", call. = FALSE)
    }
	## default baseLevel is the first coniditon
	if (is.na(baseLevel)) {
	    baseLevel <- condition[1]
	}
	## if the baseLevel is in conditon
	if (length(baseLevel) != 1L || !(baseLevel %in% condition)){
	    stop("baseLevel is not in condition", call. = FALSE)
	}
	## if the colnames of mrna and rpf are matched
    if(is.null(rownames(mrna)) || is.null(rownames(rpf))){
        stop("rownames of 'mrna' and 'rpf' must be set.",
             call. = FALSE)
    }
	if (ncol(mrna) != ncol(rpf)) {
	    stop("'mrna' and 'rpf' must have same number of columns.",
	         call. = FALSE)
	}
	## if the colnames and condition are match
	if(length(condition) != ncol(mrna)) {
	    stop("condition must have same length as the number of columns of rna.",
	         call. = FALSE)
	}
	## check confidence interval
	if (length(ci) != 1L || ci<0 || ci > 1 ) {
	    stop("ci must be a single numeric value between 0 and 1.",
	         call. = FALSE)
	}
	## check pvalue.adjust method
	supported.padj <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr")
	if (length(p.adjust.method) != 1L || !(p.adjust.method %in% supported.padj)) {
	    stop(paste0("The adjustment methods must be one of '",
	                paste(supported.padj,collapse="', '")),
	         "'.",
	         call. = FALSE)
	}
	if (length(minMeanCount) != 1L || minMeanCount < 1){
	    stop("minMeanCount must be single numeric value and be larger than 1.",
	         call. = FALSE)
	}
    # norm the input matrices and
	mrna_rpf_normed <- .norm_mrna_rpf(mrna, rpf, minMeanCount)
    mrna <- mrna_rpf_normed$mrna
    rpf <- mrna_rpf_normed$rpf
    size_factors <- .get_size_factors(mrna, rpf, normalize)
    mrna_sizeFactor <- size_factors$mrna_sizeFactor
    rpf_sizeFactor <- size_factors$rpf_sizeFactor
    if(!normalize){
        mrna <- round(mrna)
        rpf <- round(rpf)
    }
    ############################################################################
	## 1. Estimate the difference of log2FC between mRNA and RPF
	#If no replicate, we assume no more than one-third of genes' expression changed.
	message ("1. Estimate the log2 fold change in mrna")
	mrna_object = .estimate(mrna,
	                        condition,
	                        baseLevel,
	                        mrna_sizeFactor)
	message ("2. Estimate the log2 fold change in rpf")
	rpf_object = .estimate(rpf,
	                       condition,
	                       baseLevel,
	                       rpf_sizeFactor)
	message ("3. Estimate the difference between two log2 fold changes")
	result_log2FC = .xtail_test(mrna_object,
	                            rpf_object,
	                            threads,
	                            bins,
	                            baseLevel,
	                            ci)
	############################################################################
	### 2. Estimate the difference of log2R between control and treatment
	cond_base <- condition == baseLevel
	condition1_mrna <- mrna[,cond_base,drop=FALSE]
	condition1_rpf <- rpf[,cond_base,drop=FALSE]
	condition2_mrna <- mrna[,!cond_base,drop=FALSE]
	condition2_rpf <- rpf[,!cond_base,drop=FALSE]
	condition1_counts <- cbind(condition1_mrna,condition1_rpf)
	condition2_counts <- cbind(condition2_mrna,condition2_rpf)

	condition1_sizeFactor <- c(mrna_sizeFactor[cond_base],rpf_sizeFactor[cond_base])
	condition2_sizeFactor <- c(mrna_sizeFactor[!cond_base],rpf_sizeFactor[!cond_base])

	condition1_disper <- cbind(assay(mrna_object,XTAIL_DISPERSION_ASSAY)[,cond_base,drop=FALSE],
	                           assay(rpf_object,XTAIL_DISPERSION_ASSAY)[,cond_base,drop=FALSE])
	condition2_disper <- cbind(assay(mrna_object,XTAIL_DISPERSION_ASSAY)[,!cond_base,drop=FALSE],
	                           assay(rpf_object,XTAIL_DISPERSION_ASSAY)[,!cond_base,drop=FALSE])
	condition1 <- c(paste0(condition[cond_base],"_mRNA"),
	                paste0(condition[cond_base],"_rpf"))
	baseLevel1 <- paste0(baseLevel,"_mRNA")
	condition2 <- c(paste0(condition[!cond_base],"_mRNA"),
	                paste0(condition[!cond_base],"_rpf"))
	baseLevel2 <- paste0(unique(condition[!cond_base]),"_mRNA")
	#
	message ("4. Estimate the log2 ratio in first condition")
	condition1_object <- .estimate(condition1_counts,
	                               condition1,baseLevel1,
	                               condition1_sizeFactor,
	                               condition1_disper)
	message ("5. Estimate the log2 ratio in second condition")
	condition2_object <- .estimate(condition2_counts,
	                               condition2,baseLevel2,
	                               condition2_sizeFactor,
	                               condition2_disper)
	message ("6. Estimate the difference between two log2 ratios")
	result_log2R = .xtail_test(condition1_object,
	                           condition2_object,
	                           threads,
	                           bins,
	                           baseLevel,
	                           ci)
	############################################################################
	## 3. combine the log2FC and log2R results and report
	xtail <- .new_xtail(result_log2FC, result_log2R, condition, baseLevel, ci,
	                    p.adjust.method)
	nums <- resultsNum(xtail)
	message("Number of the log2FC and log2R used in determining the final p-value")
	message(paste0(" log2FC: ", nums["numFoldChange"]))
	message(paste0(" log2R: ", nums["numRatio"]))
	xtail
}

#' @importFrom locfit locfit
#' @importFrom stats predict
.get_xtail_dispersion_fit_function <- function(rawmeans,
                                               rawdisps,
                                               quantilePercent = 0.35,
                                               binnum = 50,
                                               minDisp = 1e-8){
  useForFit <- rawdisps > 100 * minDisp
  if(sum(useForFit) == 0){
        stop("all gene-wise dispersion estimates are within 2 orders of ",
             "magnitude from the minimum value, and so the gene-wise ",
             "estimates as final estimates.")
  }
  usedMeans <- rawmeans[useForFit]
  usedDisps <- rawdisps[useForFit]
  sortMeans <- usedMeans[order(usedMeans)]
  sortDisps <- usedDisps[order(usedMeans)]
  genenum <- length(sortMeans)
  quantile_means <- c()
  quantile_disps <- c()
  for(i in seq(1,genenum,binnum)){
    num = min(binnum,genenum-i+1)
    if(num<(binnum)){
      next
    }
    curbin_means <- sortMeans[i:(i+num-1)]
    curbin_disps <- sortDisps[i:(i+num-1)]

    m <- curbin_means[order(curbin_disps)][1:floor(num*quantilePercent)]
    d <- curbin_disps[order(curbin_disps)][1:floor(num*quantilePercent)]
    quantile_means <- c(quantile_means, m)
    quantile_disps <- c(quantile_disps, d)
  }
  dat <- data.frame(logDisps = log(quantile_disps), logMeans=log(quantile_means))
  fit <- locfit(logDisps~logMeans, data=dat[quantile_disps>=minDisp*10,,drop=FALSE])
  dispFit <- function(rawmeans){
      exp(predict(fit,data.frame(logMeans=log(rawmeans))))
  }
  dispFit
}

#' @import DESeq2
#' @importFrom stats relevel formula
#' @importFrom S4Vectors mcols
#' @importFrom SummarizedExperiment assay assay<-
.estimate <- function(countData, condition, baseLevel, libsize, dispers = NULL){
    ## using the DESeqDataSet to store data
    colData <- data.frame(row.names = colnames(countData),
                          condition = relevel(factor(condition), baseLevel))
    dataSet <- DESeqDataSetFromMatrix(countData = countData,
                                      colData = colData,
                                      design = ~ condition)
    # colData(dataSet)$condition <- relevel(colData(dataSet)$condition, baseLevel)

    # normalization
    if (missing(libsize)){
        dataSet <- suppressMessages(estimateSizeFactors(dataSet))
    } else {
        sizeFactors(dataSet) <- libsize
    }

    # estimateDispersion
    if (is.null(dispers)){
        if (ncol(countData) == 2){
            object <- dataSet #creat a copy of a dataSet for estimating dispersion
            design(object) <- formula(~1) # take two conditon samples as replicates for estimating dispersion
            object <- estimateDispersionsGeneEst(object)
            #fit
            dispFitFun <- .get_xtail_dispersion_fit_function(
                mcols(object)$baseMean,
                mcols(object)$dispGeneEst)
            #attr(dispFitFun, "fitType") = "bin quantile"
            fittedDisps <- dispFitFun(mcols(object)$baseMean)
            dispers <- matrix(rep(fittedDisps, ncol(dataSet)),
                              ncol = ncol(dataSet),
                              byrow = FALSE)
        } else {
            dataSet <- suppressMessages(estimateDispersions(dataSet))
            dispers <- matrix(rep(dispersions(dataSet), ncol(dataSet)),
                              ncol = ncol(dataSet),
                              byrow = FALSE)
        }
    }
    assay(dataSet, XTAIL_DISPERSION_ASSAY, withDimnames = FALSE) <- dispers
    # estimate log2FC or log2R
    dataSet <- suppressMessages(
        .estimate_MLE_for_Beta(dataSet, modelMatrixType = "standard"))
    dataSet
}


#' @importFrom parallel makeCluster clusterExport detectCores parLapply
#'   stopCluster
#' @importFrom S4Vectors mcols
#' @importFrom SummarizedExperiment assay colData
.xtail_test <- function(object1, object2,threads,bins,baseLevel, ci){
    intersect.genes <- intersect(rownames(object1), rownames(object2))
    object1 <- object1[intersect.genes,,drop=FALSE]
    object2 <- object2[intersect.genes,,drop=FALSE]
    counts1 <- counts(object1)
    counts1[which(counts1==0)] = 1
    counts2 <- counts(object2)
    counts2[which(counts2==0)] = 1
    intercept1 <- mcols(object1)$Intercept
    intercept2 <- mcols(object2)$Intercept
    dispersion1 <- assay(object1,XTAIL_DISPERSION_ASSAY)
    dispersion2 <- assay(object2,XTAIL_DISPERSION_ASSAY)
    resName1 <- resultsNames(object1)[2]
    resName2 <- resultsNames(object2)[2]
    log2Ratio1 <- mcols(object1)[[resName1]]
    log2Ratio2 <- mcols(object2)[[resName2]]
    sizefactor1 <- sizeFactors(object1)
    sizefactor2 <- sizeFactors(object2)
    cond1 <- as.numeric(colData(object1)$condition) - 1L
    cond2 <- as.numeric(colData(object2)$condition) - 1L
    ## Estimate the probabilities of log2R and log2FC
    ## parallel used for speeding up
    rowNo <- nrow(counts1)
    #
    if(is.na(threads)){
        ## automaticly detect the number of CPU cores.
        cluster <- makeCluster(detectCores())
    } else {
        cluster <- makeCluster(threads)
    }
    on.exit(stopCluster(cluster))
    clusterExport(cl = cluster,
                  varlist = c("counts1","counts2","intercept1","intercept2",
                              "log2Ratio1","log2Ratio2","dispersion1",
                              "dispersion2","sizefactor1","sizefactor2","cond1",
                              "cond2","bins","ci"),
                  envir = environment())
    res <- parLapply(cluster,
                     X = seq_len(rowNo),
                     fun = .xtail_test_wrapper)
    res <- matrix(unlist(res),
                  ncol = 4L,
                  byrow = TRUE,
                  dimnames = list(intersect.genes,
                                  c("deltaTE","Pval","lowci","highci")))
    res <- data.frame(res)
    res$log2Ratio1 <- log2Ratio1
    res$log2Ratio2 <- log2Ratio2
    res
}
