#' @name xtail
#'
#' @author Zhengtao xiao
#'
#' @title A tool to quantitatively assess differential translations with ribosome profiling data.
#'
#' @description By pairwise comparisons of ribosome profiling data, Xtail
#' identifies differentially translated genes across two experimental or
#' physiological conditions.
#'
#' @details No missing values are allowed in input data mrna and rpf.
#'
#'   Duplicate row names (gene names or gene ids) are not allowed.
#'
#'   \code{Xtail} takes in raw read counts of RPF and mRNA, and performs
#'   median-of-ratios normalization. Alternatively, users can provide normalized
#'   read counts and skip the built-in normal by setting "normalize" to FALSE.
#'
#'   The step of estimation of the probability distributions, for log2FC or
#'   log2R, will execute slowly in the current implementation, but can be
#'   speeded up by running on multiple cores using the parallel library. By
#'   default, the "detectCores" function in parallel library is used to
#'   determine the number of CPU cores in the machine on which R is running. To
#'   adjust the number of cores used, use "threads" argument to assign.
#'
#' @seealso The \code{RNAmodR.RiboMethSeq} and \code{RNAmodR.AlkAnilineSeq}
#'   package.
#'
#' @references Zhengtao Xiao, Qin Zou, Yu Liu, and Xuerui Yang: Genome-wide
#'   assessment of differential translations with ribosome profiling data.
#'
#' @docType package
#'
#' @usage xtail(mrna,rpf,condition,baseLevel=NA,minMeanCount=1,...)
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
#' @param method.adjust The method to use for adjusting multiple comparisons, by
#'   default "BH",  see \code{?p.adjust}
#' @param threads The number of CPU cores used. By default, all available  cores
#'   are used.
#' @param bins The number of bins used for calculating the probability density
#'   of log2FC or log2R (default is 10000). This paramater will determine
#'   accuracy of pvalue. Set it small for a very quick test run.
#'
#' @examples
#' #load the data
#' data(xtaildata)
#' #Get the mrna count data and rpf count data
#' test.mrna <- xtaildata$mrna
#' test.rpf <- xtaildata$rpf
#' #Assign condition labels to samples.
#' condition <- c("control","control","treat","treat")
#' #run xtail
#' test.results <- xtail(test.mrna,test.rpf,condition)
NULL

#'	@useDynLib xtail
xTest <- function(object1, object2,threads,bins,baseLevel, ci){
	intersect.genes <- intersect(rownames(object1), rownames(object2))
	object1 <- object1[intersect.genes,,drop=FALSE]
	object2 <- object2[intersect.genes,,drop=FALSE]
	counts1 <- counts(object1)
	counts1[which(counts1==0)] = 1
	counts2 <- counts(object2)
	counts2[which(counts2==0)] = 1
	intercept1 <- mcols(object1)$Intercept
	intercept2 <- mcols(object2)$Intercept
	dispersion1 <- dispersionMatrix(object1)
	dispersion2 <- dispersionMatrix(object2)
	resName1 <- resultsNames(object1)[2]
	resName2 <- resultsNames(object2)[2]
	log2Ratio1 <- mcols(object1)[[resName1]]
	log2Ratio2 <- mcols(object2)[[resName2]]
	sizefactor1 <- sizeFactors(object1)
	sizefactor2 <- sizeFactors(object2)
	cond1 <- as.numeric(colData(object1)$condition) - 1
	cond2 <- as.numeric(colData(object2)$condition) - 1
	## Estimate the probabilities of log2R and log2FC
	## parallel used for speeding up
	rowNo <- nrow(counts1)
	if (is.na(threads)){
		## automaticly detect the number of CPU cores.
		cluster <- makeCluster(detectCores())
	}else{
		cluster <- makeCluster(threads)
	}
	clusterExport(cluster, c('counts1','counts2','intercept1','intercept2','log2Ratio1','log2Ratio2','dispersion1','dispersion2','sizefactor1','sizefactor2','cond1','cond2','bins','ci'),envir=environment())
	res <- clusterMap(cluster, xTestWrapper, i = c(1:rowNo));
	stopCluster(cluster)
	res <- matrix(unlist(res),ncol=4,byrow=T,dimnames=list(intersect.genes, c("deltaTE","Pval","lowci","highci")))
	res <- data.frame(res)
	res$log2Ratio1 <- log2Ratio1
	res$log2Ratio2 <- log2Ratio2
	res
}


xTestWrapper <- function(i){
	res <- xtail_test(counts1[i,],counts2[i,],intercept1[i],intercept2[i],log2Ratio1[i],log2Ratio2[i],dispersion1[i,],dispersion2[i,],sizefactor1,sizefactor2,cond1,cond2,bins,ci)
}

estimateFun <- function(countData, condition, baseLevel, libsize, dispers){
	## using the DESeqDataSet to store data
	colData <- data.frame(row.names = colnames(countData), condition=condition)
	dataSet <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~condition)
	colData(dataSet)$condition <- relevel(colData(dataSet)$condition, baseLevel)

	# normalization
	if (missing(libsize)){
		dataSet <- suppressMessages(estimateSizeFactors(dataSet))
	}else{
		sizeFactors(dataSet) <- libsize
	}

	# estimateDispersion
	if (missing(dispers))
	{
		if (ncol(countData) == 2){
			object = dataSet #creat a copy of a dataSet for estimating dispersion
			design(object) = formula(~1) # take two conditon samples as replicates for estimating dispersion
			object = estimateDispersionsGeneEst(object)
			#fit
			dispFitFun <- xtailDispersionFit(mcols(object)$baseMean, mcols(object)$dispGeneEst)
			#attr(dispFitFun, "fitType") = "bin quantile"
			fittedDisps <- dispFitFun(mcols(object)$baseMean)
			dispersionMatrix(dataSet) = matrix(rep(fittedDisps, ncol(dataSet)), ncol=ncol(dataSet), byrow=FALSE)
		}else{
			dataSet <- suppressMessages(estimateDispersions(dataSet))
			dispersionMatrix(dataSet) <- matrix(rep(dispersions(dataSet), ncol(dataSet)),ncol=ncol(dataSet), byrow=FALSE)
		}
	}else{
		dispersionMatrix(dataSet) <- dispers
	}

	# estimate log2FC or log2R
	dataSet <- suppressMessages(estimateMLEForBeta(dataSet, modelMatrixType="standard"))
	dataSet
}

#' @rdname xtail
#' @export
xtail <- function(mrna, rpf, condition, baseLevel = NA, minMeanCount = 1, normalize = TRUE, p.adjust.method ="BH", threads=NA,bins=10000,ci = 0)
{
	## default baseLevel is the first coniditon
	if (is.na(baseLevel)) {baseLevel <- condition[1]}
	## if the baseLevel is in conditon
	if (!is.element(baseLevel, condition)) stop("baseLevel is not in condition")
	## if the colnames of mrna and rpf are matched
	if (ncol(mrna) != ncol(rpf)) stop("the mrna and rpf must have same number of columns")
	## if the colnames and condition are match
	if(length(condition)!=ncol(mrna)) stop("condition must have same length as the number of columns of rna")
	## ther must be exactly two different conditon
	if (length(unique(condition))!=2) stop("There must be exactly two different conditions")
	## check confidence interval
	if (ci<0 ) stop("ci must be non-negative")
	## check pvalue.adjust method
	supported.padj <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr")
	if (!is.element(p.adjust.method, supported.padj)) stop(paste0("The adjustment methods must be one of ",paste(supported.padj,collapse=", ")))

	## filter genes with low count number
	if (minMeanCount<1) stop("minMeanCount needs to be at least 1")
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
	#normalize together
	if (normalize){
		message("Calculating the library size factors")
		pool_sizeFactor <- estimateSizeFactorsForMatrix(cbind(mrna,rpf))
		mrna_sizeFactor <- pool_sizeFactor[1:ncol(mrna)]
		rpf_sizeFactor <- pool_sizeFactor[(ncol(mrna)+1):(ncol(mrna)+ncol(rpf))]
	}else{
		mrna <- round(mrna)
		rpf <- round(rpf)
		mrna_sizeFactor <- rep(1, ncol(mrna))
		rpf_sizeFactor <- rep(1,ncol(rpf))
	}

	## 1. Estimate the difference of log2FC between mRNA and RPF
	#If no replicate, we assume no more than one-third of genes' expression changed.
	message ("1. Estimate the log2 fold change in mrna")
	mrna_object = estimateFun(mrna,condition,baseLevel,mrna_sizeFactor)
	message ("2. Estimate the log2 fold change in rpf")
	rpf_object = estimateFun(rpf,condition,baseLevel,rpf_sizeFactor)
	message ("3. Estimate the difference between two log2 fold changes")
	result_log2FC = xTest(mrna_object,rpf_object,threads,bins,baseLevel,ci)

	### 2. Estimate the difference of log2R between control and treatment
	condition1_mrna <- mrna[,condition==baseLevel,drop=FALSE]
	condition1_rpf <- rpf[,condition==baseLevel,drop=FALSE]
	condition2_mrna <- mrna[,condition!=baseLevel,drop=FALSE]
	condition2_rpf <- rpf[,condition!=baseLevel,drop=FALSE]
	condition1_counts <- cbind(condition1_mrna,condition1_rpf)
	condition2_counts <- cbind(condition2_mrna,condition2_rpf)

	condition1_sizeFactor <- c(mrna_sizeFactor[condition==baseLevel],rpf_sizeFactor[condition==baseLevel])
	condition2_sizeFactor <- c(mrna_sizeFactor[condition!=baseLevel],rpf_sizeFactor[condition!=baseLevel])

	condition1_disper <- cbind(dispersionMatrix(mrna_object)[,condition==baseLevel,drop=FALSE], dispersionMatrix(rpf_object)[,condition==baseLevel,drop=FALSE])
	condition2_disper <- cbind(dispersionMatrix(mrna_object)[,condition!=baseLevel,drop=FALSE], dispersionMatrix(rpf_object)[,condition!=baseLevel,drop=FALSE])
	condition1 <- c(paste0(condition[condition == baseLevel],"_mRNA"), paste0(condition[condition == baseLevel],"_rpf"))
	baseLevel1 <- paste0(baseLevel,"_mRNA")
	condition2 <- c(paste0(condition[condition != baseLevel],"_mRNA"), paste0(condition[condition != baseLevel],"_rpf"))
	baseLevel2 <- paste0(unique(condition)[2],"_mRNA")
	#
	message ("4. Estimate the log2 ratio in first condition")
	condition1_object <- estimateFun(condition1_counts,condition1,baseLevel1,condition1_sizeFactor,condition1_disper)
	message ("5. Estimate the log2 ratio in second condition")
	condition2_object <- estimateFun(condition2_counts,condition2,baseLevel2,condition2_sizeFactor,condition2_disper)
	message ("6. Estimate the difference between two log2 ratios")
	result_log2R = xTest(condition1_object,condition2_object,threads,bins,baseLevel,ci)

	## 3. combine the log2FC and log2R results and report

	intersect.genes <- intersect(rownames(result_log2FC), rownames(result_log2R))
	result_log2R <- result_log2R[intersect.genes,]
	result_log2FC <- result_log2FC[intersect.genes,]

	#result data frame
	condition1_TE <- paste0(baseLevel,"_log2TE")
	condition2_TE <- paste0(unique(condition)[2], "_log2TE")
	final_result <- cbind(result_log2FC[,c("log2Ratio1","log2Ratio2","deltaTE","Pval")],result_log2R[,c("log2Ratio1","log2Ratio2","deltaTE","Pval")])
	colnames(final_result) <- c("mRNA_log2FC","RPF_log2FC","log2FC_TE_v1","pvalue_v1",condition1_TE,condition2_TE,"log2FC_TE_v2","pvalue_v2")
	final_result <- as.data.frame(final_result)
	final_result$log2FC_TE_final <- 0
	final_result$pvalue_final <- 0
	final_result$pvalue.adjust <- 0
	if (ci>0){
		resultCI <- NA
	}
	log2FC_determine_num <- 0
	log2R_determine_num <- 0
	for (i in 1:nrow(final_result)){
		if (is.na(final_result[i,"pvalue_v1"]) || is.na(final_result[i,"pvalue_v2"])){
		  final_result$log2FC_TE_final[i] <- NA
		  final_result$pvalue_final[i] <- NA
			if (ci>0) resultCI[i] <- NA
		}else if(final_result[i,"pvalue_v1"] > final_result[i,"pvalue_v2"]){
		  final_result$log2FC_TE_final[i] <- final_result[i,"log2FC_TE_v1"]
			final_result$pvalue_final[i] <- final_result[i,"pvalue_v1"]
			if (ci>0) resultCI[i] <- paste0("[",round(result_log2FC[i,"lowci"],2),",",round(result_log2FC[i,"highci"],2),"]")
			log2FC_determine_num <- log2FC_determine_num + 1
		}else{
		  final_result$log2FC_TE_final[i] <- final_result[i,"log2FC_TE_v2"]
		  final_result$pvalue_final[i] <- final_result[i,"pvalue_v2"]
			if (ci>0) resultCI[i] <- paste0("[",round(result_log2R[i,"lowci"],2),",",round(result_log2R[i,"highci"],2),"]")
			log2R_determine_num <- log2R_determine_num + 1
		}
	}

	final_result$pvalue.adjust = p.adjust(final_result$pvalue_final,method=p.adjust.method)
	if (ci>0){
		CI_string = paste0("CI(",100*ci,"%)")
		final_result[[CI_string]] = resultCI
	}
	message("Number of the log2FC and log2R used in determining the final p-value")
	message(paste0(" log2FC: ", log2FC_determine_num))
	message(paste0(" log2R: ", log2R_determine_num))
	xtail_results <- list(resultsTable = final_result, log2FC_determine_num = log2FC_determine_num,
						log2R_determine_num=log2R_determine_num,condition1=baseLevel,
						condition2=unique(condition)[2] )
	class(xtail_results) <- c("xtailResults","list")
	xtail_results
}
