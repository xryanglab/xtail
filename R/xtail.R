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

#' @export
xtail <- function(mrna, rpf, condition, baseLevel = NA, minMeanCount = 1, normalize = TRUE, method.adjust="BH", threads=NA,bins=10000,ci = 0)
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
	#so only ... were used for fit dispersion.
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
	condition1_disper <- cbind(dispersionMatrix(mrna_object)[,1:ncol(condition1_mrna),drop=FALSE], dispersionMatrix(rpf_object)[,1:ncol(condition1_rpf),drop=FALSE])
	condition2_disper <- cbind(dispersionMatrix(mrna_object)[,(ncol(condition1_mrna)+1):ncol(mrna),drop=FALSE], dispersionMatrix(rpf_object)[,(ncol(condition1_rpf)+1):ncol(rpf),drop=FALSE])

	#
	message ("4. Estimate the log2 ratio in first condition")
	condition1_object <- estimateFun(condition1_counts,condition,baseLevel,condition1_sizeFactor,condition1_disper)
	message ("5. Estimate the log2 ratio in second condition")
	condition2_object <- estimateFun(condition2_counts,condition,baseLevel,condition2_sizeFactor,condition2_disper)
	message ("6. Estimate the difference between two log2 ratios")
	result_log2R = xTest(condition1_object,condition2_object,threads,bins,baseLevel,ci)

	## 3. combine the log2FC and log2R results and report
	###!!!!
	intersect.genes <- intersect(rownames(result_log2FC), rownames(result_log2R))
	result_log2R <- result_log2R[intersect.genes,]
	result_log2FC <- result_log2FC[intersect.genes,]

	##!!!return ("deltaTE","Pval","lowci","highci")
	final_result <- cbind(result_log2FC[,1:2],result_log2R[,1:2])
	colnames(final_result) <- c("log2FC_TE_v1","pvalue_v1","log2FC_TE_v2","pvalue_v2")
	final_result <- as.data.frame(final_result)

	final_result$pvalue_final <- rep(0,nrow(final_result))
	final_result$pvalue.adjust <- rep(0,nrow(final_result))
	final_result$log2FC_TE_final <- rep(0,nrow(final_result))
	final_result_CI <- rep(0,nrow(final_result))
	log2FC_determine <- 0
	log2R_determine <- 0
	for (i in 1:nrow(final_result)){
		if (is.na(final_result[i,2]) || is.na(final_result[i,4])){
			final_result$pvalue_final[i] <- NA
			final_result$log2FC_TE_final[i] <- NA
			final_result_CI[i] <- NA
		}else if(final_result[i,2] > final_result[i,4]){
			final_result$pvalue_final[i] <- final_result[i,2]
			final_result$log2FC_TE_final[i] <- final_result[i,1]
			final_result_CI[i] <- paste0("[",round(result_log2FC[i,3],2),",",round(result_log2FC[i,4],2),"]")
			log2FC_determine <- log2FC_determine + 1
		}else{
			final_result$pvalue_final[i] <- final_result[i,4]
			final_result$log2FC_TE_final[i] <- final_result[i,3]
			final_result_CI[i] <- paste0("[",round(result_log2R[i,3],2),",",round(result_log2R[i,4],2),"]")
			log2R_determine <- log2R_determine + 1
		}
	}
	final_result$pvalue.adjust = p.adjust(final_result$pvalue_final,method=method.adjust)
	if (ci>0){
		CI_string = paste0("CI(",100*ci,"%)")
		final_result[[CI_string]] = final_result_CI
	}
	message("Number of times the log2FC and log2R used in determining the final p-value")
	message(paste0(" log2FC: ", log2FC_determine))
	message(paste0(" log2R: ", log2R_determine))

	final_result
}
