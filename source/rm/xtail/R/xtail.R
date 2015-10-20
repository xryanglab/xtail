#'	@useDynLib xtail
doTest <- function(object1, object2,threadsNo,nbins,baseLevel){
	intersectionGenes <- intersect(rownames(object1), rownames(object2))
	object1 <- object1[intersectionGenes,]
	object2 <- object2[intersectionGenes,]
	counts1 <- counts(object1)
	counts2 <- counts(object2)
	baseCondition <- paste("condition",baseLevel,sep="")
	intercept1 <- mcols(object1)$Intercept
	intercept2 <- mcols(object2)$Intercept
	dispersion1 <- dispersions(object1)
	dispersion2 <- dispersions(object2)
	resNms1 = resultsNames(object1)
	resNms2 = resultsNames(object2)
	log2Ratio1 <- mcols(object1)[[resNms1[length(resNms1)]]]
	log2Ratio2 <- mcols(object2)[[resNms2[length(resNms2)]]]
	sizefactor1 <- sizeFactors(object1)
	sizefactor2 <- sizeFactors(object2)
	betaPriorVar1 <- attr(object1, "betaPriorVar")
	betaPriorVar2 <- attr(object2, "betaPriorVar")
	priorSigma1 <- sqrt(betaPriorVar1[2]) * sqrt(2)
	priorSigma2 <- sqrt(betaPriorVar2[2]) * sqrt(2)
	cond1 <- as.numeric(colData(object1)$condition) - 1
	cond2 <- as.numeric(colData(object2)$condition) - 1
	
	## Estimate the probabilities of log2R and log2FC
	## parallel used for speeding up


	rowNo <- nrow(counts1)
	if (is.na(threadsNo)){
		cluster <- makeCluster(detectCores())
	}else{
		cluster <- makeCluster(threadsNo)
	}
	#####!!!!!!!# cpp


	clusterExport(cluster, c('counts1','counts2','intercept1','intercept2','log2Ratio1','log2Ratio2','priorSigma1','priorSigma2','dispersion1','dispersion2','sizefactor1','sizefactor2','cond1','cond2','nbins'),envir=environment())
	res <- clusterMap(cluster, doTestWrapper, i = c(1:rowNo));
	stopCluster(cluster)
	res = unlist(res)
	tmp = length(res)
	ovl = res[1:tmp %% 2 == 1]
	pvl = res[1:tmp %% 2 == 0]
	res = cbind(log2Ratio1,log2Ratio2,ovl,pvl)
	rownames(res) = rownames(object1)
	res 
}
doTestWrapper <- function(i){
	res <- xtail_test(counts1[i,],counts2[i,],intercept1[i],intercept2[i],log2Ratio1[i],log2Ratio2[i],priorSigma1,priorSigma2,dispersion1[i],dispersion2[i],sizefactor1,sizefactor2,cond1,cond2,nbins)
}
callDESeq2 <- function(data, condition, baseLevel){
	##using the DESeqDataSet to store the counts
	colData <- data.frame(row.names = colnames(data), condition=condition)
	dataSet <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~condition)
	colData(dataSet)$condition <- relevel(colData(dataSet)$condition, baseLevel)


	# normalization
	message ("normalize counts by size factors") 
	dataSet <- suppressMessages(estimateSizeFactors(dataSet))
	# estimateDispersion
	message ("estimate dispersion paramater of NB model")
	dataSet <- suppressMessages(estimateDispersions(dataSet))
	# estimatelog2FC
	message("estimate the log2 fold change")
	dataSet <- suppressMessages(nbinomWaldTest(dataSet,modelMatrixType="standard"))
	dataSet
}
#' @export
xtail <- function(mrna, rpf, condition, baseLevel = NA, minMeanCount = 10, method.adjust="BH", threadsNo=NA,bins=10000)
{
	## default baseLevel is the first coniditon
	if (is.na(baseLevel)) {baseLevel <- condition[1]}
	## if the baseLevel is in conditon
	if (!is.element(baseLevel, condition)) stop("baseLevel is not in mrna_condition")
	## if the colnames of mrna and rpf are matched
	if (dim(mrna)[2] != dim(rpf)[2]) stop("the mrna and rpf must have same number of columns")
	## if the condition must have the same length as the number of the columns of mrna
	if (length(condition)!=dim(mrna)[2]) stop("condition must have same length as the number of columns of mrna")
	## ther must be exactly two different conditon
	if (length(unique(condition))!=2) stop("There must be eaxctly two different condition")
	## filter genes with low count number
	if (minMeanCount<1) stop("minMeanCount needs to be at least 1")
	keep <- rowMeans(mrna) >= minMeanCount
	mrna_keep <- rownames(mrna)[keep]
	keep <- rowMeans(rpf) >= minMeanCount
	rpf_keep <- rownames(rpf)[keep]
	## merge genes in mrna and rpf
	keep_genes <- intersect(mrna_keep,rpf_keep)
	mrna <- mrna[keep_genes,]
	rpf <- rpf[keep_genes,]

	## if the colnames of mrna and rpf are same, add characters to distinguish them .
	if (sum(colnames(mrna) == colnames(rpf)) >= 1){
		colnames(mrna) = paste("mrna_",colnames(mrna),sep="")
		colnames(rpf) = paste("rpf_",colnames(rpf),sep="")
	}
	#else{
	#	stop ("colnames of mrna and rpf must match")
	#}

	## 1. Estimate the difference of log2FC between mRNA and RPF
	## normalization, estimate dispersion of NB model, estimate the log2 fold change
	message ("1. Estimate the log2 fold change in mrna")
	mrna_object = callDESeq2(mrna,condition,baseLevel)
	message ("2. Estimate the log2 fold change in rpf")
	rpf_object = callDESeq2(rpf,condition,baseLevel)
	message ("3. Estimate the log2FC difference between mrna and rpf")
	result_log2FC = doTest(mrna_object,rpf_object,threadsNo,bins,baseLevel)

	### 2. Estimate the difference of log2R between control and treatment 
	## normalization, estimate dispersion of NB model, estimate the log2 fold change
	condition1_mrna <- mrna[,condition==baseLevel]
	condition1_rpf <- rpf[,condition==baseLevel]
	condition2_mrna <- mrna[,condition!=baseLevel]
	condition2_rpf <- rpf[,condition!=baseLevel]
	condition1_counts <- cbind(condition1_mrna,condition1_rpf)
	condition2_counts <- cbind(condition2_mrna,condition2_rpf)
	##
	message ("4. Estimate the log2 ratio in first condition")
	condition1_object <- callDESeq2(condition1_counts,condition,baseLevel)
	message ("5. Estimate the log2 ratio in second condition")
	condition2_object <- callDESeq2(condition2_counts,condition,baseLevel)
	message ("6. Estimate the log2R difference between two conditions")
	result_log2R = doTest(condition1_object,condition2_object,threadsNo,bins,baseLevel)	

	## 3. combine the log2FC and log2R results and report
	###!!!!
	intersectionGenes <- intersect(rownames(result_log2FC), rownames(result_log2R))
	result_log2FC <- result_log2FC[intersectionGenes,]
	result_log2R <- result_log2R[intersectionGenes,]
	result_combine = cbind(result_log2FC,result_log2R)
	baseLevelName <- paste(baseLevel,"_log2R",sep="")
	secondLevelName <- paste(unique(condition[condition!=baseLevel]), "_log2R",sep="")
	colnames(result_combine) <- c("mRNA_log2FC", "RPF_log2FC", "OVL_log2FC","pvalue_log2FC",baseLevelName,secondLevelName,"OVL_log2R","pvalue_log2R")
	result_combine <- as.data.frame(result_combine)
	pv_log2FC <- result_combine$pvalue_log2FC
	pv_log2R <- result_combine$pvalue_log2R
	log2FC_diff <- result_combine$RPF_log2FC - result_combine$mRNA_log2FC
	log2R_diff <- result_combine[[baseLevelName]] - result_combine[[secondLevelName]]

	final_pvalues <- rep(0,nrow(mrna))
	final_diff <- rep(0,nrow(mrna))
	for (i in 1:nrow(mrna)){
		if (is.na(pv_log2FC[i]) || is.na(pv_log2R[i])){
			final_pvalues[i] <- NA
			final_diff[i] <- NA
		}else if(pv_log2FC[i] > pv_log2R[i]){
			final_pvalues[i] <- pv_log2FC[i]
			final_diff[i] <- log2FC_diff[i]
		}else{
			final_pvalues[i] <- pv_log2R[i]
			final_diff[i] <- log2R_diff[i]
		}
	}

	result_combine$final_pvalues <- final_pvalues
	fdr <- p.adjust(final_pvalues, method = method.adjust)
	result_combine$FDR = fdr
	# result_combine[["diff Ratios"]] = final_diff
	result_combine
}
