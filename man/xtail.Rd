\name{Xtail}
\alias{xtail}
\docType{package}
\title{
A tool to quantitatively assess differential translations with ribosome profiling data.
}
\description{
By pairwise comparisons of ribosome profiling data, Xtail identifies differentially translated genes across two experimental or physiological conditions.
}
\usage{
xtail(mrna,rpf,condition,baseLevel=NA,minMeanCount=1,...)
}
\arguments{
  \item{mrna}{a matrix or data frame of raw mRNA count data whose rows correspond to genes and columns correspond to samples. The column names should be non-empty, and in same order with condition.}
  \item{rpf}{a matrix or data frame of raw RPF count data whose rows correspond to genes and columns correspond to samples.The column names should be non-empty, and in same order with condition.}
  \item{condition}{condition labels corresponding to the order of samples in mrna and rpf. There must be exactly two unique values.}
  \item{baseLevel}{The baseLevel indicates which one of the two conditions will be compared against by the other one. If not specified, \code{Xtail} will return results of comparing the second condition over the first one.}
  \item{minMeanCount}{\code{Xtail} uses the average expression level of each gene, across all samples as filter criterion and it omits all genes with mean counts below minMeanCount.}
  \item{ci}{The level of confindence to get credible intervals of log2 fold change of translational efficiency (TE), for example 0.95.}
  \item{normalize}{Whether normalization should be done (TRUE \\ FALSE). If missing, \code{Xtail} will perform median-of-ratios normlazation by default.}
  \item{method.adjust}{The method to use for adjusting multiple comparisons, by default "BH",  see \code{?p.adjust} }
  \item{threads}{The number of CPU cores used. By default, all available  cores are used.}
  \item{bins}{The number of bins used for calculating the probability density of log2FC or log2R (default is 10000). This paramater will determine accuracy of pvalue. Set it small for a very quick test run.}
}
\details{
  No missing values are allowed in input data mrna and rpf.

  Duplicate row names (gene names or gene ids) are not allowed.

 \code{Xtail} takes in raw read counts of RPF and mRNA, and performs median-of-ratios normalization. Alternatively, users can provide normalized read counts and skip the built-in normal by setting "normalize" to FALSE.

  The step of estimation of the probability distributions, for log2FC or log2R, will execute slowly in the current implementation, but can be speeded up by running on multiple cores using the parallel library. By default, the "detectCores" function in parallel library is used to determine the number of CPU cores in the machine on which R is running. To adjust the number of cores used, use "threads" argument to assign.
}
\examples{
  #load the data
  data(xtaildata)
  #Get the mrna count data and rpf count data
  test.mrna <- xtaildata$mrna
  test.rpf <- xtaildata$rpf
  #Assign condition labels to samples.
  condition <- c("control","control","treat","treat")
  #run xtail
  test.results <- xtail(test.mrna,test.rpf,condition)
}
\author{
Zhengtao xiao
}
\references{
Zhengtao Xiao, Qin Zou, Yu Liu, and Xuerui Yang: Genome-wide assessment of differential translations with ribosome profiling data.
}
