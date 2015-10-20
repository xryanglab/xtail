\name{xtail}
\alias{xtail}
\docType{package}
\title{
A tool to analyze differential translation.
}
\description{
From paird mRNA and RPF count data, run xtail to asessment of translations.
}
\usage{
xtail(mrna,rpf,condition,baseLevel=NA,minMeanCount=10,method.adjust="BH",
threadsNo=NA,bins=10000)
}
\arguments{
  \item{mrna}{a matrix or data frame of raw (not normalized) mRNA count data whose rows correspond to genes and columns correspond to samples.}
  \item{rpf}{a matrix or data frame of raw (not normalized) RPF count data whose rows correspond to genes and columns correspond to samples.}
  \item{condition}{condition lavels corresponding to the samples in mrna and rpf. There must be exactly two unique values.}
  \item{baseLevel}{The baseLevel is a label that assign the condition be compared. If \code{xtail} is run without specifying baseLevel, it will return the comparison of the second condition over the first condition.}
  \item{minMeanCount}{\code{xtail} uses the average expression level of each gene, across all samples as filter criterion and it omits all genes with mean counts below minMeanCount}
  \item{method.adjust}{The method to use for adjusting multiple comparisons, see \code{?p.adjust} }
  \item{threadsNo}{The number of CPU cores used. By default, all available  cores are used.}
  \item{bins}{The number of bins used for calculate the probability density of log2FC or log2R. This paramater will determine the sampling probability density function.Set it small for a very quick test run.}
}
\details{
  
  The step of estimation of the probability distribution, for log2FC and log2R, will execute slowly in the current implementation, but can be speeded up by running on multiple cores using the parallel library. By default, the 'detectCores' function in parallel library is used to determine the number of CPU cores in the machine on which R is running. To adjust the number of cores used, use 'threadsNo' argument to assign.
}
\examples{
  data(testdata)
  #Get the mrna count data and rpf count data
  test.mrna <- testdata$mrna
  test.rpf <- testdata$rpf
  #Assign condition labels to samples. 
  condition <- c("control","control","treatment","treatment")
  #run xtail 
  test.results <- xtail(test.mrna,test.rpf,condition)
  #save the results
  write.table(test.results,"test_results.txt",quote=F,sep="\t")
}
\author{
Zhengtao xiao
}
\references{
Zhengtao Xiao, Qin Zou, Yu Liu, and Xuerui Yang: Genome-wide assessment of differential translations with ribosome profiling data.
}
