## ----knitr, echo=FALSE---------------------------------------------------
library("knitr")

## ----<style,eval=TRUE,echo=FALSE,results='asis'-----------------------------------------
BiocStyle::latex()

## ----include=FALSE----------------------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----begain,results="hold",message=FALSE------------------------------------------------
library(xtail)
data(xtaildata)

## ---------------------------------------------------------------------------------------
mrna <- xtaildata$mrna
rpf <- xtaildata$rpf
head(mrna,5)
head(rpf,5)

## ---------------------------------------------------------------------------------------
condition <- c("control","control","treat","treat")

## ---------------------------------------------------------------------------------------
test.results <- xtail(mrna,rpf,condition,bins=1000)

## ----inspectData,echo=TRUE--------------------------------------------------------------
head(test.results,5)

## ----writeResult,eval=FALSE-------------------------------------------------------------
#  write.table(test.results,"test_results.txt",quote=F,sep="\t")

## ----sessInfo---------------------------------------------------------------------------
sessionInfo()

