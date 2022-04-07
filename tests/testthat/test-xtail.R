context("xtail")
test_that("xtail", {
    data(xtaildata)
    mrna <- xtaildata$mrna[1:100,]
    rpf <- xtaildata$rpf[1:100,]
    condition <- c("control","control","treat","treat")
    xtail <- xtail(mrna,rpf,condition,bins=100, threads = 2)
    expect_s4_class(xtail,"xtail")
    # accessors
    expect_s4_class(resultsTable(xtail),"DataFrame")
    expect_type(conditions(xtail),"character")
    # plots
    plot <- plotFCs(xtail)
    expect_s3_class(plot,"ggplot")
    plot <- plotRs(xtail)
    expect_s3_class(plot,"ggplot")
    plot <- volcanoPlot(xtail)
    expect_s3_class(plot,"ggplot")

    # SummarizedExperiment wrapper
    se <- SummarizedExperiment(assays = list(mrna = mrna, rpf = rpf),
                               colData = DataFrame(condition = condition))
    xtail2 <- runXtail(se, "mrna","rpf",condition = colData(se)$condition,
                       bins=100, threads = 2)
    expect_equal(xtail,xtail2)
    se <- addXtail(se, "mrna","rpf",condition = colData(se)$condition,
                   bins=100, threads = 2)
    expect_s4_class(se,"SummarizedExperiment")
})
