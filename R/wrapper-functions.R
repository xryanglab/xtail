# Fit beta coefficients for negative binomial GLM
#
# This function estimates the coefficients (betas) for negative binomial generalized linear models.
#
# ySEXP n by m matrix of counts
# xSEXP m by k design matrix
# nfSEXP n by m matrix of normalization factors
# alpha_hatSEXP n length vector of the disperion estimates
# contrastSEXP a k length vector for a possible contrast
# beta_matSEXP n by k matrix of the initial estimates for the betas
# lambdaSEXP k length vector of the ridge values
# tolSEXP tolerance for convergence in estimates
# maxitSEXP maximum number of iterations
# useQRSEXP whether to use QR decomposition
#
# Note: at this level the betas are on the natural log scale
.fit_Beta <- function(ySEXP,
                      xSEXP,
                      nfSEXP,
                      alpha_hatSEXP,
                      contrastSEXP,
                      beta_matSEXP,
                      lambdaSEXP,
                      tolSEXP,
                      maxitSEXP,
                      useQRSEXP) {
    if ( missing(contrastSEXP) ) {
        # contrast is not required, just give 1,0,0,...
        contrastSEXP <- c(1,rep(0,ncol(xSEXP)-1))
    }
    # test for any NAs in arguments
    arg.names <- names(formals(.fit_Beta))
    na.test <- vapply(mget(arg.names), function(x) any(is.na(x)),
                      logical(1))
    if (any(na.test)) {
        stop(paste("in call to fitBeta, the following arguments contain NA:",
                   paste(arg.names[na.test],collapse=", ")),
             call. = FALSE)
    }
    fitBeta2(ySEXP=ySEXP, xSEXP=xSEXP, nfSEXP=nfSEXP, alpha_hatSEXP=alpha_hatSEXP,
             contrastSEXP=contrastSEXP, beta_matSEXP=beta_matSEXP,
             lambdaSEXP=lambdaSEXP, tolSEXP=tolSEXP, maxitSEXP=maxitSEXP,
             useQRSEXP=useQRSEXP)
}

.xtail_test_wrapper <- function(X){
    res <- xtail_test(counts1[X,],
                      counts2[X,],
                      intercept1[X],
                      intercept2[X],
                      log2Ratio1[X],
                      log2Ratio2[X],
                      dispersion1[X,],
                      dispersion2[X,],
                      sizefactor1,
                      sizefactor2,
                      cond1,
                      cond2,
                      bins,
                      ci)
    res
}
