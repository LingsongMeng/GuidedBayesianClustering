##' Get Bayesian FDR
##'
##' Obtain Bayesian FDR after using GBclustering.
##' @title BayesianFDR
##' @param belief Estilated probability of each gene to be intrinsic.
##' @return A vector
##' \item{qvalue}{A vector of qvalues}
##' @export
##' @author Lingsong Meng


BayesianFDR <- function(belief) {
    pvalue <- 1 - belief
    pvalue_order <- order(pvalue)
    sortedP <- pvalue[pvalue_order]
    sortedQ <- cumsum(sortedP)/(1:length(sortedP))
    qvalue <- sortedQ[match(1:length(pvalue_order), pvalue_order)]
    return(qvalue)
}

