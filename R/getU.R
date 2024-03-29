##' Get adjusted R-squared or adjusted pseudo R-squared
##'
##' Obtain adjusted R-squared or adjusted pseudo R-squared between phenotypic variable and expression value of each gene,
##' which will be used in GBclustering function.
##' @title getU
##' @param x Gene expression matrix, n*p (rows for subjects and columns for genes).
##' @param z One phenotypic variable from clinical dataset, a vector.
##' @param model The model fitted to obtain R2, please select model from 'linear', 'logit', 'exp', 'polr','cox'.
##' @return A vector
##' \item{U}{adjusted R-squared or adjusted pseudo R-squared between phenotypic variable and expression value of each gene, a vector.}
##' @export
##' @author Lingsong Meng


getU <- function(x, z, model) {
    
    if (is.null(model) || model %in% c("linear", "logit", "exp", "polr", "cox") != TRUE) 
        stop("Must select one from 'linear', 'logit', 'exp', 'polr','cox'.")
    
    if (model == "linear") {
        R2.per <- as.vector(cor(x, z))^2
    } else if (model == "logit") {
        R2.per <- apply(x, 2, function(g) PseudoR2(glm(z ~ g, family = binomial(link = "logit")), which = "CoxSnell"))
    } else if (model == "exp") {
        R2.per <- apply(x, 2, function(g) PseudoR2(glm(z ~ g, family = poisson), which = "CoxSnell"))
    } else if (model == "polr") {
        R2.per <- apply(x, 2, function(g) PseudoR2(polr(z ~ g, method = "logistic"), which = "CoxSnell"))
    } else {
        cox <- apply(x, 2, function(g) coxph(z ~ g, method = "breslow"))
        R2.per <- unlist(lapply(cox, function(g) 1 - exp(2 * (g$loglik[1] - g$loglik[2])/g$n * 2)))
    }
    
    U <- (R2.per - min(R2.per))/(max(R2.per) - min(R2.per))
    
    return(U)
}
