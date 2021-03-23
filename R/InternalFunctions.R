# Internal functions in GBclustering function

# A function for sampling from the Dirichlet distribution
Dirichlet <- function(alpha_vec) {
    # input is a parameter vector output is one sample from Dirichlet distribution
    
    n <- length(alpha_vec)
    gamma_rvs <- vapply(seq_len(n), function(i) {
        rgamma(1, shape = alpha_vec[i], rate = 1)
    }, FUN.VALUE = 1)
    return(gamma_rvs/sum(gamma_rvs))
}

# A function for sampling from the Inverse Gamma distribution
Inv_Gamma <- function(a, b) {
    # input is shape a and scale b of Inv-Gamma distribution output is one sample from the inverse gamma distribution
    tmp <- rgamma(1, shape = a, rate = b)
    return(1/tmp)
}



# 1. Update the differentially expressed (DE) indicator proportion
update_DE_prop <- function(L_t, a_p, b_p) {
    # a_p and b_p are the hyperparameters for the DE proportion parameter
    G <- length(L_t)
    s <- sum(L_t)
    tmp <- rbeta(1, shape1 = s + a_p, shape2 = G - s + b_p)
    return(tmp)
}



# 2. Update the standard deviation for the spike and slab part
update_tau_mu <- function(L_t, mu_t, a_mu0, b_mu0, a_mu1, b_mu1) {
    K <- ncol(mu_t)
    ind <- which(L_t == 0)
    s1 <- length(ind)
    s2 <- sum(mu_t[ind, ]^2)
    tau_mu0 <- sqrt(Inv_Gamma(a_mu0 + K/2 * s1, b_mu0 + 1/2 * s2))
    ind <- which(L_t == 1)
    s1 <- length(ind)
    s2 <- sum(mu_t[ind, ]^2)
    tau_mu1 <- sqrt(Inv_Gamma(a_mu1 + K/2 * s1, b_mu1 + 1/2 * s2))
    return(c(tau_mu0, tau_mu1))
}



# 3. Update the sigma for Gaussian mixture model
update_Guidance <- function(U, L_t, a_U0, b_U0, a_U1, b_U1) {
    U1 <- U[which(L_t == 1)]  ## Guidance for intrinsic genes
    U0 <- U[which(L_t == 0)]  ## Guidance for non-intrinsic genes
    tau_U0 <- sqrt(Inv_Gamma(a_U0 + length(U0)/2, b_U0 + sum((U0 - 0)^2)/2))
    tau_U1 <- sqrt(Inv_Gamma(a_U1 + length(U1)/2, b_U1 + sum((U1 - 1)^2)/2))
    return(c(tau_U0, tau_U1))
}



# 4. (fast) Update the DE indicators
update_L <- function(mu_t, U, DE_prop_t, tau_mu0_t, tau_mu1_t, tau_U0_t, tau_U1_t) {
    # update L_g, L_g = 1 if the expression level of gene g is differentially expressed output is a vector with G components
    G <- nrow(mu_t)
    prob1 <- DE_prop_t * apply(dnorm(mu_t, 0, tau_mu1_t), 1, prod) * dnorm(U, 1, tau_U1_t)
    prob2 <- (1 - DE_prop_t) * apply(dnorm(mu_t, 0, tau_mu0_t), 1, prod) * dnorm(U, 0, tau_U0_t)
    r <- runif(G)
    L <- as.numeric(r < prob1/(prob1 + prob2))
    return(L)
}



# 5. Update the subtype proportions
update_pi <- function(Z_t, c, K) {
    # update the subtype proportion Z_t is a vector containing n-component where the entry Z_t[j] stands for the subtype indicator of jth subject c is the
    # hyperparameter vector in the Dirichlet prior output is a vector with K components
    prob <- vapply(seq_len(K), function(k) sum(Z_t == k), FUN.VALUE = 1)
    tmp <- Dirichlet(prob + c)
    return(tmp)
}



# 6. (fast) Update the subtype indicator for samples
update_Z_v2 <- function(Z_t, mu_t, sigma_sq_t, pi_t, Y) {
    # update Z by Metropolis-Hasting step Z_t is a n-component vector where the entry Z_t[j] stands for the subtype indicator of jth subject mu_t is the G
    # by K subtype effect matrix sigma_sq_t is the variance vector with G components pi_t is a vector with K components Y is the gene expression data
    # matrix output is a vector containing n-component
    K <- ncol(mu_t)
    n <- nrow(Y)
    Z_tmp <- Z_t
    # get proposal
    proposal <- sample(seq_len(K), n, replace = T)
    # calculate posterior ratio in log scale
    aa <- colSums(-(t(Y) - mu_t[, proposal])^2/(2 * sigma_sq_t))
    bb <- colSums(-(t(Y) - mu_t[, Z_tmp])^2/(2 * sigma_sq_t))
    # calculate the acceptance prob
    prob <- exp(aa - bb) * pi_t[proposal]/pi_t[Z_tmp]
    prob[prob > 1] <- 1
    tmp <- runif(n)
    ind <- tmp <= prob
    Z_tmp[ind] <- proposal[ind]
    return(Z_tmp)
}



# 7. (fast) Update the subtype effects
update_mu <- function(K, L_t, Z_t, sigma_sq_t, tau_mu0_t, tau_mu1_t, Y) {
    # update subtype effect L_t is a vector with G components Z_t is a vector containing n components Z_t[j] stands for the subtype indicator of jth
    # subject sigma_sq_t is the variance vector with G components output is a G by K matrix
    n <- nrow(Y)
    G <- ncol(Y)
    g0 <- which(L_t == 0)
    g1 <- which(L_t == 1)
    mu <- matrix(nrow = G, ncol = K)
    for (k in seq_len(K)) {
        mean.sd <- matrix(nrow = G, ncol = 2)
        j <- which(Z_t == k)
        # intrinsic genes
        s1 <- tau_mu1_t^2 * rowSums(t(Y[j, g1])/sigma_sq_t[g1])
        s2 <- tau_mu1_t^2 * length(j)/sigma_sq_t[g1] + 1
        mean.sd[g1, 1] <- s1/s2
        mean.sd[g1, 2] <- sqrt(tau_mu1_t^2/s2)
        # non-intrinsic genes
        s1 <- tau_mu0_t^2 * rowSums(t(Y[j, g0])/sigma_sq_t[g0])
        s2 <- tau_mu0_t^2 * length(j)/sigma_sq_t[g0] + 1
        mean.sd[g0, 1] <- s1/s2
        mean.sd[g0, 2] <- sqrt(tau_mu0_t^2/s2)
        mu[, k] <- apply(mean.sd, 1, function(x) rnorm(1, x[1], x[2]))
    }
    return(mu)
}



# 8. (fast) Update the variances
update_sigma_sq <- function(mu_t, Z_t, a_sigma, b_sigma, Y) {
    # update variance mu_t is the G by K subtype effect matrix Z_t is a vector containing n-component where Z_t[j] stands for the subtype indicator of jth
    # subject output is a vector with G components
    n <- nrow(Y)
    shape <- a_sigma + n/2
    rate <- b_sigma + rowSums((t(Y) - mu_t[, Z_t])^2)/2
    sigma_sq <- vapply(rate, function(x) Inv_Gamma(shape, x), FUN.VALUE = 1)
    return(sigma_sq)
}

