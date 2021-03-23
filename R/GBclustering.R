##' Guided Bayesian Clustering
##'
##' Guided Bayesian clustering integrating clinical dataset with gene expression dataset.
##' @title GBclustering
##' @param Y Gene expression matrix, n*p (rows for subjects and columns for genes).
##' @param U Adjusted R-squared or adjusted pseudo R-squared between phenotypic variable and expression value of each gene, a vector.
##' @param K Number of clusters.
##' @param n.iter Number of iteration.
##' @param hyperparameters A vector including the following hyperparameters in order: c, a_p, b_p, a_sigma, b_sigma, a_U0, b_U0, a_U1, b_U1, a_mu0, b_mu0, a_mu1, b_mu1.
##' @param showIteration Output progress or not.
##' @return a list consisting of
##' \item{L_PosterSamp}{Posterior samples of differential expression indicator L}
##' \item{Subtypes}{Clustering results}
##' \item{p}{Posterior mean of DE proportion p}
##' \item{tau_mu0}{Posterior mean of tau_mu0}
##' \item{tau_mu1}{Posterior mean of tau_mu1}
##' \item{tau_U0}{Posterior mean of tau_U0}
##' \item{tau_U1}{Posterior mean of tau_U1}
##' \item{pi}{Posterior mean of subtype proportions pi}
##' \item{mu}{Posterior mean of subtype mean mu}
##' \item{sigma_sq}{Posterior mean of subtype variance sigma_sq}
##' \item{tau_mu0_record}{Iteration record for tau_mu0}
##' \item{tau_mu1_record}{Iteration record for tau_mu1}
##' \item{tau_U0_record}{Iteration record for tau_U0}
##' \item{tau_U1_record}{Iteration record for tau_U1}
##' \item{mu_record}{Iteration record for mu}
##' \item{sigma_sq_record}{Iteration record for sigma_sq}
##' \item{BIC}{BIC value}
##' @export
##' @author Lingsong Meng


GBclustering <- function(Y, U, K, n.iter = 500, hyperparameters = c(1, 1, 1, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 2, 0.005, 4, 450), showIteration = T) {
    
    n.record <- floor(n.iter/2)
    n <- nrow(Y)
    G <- ncol(Y)
    
    # specify hyperparameters
    c <- hyperparameters[1]
    a_p <- hyperparameters[2]
    b_p <- hyperparameters[3]
    a_sigma <- hyperparameters[4]
    b_sigma <- hyperparameters[5]
    a_U0 <- hyperparameters[6]
    b_U0 <- hyperparameters[7]
    a_U1 <- hyperparameters[8]
    b_U1 <- hyperparameters[9]
    a_mu0 <- hyperparameters[10]
    b_mu0 <- hyperparameters[11]
    a_mu1 <- hyperparameters[12]
    b_mu1 <- hyperparameters[13]
    
    # set initial values
    DE_prop_t <- runif(1, 0, 1/2)
    tau_mu0_t <- 0.005
    tau_mu1_t <- 5
    tau_U0_t <- 0.005
    tau_U1_t <- 5
    pi_t <- NULL
    
    ## initial clusters
    ind.init.zero <- order(U, decreasing = T)[601:G]
    U.new <- U
    U.new[ind.init.zero] <- 0
    w <- abs(U.new/sum(U.new))
    Y.new <- sweep(Y[, w != 0], 2, sqrt(w[w != 0]), "*")
    Z_t <- kmeans(Y.new, K, nstart = 20)$cluster
    
    ## initial value of mu_t
    raw_Means <- matrix(nrow = K, ncol = G)
    for (k in seq_len(K)) {
        sumzt_k <- sum(Z_t == k)
        if (sumzt_k > 1) {
            raw_Means[k, ] <- colMeans(Y[Z_t == k, ])
        } else {
            raw_Means[k, ] <- Y[Z_t == k, ]
        }
    }
    mu_t <- t(raw_Means)
    
    ## initial value of L_t
    L_t <- update_L(mu_t, U, DE_prop_t, tau_mu0_t, tau_mu1_t, tau_U0_t, tau_U1_t)
    
    ## initial value of sigma_sq_t
    sigma_sq_t <- update_sigma_sq(mu_t, Z_t, a_sigma, b_sigma, Y)
    
    
    # record samples
    DE_prop_t_record <- array(NA, dim = c(n.record))
    tau_mu0_t_record <- array(NA, dim = c(n.record))
    tau_mu1_t_record <- array(NA, dim = c(n.record))
    tau_U0_t_record <- array(NA, dim = c(n.record))
    tau_U1_t_record <- array(NA, dim = c(n.record))
    L_t_record <- array(NA, dim = c(G, n.record))
    pi_t_record <- array(NA, dim = c(K, n.record))
    Z_t_record <- array(NA, dim = c(n, n.record))
    mu_t_record <- array(NA, dim = c(G, K, n.record))
    sigma_sq_t_record <- array(NA, dim = c(G, n.record))
    
    
    # Gibbs sampler
    t1 <- Sys.time()
    message("  running the Gibbs sampler ...\n")
    for (t in seq_len(n.iter)) {
        ## sequetially update each parameter
        DE_prop_t <- update_DE_prop(L_t, a_p, b_p)
        tau_mu_t <- update_tau_mu(L_t, mu_t, a_mu0, b_mu0, a_mu1, b_mu1)
        tau_mu0_t <- tau_mu_t[1]
        tau_mu1_t <- tau_mu_t[2]
        Guidance_t <- update_Guidance(U, L_t, a_U0, b_U0, a_U1, b_U1)
        tau_U0_t <- Guidance_t[1]
        tau_U1_t <- Guidance_t[2]
        L_t <- update_L(mu_t, U, DE_prop_t, tau_mu0_t, tau_mu1_t, tau_U0_t, tau_U1_t)
        pi_t <- update_pi(Z_t, c, K)
        Z_t <- update_Z_v2(Z_t, mu_t, sigma_sq_t, pi_t, Y)
        mu_t <- update_mu(K, L_t, Z_t, sigma_sq_t, tau_mu0_t, tau_mu1_t, Y)
        sigma_sq_t <- update_sigma_sq(mu_t, Z_t, a_sigma, b_sigma, Y)
        
        if (showIteration == TRUE) {
            message(c("  Iteration ", t, "\n"))
        }
        
        ## store sampler
        if (t > n.iter - n.record) {
            DE_prop_t_record[t - (n.iter - n.record)] <- DE_prop_t
            tau_mu0_t_record[t - (n.iter - n.record)] <- tau_mu0_t
            tau_mu1_t_record[t - (n.iter - n.record)] <- tau_mu1_t
            tau_U0_t_record[t - (n.iter - n.record)] <- tau_U0_t
            tau_U1_t_record[t - (n.iter - n.record)] <- tau_U1_t
            L_t_record[, t - (n.iter - n.record)] <- L_t
            pi_t_record[, t - (n.iter - n.record)] <- pi_t
            Z_t_record[, t - (n.iter - n.record)] <- Z_t
            mu_t_record[, , t - (n.iter - n.record)] <- mu_t
            sigma_sq_t_record[, t - (n.iter - n.record)] <- sigma_sq_t
        }
    }
    t2 <- Sys.time()
    message(paste0("  The Gibbs sampler takes: ", round(difftime(t2, t1, units = "mins"), 3), " mins", "\n"))
    
    
    # posterior
    output <- list()
    
    ## posterior samples of DE indicators L
    output[[1]] <- L_t_record
    
    ## posterior mode of Z
    Z_t_post <- vapply(seq_len(n), function(j) {
        temp <- table(Z_t_record[j, ])
        as.numeric(names(temp)[which.max(temp)])
    }, FUN.VALUE = 1)
    output[[2]] <- Z_t_post
    
    ## posterior mean of other parameters
    DE_prop_t_post <- mean(DE_prop_t_record)
    tau_mu0_t_post <- mean(tau_mu0_t_record)
    tau_mu1_t_post <- mean(tau_mu1_t_record)
    tau_U0_t_post <- mean(tau_U0_t_record)
    tau_U1_t_post <- mean(tau_U1_t_record)
    pi_t_post <- rowMeans(pi_t_record)
    mu_t_post <- matrix(nrow = G, ncol = K)
    for (g in seq_len(G)) {
        mu_t_post[g, ] <- rowMeans(mu_t_record[g, seq_len(K), ])
    }
    sigma_sq_t_post <- rowMeans(sigma_sq_t_record)
    
    output[[3]] <- DE_prop_t_post
    output[[4]] <- tau_mu0_t_post
    output[[5]] <- tau_mu1_t_post
    output[[6]] <- tau_U0_t_post
    output[[7]] <- tau_U1_t_post
    output[[8]] <- pi_t_post
    output[[9]] <- mu_t_post
    output[[10]] <- sigma_sq_t_post
    
    output[[11]] <- tau_mu0_t_record
    output[[12]] <- tau_mu1_t_record
    output[[13]] <- tau_U0_t_record
    output[[14]] <- tau_U1_t_record
    output[[15]] <- mu_t_record
    output[[16]] <- sigma_sq_t_record
    
    logL <- sum(log(pi_t_post[Z_t_post])) + sum(log(dnorm(t(Y), mu_t_post[, Z_t_post], sqrt(sigma_sq_t_post))))
    BIC <- (-2) * logL + K * G * log(n * G)
    output[[17]] <- BIC
    
    
    names(output) <- c("L_PosterSamp", "Subtypes", "p", "tau_mu0", "tau_mu1", "tau_U0", "tau_U1", "pi", "mu", "sigma_sq", "tau_mu0_record", "tau_mu1_record", 
        "tau_U0_record", "tau_U1_record", "mu_record", "sigma_sq_record", "BIC")
    
    return(output)
}

