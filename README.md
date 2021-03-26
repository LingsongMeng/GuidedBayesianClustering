# GuidedBayesianClustering
Github repository for Outcome-guided Bayesian Clustering Method (GuidedBayesianClustering)


## Install This Package from Github
In R console

```{R}
library(devtools)
install_github("LingsongMeng/GuidedBayesianClustering") 
```


## Short Tutorial
This short tutorial is about how to apply GuidedBayesianClustering to obtain subtype clustering result and gene selection result.

```{R}
library(GuidedSparseKmeans)
n <- 100 ## number of subjects
g <- 5000 ## number of gene features
set.seed(123)
epxression <- matrix(rnorm(n*g), nrow=g, ncol=n) ## briefly simulate gene expression data
outcome <- rnorm(n) ## simulate a continuous outcome variable


# First, obtain outcome guidance
U <- getU(t(epxression), outcome, model="linear")

# Second, Gibbs sampling and obtain clustering result
output <- GBclustering(Y=t(epxression), U=U, K=3, n.iter = 500,  
                       hyperparameters = c(1,1,1,0.001,0.001,0.001,0.001,0.001,0.001,2,0.005,4,450),
                       showIteration = T)
clusters <- output$Subtypes

# Then, obtain gene selection result 
belief <- rowSums(output$L_PosterSamp)/250
qval <- BayesianFDR(belief)
ind.gene <- which(qval < 0.05)

```
