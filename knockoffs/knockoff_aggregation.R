# Algorithm 1 AKO Aggregation of multiple knockoffs
library(Matrix)
library(doParallel)
library(glmnet)
library(knockoff)

# Definition 3 (Intermediate p-value)

empirical_pval <- function(test_score, offset = 1) {
  # Function of pi_j (Eq (5) intermediate p-value)
  
  pvals <- vector()
  n_features <- length(test_score)
  
  if (offset != 0 && offset != 1) {
    stop("'offset' must be either 0 or 1")
  }
  
  test_score_inv <- -test_score
  for (i in 1:n_features) {
    if (test_score[i] <= 0) {
      pvals <- c(pvals, 1)
    } else {
      pvals <- c(pvals,
                 (offset + sum(test_score_inv >= test_score[i])) /
                   n_features
      )
    }
  }
  
  return(pvals)
}
