# gaussian_kcockoff.py

## Generation of model-x knockoff following equi-correlated method or
## optimization scheme following Barber et al. (2015). Requires cvxopt.

is_posdef <- function(X, tol = 1e-14){
  # Check a matrix is positive definite by calculating eigenvalue of the matrix
  # Parameters
  # ----------
  # X : matrix (n_samples x n_features)
  #     matrix to check
  # tol : float, optional
  #     minimum threshold for eigenvalue
  #
  # Returns
  # -------
  # True or False 
  #
  eig_value <- eigen(X)$values
  return(all(eig_value > tol))
}


# The _s_equi function is used to estimate the diagonal matrix of correlation between real and knockoff variables using the equi-correlated equation. This estimation is a crucial step in the Model-X knockoff inference procedure.
# The equi-correlated equation is derived from the assumption that the correlation structure between the original features and their knockoff counterparts is approximately the same for all features. The goal is to estimate the diagonal elements of a matrix S, where S[i, i] represents the correlation between the i-th original feature and its knockoff counterpart.
# Here's a step-by-step explanation of the theory behind the _s_equi function:
# The input Sigma is the empirical covariance matrix calculated from the original design matrix.
# The empirical covariance matrix Sigma is transformed into a correlation matrix G using the cov2cor function. The correlation matrix G has the same pairwise correlations as Sigma, but with all variances normalized to 1.
# The eigenvalues of the correlation matrix G are computed using the eigen function. The eigenvalues represent the variances explained by each eigenvector.
# The minimum eigenvalue lambda_min is extracted from the eigenvalues. It provides a lower bound on the correlation between any original feature and its knockoff counterpart.
# The vector S is initialized with all elements set to min(2 * lambda_min, 1). This step ensures that the estimated correlations are non-negative and do not exceed the range of [-1, 1].
# The variable psd is set to FALSE to indicate that the matrix is not positive semi-definite (psd) initially. The variable s_eps is initialized to 0, representing a small increment factor for adjusting S.
# A loop is initiated to adjust the elements of S until the matrix 2 * G - diag(S * (1 - s_eps)) becomes positive semi-definite. This adjustment is necessary to satisfy the constraint that the estimated matrix diag(s) should be positive semi-definite.
# Inside the loop, the is.positive.definite function is used to check if the adjusted matrix is positive semi-definite. If it is not, the adjustment factor s_eps is incremented to perturb the elements of S further.
# Once the adjusted matrix is positive semi-definite, the loop exits, and the final estimation of S * diag(Sigma) is returned.
# In summary, the _s_equi function uses the equi-correlated equation to estimate the diagonal matrix of correlations between real and knockoff variables. The estimation ensures that the correlations satisfy the constraints of positive semi-definiteness and the range of [-1, 1]. This estimation is an important step in the Model-X knockoff inference procedure, which aims to control the false discovery rate in high-dimensional variable selection problems."

s_equi <- function(Sigma, max_iterations = 10000) {
  # Estimate diagonal matrix of correlation between real and knockoff variables using equi-correlated equation
  # Parameters
  #   Sigma: 2D ndarray (n_features, n_features)
  #     Empirical covariance matrix calculated from original design matrix
  # Returns
  #   1D ndarray (n_features, )
  #     Vector of diagonal values of estimated matrix diag{s}
  
  n_features <- nrow(Sigma)
  
  G <- cov2cor(Sigma) # Convert cov matrix to cor matrix
  eig_value <- eigen(G)$values
  lambda_min <- min(eig_value)
  S <- rep(min(2 * lambda_min, 1), n_features)
  
  psd <- FALSE
  s_eps <- 0
  iterations <- 0
  
  while (!psd && iterations < max_iterations) {
    # if all eigval > 0 then the matrix is psd
    psd <- is_posdef(2 * G - diag(S * (1 - s_eps)))
    if (!psd) {
      if (s_eps == 0) {
        s_eps <- 1e-08
      } else {
        s_eps <- s_eps * 10
      }
    }
    iterations <- iterations + 1
  }
  
  S <- S * (1 - s_eps)
  
  return(S * diag(Sigma))
}

### Alternative cov2cor
cov_to_corr <- function(Sigma) {
  features_std <- sqrt(diag(Sigma))
  Scale <- outer(features_std, features_std)
  Corr_matrix <- Sigma / Scale
  return(Corr_matrix)
}


library(glasso)
#The _estimate_distribution function is used to estimate the distribution of the
# data by estimating the mean and covariance matrix. 
# It provides two options for covariance matrix estimation: Ledoit-Wolf shrinkage
# and graphical LASSO.
# Ledoit-Wolf Shrinkage:
# Ledoit-Wolf shrinkage is a method for estimating the covariance matrix that 
# improves the accuracy of estimation, especially when the number of variables 
# is large compared to the number of samples. It incorporates shrinkage 
# intensity to balance between the sample covariance matrix and a target matrix 
# (usually a diagonal matrix with the sample variances). 
# This helps to reduce the bias in the estimation and improve the stability of 
# the covariance matrix. The cov.shrink function from the corpcor package is used to perform Ledoit-Wolf shrinkage.
#    Graphical LASSO:
#    Graphical LASSO is a method for estimating the precision matrix 
# (inverse of the covariance matrix) by imposing sparsity. 
# It assumes that the underlying graph structure of the variables is sparse,
# meaning that most of the entries in the precision matrix are zero. 
# The glasso function from the glasso package is used for graphical LASSO 
# estimation. The rho parameter in glasso controls the sparsity level of the 
# estimated precision matrix. Higher values of rho lead to sparser precision matrices.

# By estimating the mean and covariance matrix (or precision matrix), the _estimate_distribution function provides the necessary information for subsequent steps in the Model-X knockoff inference procedure.

estimate_distribution <- function(X, shrink=FALSE, cov_estimator="ledoit_wolf") {
  alphas <- c(1e-3, 1e-2, 1e-1, 1)
  
  mu <- colMeans(X)
  Sigma <- cov(X)
  
  if (shrink || !is_posdef(Sigma)) {
    if (cov_estimator == "ledoit_wolf") {
      Sigma_shrink <- corpcor::cov.shrink(Sigma)[,]
    } else if (cov_estimator == "graph_lasso") {
      model <- glasso(Sigma, rho = 0.5)
      Sigma_shrink <- solve(model$wi)
    } else {
      stop(paste(cov_estimator, "is not a valid covariance estimation method"))
    }
    
    return(list(mu=mu, Sigma_shrink=Sigma_shrink))
  }
  
  return(list(mu=mu, Sigma=Sigma))
}

gaussian_knockoff_generation <- function(X, mu, Sigma, method = "equi", seed = NULL) {
  n_samples <- nrow(X)
  n_features <- ncol(X)
  
  if (method == "equi") {
    Diag_s <- diag(s_equi(Sigma))
  } else {
    stop(paste(method, "is not a valid knockoff construction method"))
  }
  
  Sigma_inv_s <- solve(Sigma, Diag_s)
  
  # First part on the RHS of equation 1.4 in Barber & Candes (2015)
  Mu_tilde <- X - (X - mu) %*% Sigma_inv_s
  
  # To calculate the Cholesky decomposition later on
  Sigma_tilde <- 2 * Diag_s - Diag_s %*% Sigma_inv_s %*% Diag_s
  while (!is_posdef(Sigma_tilde)) {
    Sigma_tilde <- Sigma_tilde + 1e-10 * diag(n_features)
    warning("The conditional covariance matrix for knockoffs is not positive definite. Adding minor positive value to the matrix.")
  }
  
  set.seed(seed)
  U_tilde <- matrix(rnorm(n_samples * n_features), n_samples, n_features)
  # Equation 1.4 in Barber & Candes (2015)
  X_tilde <- Mu_tilde + U_tilde %*% chol(Sigma_tilde)
  
  return(X_tilde)
}

