# stat_coed_diff.py

stat_coef_diff <- function (X, X_k, y, family = "gaussian", cores = 2, ...) 
{
  if (!requireNamespace("glmnet", quietly = T)) 
    stop("glmnet is not installed", call. = F)
  parallel = T
  if (!requireNamespace("doParallel", quietly = T)) {
    warning("doParallel is not installed. Without parallelization, the statistics will be slower to compute", 
            call. = F, immediate. = T)
    parallel = F
  }
  if (!requireNamespace("parallel", quietly = T)) {
    warning("parallel is not installed. Without parallelization, the statistics will be slower to compute.", 
            call. = F, immediate. = T)
    parallel = F
  }
  if (parallel) {
    ncores = parallel::detectCores(all.tests = TRUE, logical = TRUE)
    if (cores == 2) {
      cores = min(2, ncores)
    }
    else {
      if (cores > ncores) {
        warning(paste("The requested number of cores is not available. Using instead", 
                      ncores, "cores"), immediate. = T)
        cores = ncores
      }
    }
    if (cores > 1) {
      doParallel::registerDoParallel(cores = cores)
      parallel = TRUE
    }
    else {
      parallel = FALSE
    }
  }
  swap = rbinom(ncol(X), 1, 0.5)
  swap.M = matrix(swap, nrow = nrow(X), ncol = length(swap), 
                  byrow = TRUE)
  X.swap = X * (1 - swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1 - swap.M)
  p = ncol(X)
  glmnet.coefs = cv_coeffs_glmnet(cbind(X.swap, Xk.swap), 
                                  y, family = family, parallel = parallel, ...)
  if (family == "multinomial") {
    Z <- abs(glmnet.coefs[[1]][2:(2 * p + 1)])
    for (b in 2:length(glmnet.coefs)) {
      Z <- Z + abs(glmnet.coefs[[b]][2:(2 * p + 1)])
    }
  }
  else if (family == "cox") {
    Z <- glmnet.coefs[1:(2 * p)]
  }
  else {
    Z <- glmnet.coefs[2:(2 * p + 1)]
  }
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig + p])
  W = W * (1 - 2 * swap)
  if (parallel) {
    if (cores > 1) {
      doParallel::stopImplicitCluster()
    }
  }
  return(W)
}

cv_coeffs_glmnet <- function (X, y, nlambda = 500, intercept = T, parallel = T, ...) 
{
  X = scale(X)
  n = nrow(X)
  p = ncol(X)
  if (!methods::hasArg(family)) 
    family = "gaussian"
  else family = list(...)$family
  if (!methods::hasArg(lambda)) {
    if (identical(family, "gaussian")) {
      if (!is.numeric(y)) {
        stop("Input y must be numeric.")
      }
      lambda_max = max(abs(t(X) %*% y))/n
      lambda_min = lambda_max/2000
      k = (0:(nlambda - 1))/nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }
  cv.glmnet.fit <- glmnet::cv.glmnet(X, y, lambda = lambda, 
                                     intercept = intercept, standardize = F, standardize.response = F, 
                                     parallel = parallel, ...)
  return(coef(cv.glmnet.fit, s = "lambda.min"))
}

coef_diff_threshold <- function (W, fdr = 0.1, offset = 1) 
{
  # Calculate the knockoff threshold based on the procedure stated in the article
  
  # Parameters 
  # ----------
  # test_score : numeric vector
  #   vector of test statistic
  # fdr : numeric, optional
  #   desired controlled FDR level
  # offset : integer, 0 or 1, optional
  #   offset equals 1 is the knockoff+ procedure
  
  if (offset != 1 && offset != 0) {
    stop("Input offset must be either 0 or 1")
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t) (offset + sum(W <= -t))/max(1, 
                                                             sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}
