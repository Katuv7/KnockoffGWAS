# Define the custom function to generate truncated random values
truncated_rnorm <- function(n, mean, sd, lower, upper) {
  samples <- numeric(n)  # Create an empty vector to store the generated values
  count <- 1  # Counter for the number of generated values
  
  while (count <= n) {
    value <- rnorm(1, mean, sd)  # Generate a random value from the normal distribution
    
    # Check if the value is outside the range (-1, 1)
    if (value < lower || value > upper) {
      samples[count] <- value  # Store the value in the vector
      count <- count + 1  # Increment the counter
    }
  }
  
  return(samples)  # Return the vector of generated values
}


simu_data <- function(n, p, rho = 0.25, snr = 2.0, sparsity = 0.06,
                      beta_mean = 0, beta_sd = 1, low_er = -1, up_per = 1, 
                      effect = 1.0, minMAF=0.05, maxMAF=0.5, seed = NULL) {
  # Function to simulate data following an autoregressive structure with Toeplitz covariance matrix
  # Parameters
  # ----------
  # n : integer
  #   number of observations
  # p : integer
  #   number of variables
  # sparsity : numeric, optional
  #   ratio of the number of variables with non-zero coefficients over total coefficients
  # rho : numeric, optional
  #   correlation parameter
  # beta_mean, beta_sd, low_er, up_per 
  #   parameters for the truncated normal distribution
  # effect : numeric
  #   effect for the coefficients
  # minMAF : numeric
  #   minimum value of the Minor Allele Frequency
  # maxMAF : numeric
  #   maximum value of the Minor Allele Frequency
  # seed : integer or NULL, optional
  #   random seed for generator
  
  # Setup seed generator
  set.seed(seed)
  
  # Number of non-null
  k <- as.integer(sparsity * p)
  
  # Generate the variables from a multivariate normal distribution
  mu <- rep(0, p)
  Sigma <- toeplitz(rho ^ (0:(p - 1)))  # covariance matrix of X
  Z <- rmvnorm(n, mean = mu, sigma = Sigma)

  # Generate minor allele frequencies for each SNP
  MAF <- runif(p, min = minMAF, max = maxMAF)  # It works asymptotically
  
  # Generate genotype vectors for each SNP
  X <- matrix(0, nrow = n, ncol = p)
  for (ii in seq_along(MAF)) {
    maf <- MAF[ii]
    c1 <- qnorm(maf/2) 
    c2 <- qnorm(2*maf)
    vect <- ifelse(Z[,ii] > c2, 0, ifelse(Z[,ii] < c1, 2, 1))
    X[, ii] <- vect
  }
  
  # Generate the response from a linear model
  non_zero <- sample(1:p, k)
  beta_true <- rep(0, p)
  causal_true <- rep(0, p)
  causal_true[non_zero] <- effect
  beta_true[non_zero] <- truncated_rnorm(n = k, mean = beta_mean, sd = beta_sd,
                                         lower = low_er, upper = up_per)
  eps <- rnorm(n)
  prod_temp <- X %*% beta_true
  noise_mag <- sqrt(sum(prod_temp ^ 2) / (snr * sum(eps ^ 2)))
  y <- prod_temp + noise_mag * eps
  
  return(list(genotypes = X, 
              phenotypes = y, 
              causal_true = causal_true, 
              non_zero = non_zero))
}
