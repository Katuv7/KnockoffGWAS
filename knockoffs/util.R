quantile_aggregation <- function(pvals, gamma = 0.5, gamma_min = 0.05, adaptive = FALSE) {
  if (adaptive) {
    return(adaptive_quantile_aggregation(pvals, gamma_min))
  } else {
    return(fixed_quantile_aggregation(pvals, gamma))
  }
}


fixed_quantile_aggregation <- function(pvals, gamma = 0.3) {
  # Quantile aggregation function based on Meinshausen et al (2008)
  # Parameters
  # ----------
  # pvals : 2D ndarray (n_bootstrap, n_test)
  #      p-value (adjusted)
  #  gamma : float
  #      Percentile value used for aggregation.
  #  Returns
  #  -------
  #  1D ndarray (n_tests, )
  #      Vector of aggregated p-value
  converted_score <- (1 / gamma) * quantile(pvals, probs = gamma, na.rm = TRUE)
  return(pmin(1, converted_score))
}


adaptive_quantile_aggregation <- function(pvals, gamma_min = 0.05) {
  # adaptive version of the quantile aggregation method, Meinshausen et al. (2008)
  gammas <- seq(gamma_min, 1, 0.05)
  list_Q <- sapply(gammas, function(gamma) fixed_quantile_aggregation(pvals, gamma))
  return(pmin(1, (1 - log(gamma_min)) * min(list_Q)))
}



fdr_threshold <- function(pvals, fdr = 0.1, method = "bhq", reshaping_function = NULL) {
  if (method == "bhq") {
    return(bhq_threshold(pvals, fdr))
  } else if (method == "bhy") {
    return(bhy_threshold(pvals, fdr, reshaping_function))
  } else {
    stop(paste(method, "is not a supported FDR control method"))
  }
}

cal_fdp_power <- function(selected, non_zero_index, r_index = TRUE) {
  if (length(selected) == 0) {
    return(list(fdp = 0.0, power = 0.0))
  }
  
  if (r_index) {
    selected <- selected - 1
  }
  
  true_positive <- selected[selected %in% non_zero_index]
  false_positive <- selected[!(selected %in% non_zero_index)]
  fdp <- length(false_positive) / max(1, length(selected))
  power <- length(true_positive) / length(non_zero_index)
  
  return(list(fdp = fdp, power = power))
}

bhq_threshold <- function(pvals, fdr = 0.1) {
  n_features <- length(pvals)
  pvals_sorted <- sort(pvals)
  selected_index <- 2 * n_features
  for (i in (n_features:1)) {
    if (pvals_sorted[i] <= fdr * (i) / n_features) {
      print(i)
      selected_index <- i
      break
    }
  }
  if (selected_index <= n_features) {
    return(pvals_sorted[selected_index])
  } else {
    return(-1.0)
  }
}

bhy_threshold <- function(pvals, reshaping_function = NULL, fdr = 0.1) {
  # Benjamini-Hochberg-Yekutieli procedure for controlling FDR, with input
  # shape function. Reference: Ramdas et al (2017)
  n_features <- length(pvals)
  pvals_sorted <- sort(pvals)
  selected_index <- 2 * n_features
  
  if (is.null(reshaping_function)) {
    temp <- seq_len(n_features)
    sum_inverse <- sum(1 / (temp + 1))
    return(bhq_threshold(pvals, fdr / sum_inverse))
  } else {
    for (i in (n_features:1)) {
      if (pvals_sorted[i] <= fdr * reshaping_function(i + 1) / n_features) {
        selected_index <- i
        break
      }
    }
    if (selected_index <= n_features) {
      return(pvals_sorted[selected_index])
    } else {
      return(-1.0)
    }
  }
}

############################################
# FDP
############################################
fdp <- function(cont_table) {
  
  # Extract values from contingency table
  TN <- cont_table[1,1]
  FN <- cont_table[1,2]
  FP <- cont_table[2,1]
  TP <- cont_table[2,2]
  
  # Calculate FDP
  if ((FP + TP) == 0) {
    fdp <- 0
  } else {
    fdp <- FP / (FP + TP)
  }
  
  return(fdp)
}
###############################################################
# present.results
# Function to compute ROC curve, with AUC and FDP from threshold
###############################################################

present.results<-function(pvalues,causal.true,threshold=0.4){
  pred <- prediction(1-pvalues, causal.true)
  par(mfrow=c(1,2))
  # ROC curve
  perf <- performance(pred, "tpr", "fpr")
  plot(perf,colorize=TRUE, lwd= 3,main= "ROC curve")
  # Precision Recall curve
  perf <- performance(pred, "prec", "rec")
  plot(perf, colorize=TRUE, lwd= 3,main= "Precision/Recall")
  
  # Confusion table for given theshold
  causal.estimated=factor(p.adjust(pvalues,"BH")< threshold, levels=c(FALSE,TRUE))
  print(test.table<-table(causal.estimated,causal.true))
  
  print(paste("FDP",fdp(test.table)),digits=3)
  # AUC
  print(paste("AUC: ",round(performance(pred,"auc")@y.values[[1]],digits=3)))
}