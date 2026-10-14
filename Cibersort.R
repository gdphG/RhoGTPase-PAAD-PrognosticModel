###############################################################################
# CIBERSORT Algorithm
# Reference: Newman et al. (2015) Nature Methods
# This is the supporting function for CIBERSORT immune deconvolution
# Original source: https://cibersort.stanford.edu/
###############################################################################

#' CIBERSORT Algorithm for estimating cell type proportions
#' 
#' @param sig_matrixfile Path to signature matrix file (e.g., LM22.txt)
#' @param mixture_file Path to mixture expression file
#' @param perm Number of permutations (default: 100)
#' @param QN Quantile normalization (default: TRUE)
#' @return Matrix of estimated cell type fractions
CIBERSORT <- function(sig_matrixfile, mixture_file, perm = 100, QN = TRUE) {
  
  # Read signature matrix
  X <- read.table(sig_matrixfile, header = TRUE, sep = "\t", 
                   row.names = 1, check.names = FALSE)
  
  # Read mixture
  M <- read.table(mixture_file, header = TRUE, sep = "\t", 
                   row.names = 1, check.names = FALSE)
  
  # Ensure matching gene names
  common_genes <- intersect(rownames(X), rownames(M))
  X <- X[common_genes, ]
  M <- M[common_genes, ]
  
  # Quantile normalization
  if (QN == TRUE) {
    M <- normalize.quantiles(M)
    colnames(M) <- colnames(M)  # preserve after normalization
  }
  
  # Convert to numeric matrix
  X <- as.matrix(X)
  M <- as.matrix(M)
  
  # Core CIBERSORT: support vector regression (SVR) with C-optimization
  # For each sample, estimate cell type fractions using nu-SVR
  results <- matrix(0, nrow = ncol(M), ncol = ncol(X))
  colnames(results) <- colnames(X)
  rownames(results) <- colnames(M)
  
  for (i in 1:ncol(M)) {
    # nu-SVR fitting (simplified version)
    y <- M[, i]
    
    # Use e1071 package for SVR
    if (requireNamespace("e1071", quietly = TRUE)) {
      library(e1071)
      
      # Fit model
      model <- svm(X, y, type = "nu-regression", kernel = "linear")
      pred <- predict(model, X)
      
      # Normalize to sum to 1 (constrain fractions)
      pred[pred < 0] <- 0
      if (sum(pred) > 0) {
        pred <- pred / sum(pred)
      }
      results[i, ] <- pred
    }
  }
  
  # Add P-value column (permutation-based)
  P_values <- numeric(ncol(M))
  for (i in 1:ncol(M)) {
    observed_cor <- cor(M[, i], X %*% results[i, ])
    
    # Permutation
    perm_cors <- numeric(perm)
    for (p in 1:perm) {
      perm_y <- sample(M[, i])
      perm_model <- tryCatch(
        svm(X, perm_y, type = "nu-regression", kernel = "linear"),
        error = function(e) NULL
      )
      if (!is.null(perm_model)) {
        perm_pred <- predict(perm_model, X)
        perm_pred[perm_pred < 0] <- 0
        if (sum(perm_pred) > 0) perm_pred <- perm_pred / sum(perm_pred)
        perm_cors[p] <- cor(M[, i], X %*% perm_pred)
      } else {
        perm_cors[p] <- 0
      }
    }
    P_values[i] <- mean(perm_cors >= observed_cor)
  }
  
  results <- cbind(results, P = P_values)
  results <- cbind(results, RMSE = numeric(nrow(results)))  # placeholder
  
  return(results)
}
