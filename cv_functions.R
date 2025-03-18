distance_cv <- function(Y, K, n_sams, n_burn, n_skip) {
     # settings
     require(doParallel)
     require(ROCR)

     # Cross-validation settings
     L <- 5
     folds <- get_folds(M = nrow(Y), J = 1, L = L)
     
     # Set up parallel computing
     cl <- makeCluster(min(detectCores(), L))
     registerDoParallel(cl)
     
     # Exportar variables necesarias al cluster
     clusterExport(cl, c("Y", "K", "n_sams", "n_burn", "n_skip", "folds"), envir = environment())
     
     # Cross-validation using parallel processing
     cv <- foreach(l = seq_len(L), .inorder = FALSE, .packages = c("Rcpp")) %dopar% {
          # Load necessary functions
          source   ("r_functions.R")
          source   ("dist_functions.R")
          sourceCpp("dist_functions.cpp")
          
          # Identify indices for the current fold
          na_indices <- folds == l
          
          # Extract true values and create masked dataset
          y_true <- Y[na_indices]
          Yna <- Y
          Yna[na_indices] <- NA
          
          # Run model with missing values
          y_ppp <- YPPP(Yna, na_indices, K, n_sams, n_burn, n_skip)
          
          # Return results as a named list
          list(l = l, K = K, y_true = y_true, y_ppp = y_ppp)
     }
     
     # Stop parallel processing
     stopCluster(cl)
     
     # Initialize AUCs vector
     aucs <- numeric(L)
     
     # Compute ROC curves and AUCs
     for (l in seq_len(L)) {
          # Compute performance metrics
          pred <- prediction(predictions = cv[[l]]$y_ppp, labels = cv[[l]]$y_true)
          perf_roc <- performance(pred, measure = "tpr", x.measure = "fpr")
          
          # Compute AUC
          perf_auc <- performance(pred, measure = "auc")
          aucs[l] <- perf_auc@y.values[[1]]
     }
     
     return(aucs)
}

class_cv <- function(Y, K, n_sams, n_burn, n_skip) {
     # settings
     require(doParallel)
     require(ROCR)
     
     # Cross-validation settings
     L <- 5
     folds <- get_folds(M = nrow(Y), J = 1, L = L)
     
     # Set up parallel computing
     cl <- makeCluster(min(detectCores(), L))
     registerDoParallel(cl)
     
     # Exportar variables necesarias al cluster
     clusterExport(cl, c("Y", "K", "n_sams", "n_burn", "n_skip", "folds"), envir = environment())
     
     # Cross-validation using parallel processing
     cv <- foreach(l = seq_len(L), .inorder = FALSE, .packages = c("Rcpp")) %dopar% {
          # Load necessary functions
          source   ("r_functions.R")
          source   ("class_functions.R")
          sourceCpp("class_functions.cpp")
          
          # Identify indices for the current fold
          na_indices <- folds == l
          
          # Extract true values and create masked dataset
          y_true <- Y[na_indices]
          Yna <- Y
          Yna[na_indices] <- NA
          
          # Run model with missing values
          y_ppp <- YPPP(Yna, na_indices, K, n_sams, n_burn, n_skip)
          
          # Return results as a named list
          list(l = l, K = K, y_true = y_true, y_ppp = y_ppp)
     }
     
     # Stop parallel processing
     stopCluster(cl)
     
     # Initialize AUCs vector
     aucs <- numeric(L)
     
     # Compute ROC curves and AUCs
     for (l in seq_len(L)) {
          # Compute performance metrics
          pred <- prediction(predictions = cv[[l]]$y_ppp, labels = cv[[l]]$y_true)
          perf_roc <- performance(pred, measure = "tpr", x.measure = "fpr")
          
          # Compute AUC
          perf_auc <- performance(pred, measure = "auc")
          aucs[l] <- perf_auc@y.values[[1]]
     }
     
     return(aucs)
}

eigen_cv <- function(Y, K, n_sams, n_burn, n_skip) {
     # settings
     require(doParallel)
     require(ROCR)
     
     # Cross-validation settings
     L <- 5
     folds <- get_folds(M = nrow(Y), J = 1, L = L)
     
     # Set up parallel computing
     cl <- makeCluster(min(detectCores(), L))
     registerDoParallel(cl)
     
     # Exportar variables necesarias al cluster
     clusterExport(cl, c("Y", "K", "n_sams", "n_burn", "n_skip", "folds"), envir = environment())
     
     # Cross-validation using parallel processing
     cv <- foreach(l = seq_len(L), .inorder = FALSE, .packages = c("Rcpp")) %dopar% {
          # Load necessary functions
          source   ("r_functions.R")
          source   ("eigen_functions.R")
          sourceCpp("eigen_functions.cpp")
          
          # Identify indices for the current fold
          na_indices <- folds == l
          
          # Extract true values and create masked dataset
          y_true <- Y[na_indices]
          Yna <- Y
          Yna[na_indices] <- NA
          
          # Run model with missing values
          y_ppp <- YPPP(Yna, na_indices, K, n_sams, n_burn, n_skip)
          
          # Return results as a named list
          list(l = l, K = K, y_true = y_true, y_ppp = y_ppp)
     }
     
     # Stop parallel processing
     stopCluster(cl)
     
     # Initialize AUCs vector
     aucs <- numeric(L)
     
     # Compute ROC curves and AUCs
     for (l in seq_len(L)) {
          # Compute performance metrics
          pred <- prediction(predictions = cv[[l]]$y_ppp, labels = cv[[l]]$y_true)
          perf_roc <- performance(pred, measure = "tpr", x.measure = "fpr")
          
          # Compute AUC
          perf_auc <- performance(pred, measure = "auc")
          aucs[l] <- perf_auc@y.values[[1]]
     }
     
     return(aucs)
}

sociality_cv <- function(Y, n_iter, n_burn, n_thin) {
  # settings
  require(doParallel)
  require(ROCR)
  
  # hyperparameters
  a_sigma <- 2 
  b_sigma <- 1/3
  a_tau   <- 2 
  b_tau   <- 1/3
  
  I <- get_I(Y)
  
  # Cross-validation settings
  L <- 5
  folds <- get_folds(M = nrow(Y), J = 1, L = L)
  
  # Set up parallel computing
  cl <- makeCluster(min(detectCores(), L))
  registerDoParallel(cl)

  # Exportar variables necesarias al cluster
  clusterExport(cl, c("Y", "I", "n_iter", "n_burn", "n_thin", "folds", "a_sigma", "b_sigma", "a_tau", "b_tau"), envir = environment())

  # Cross-validation using parallel processing
  cv <- foreach(l = seq_len(L), .inorder = FALSE) %dopar% {
       # Load necessary functions
       source("MCMC.R")
       source("r_functions.R")

       # Identify indices for the current fold
       na_indices <- folds == l

       # Extract true values and create masked dataset
       y_true <- Y[na_indices]
       Yna <- Y
       Yna[na_indices] <- NA

       # convert Yna and na_indices to matrices for compability
       Yna_matrix <- matrix(0, I, I)
       na_indices_matrix <- matrix(FALSE, I, I)
       for (i in 1:(I-1)) {
            for (ii in (i+1):I) {
                 k <- get_k(i, ii, I)
                 na_indices_matrix [i, ii] <- na_indices[k]
                 na_indices_matrix [ii, i] <- na_indices[k]
                 Yna_matrix[i, ii] <- Yna[k]
                 Yna_matrix[ii, i] <- Yna[k]
            }
       }

       # Run model with missing values
       y_ppp <- YPPP(Yna_matrix, na_indices_matrix, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

       # Return results as a named list
       list(l = l, y_true = y_true, y_ppp = y_ppp)
  }
  
  # Stop parallel processing
  stopCluster(cl)
  
  # Initialize AUCs vector
  aucs <- numeric(L)
  
  # Compute ROC curves and AUCs
  for (l in seq_len(L)) {
    # Compute performance metrics
    pred <- prediction(predictions = cv[[l]]$y_ppp, labels = cv[[l]]$y_true)
    perf_roc <- performance(pred, measure = "tpr", x.measure = "fpr")
    
    # Compute AUC
    perf_auc <- performance(pred, measure = "auc")
    aucs[l] <- perf_auc@y.values[[1]]
  }
  
  return(aucs)
}