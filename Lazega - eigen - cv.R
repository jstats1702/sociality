# Clear the global environment
rm(list = ls(all.names = TRUE))

# Set working directory
setwd("C:/Users/User/Dropbox/PAPERS/projects/sociality")

# Load required libraries and source files
library(igraph)
library(sand)
library(doParallel)
library(ROCR)

source("r_functions.R")

# Load data
g <- graph_from_data_frame(d = elist.lazega, directed = "F")
Y <- as.matrix(as_adjacency_matrix(graph = g, names = F))
I <- dim(Y)[1]

zro <- (1:I)[ Y %*% rep(1, I) == 0 ]
if (length(zro) > 0) {
     Y <- Y[-zro, -zro]
     Ylabs <- Ylabs[-zro]
}

Y[is.na(Y)] <- 0
Y[Y != 0 ]  <- 1

I <- dim(Y)[1]
Ycube <- Y

Y <- as.matrix(rep(NA, I*(I-1)/2))
for (i in 1:(I-1))
     for (ii in (i+1):I)
          Y[get_k(i, ii, I)] <- Ycube[i, ii]

# CV

# Cross-validation settings
L <- 5
set.seed(42)
folds <- get_folds(M = nrow(Y), J = 1, L = L)

# Model fitting settings
K      <- 4
n_sams <- 25000
n_burn <- 10000
n_skip <- 10

# # Set up parallel computing
# cl <- makeCluster(min(detectCores(), L))
# registerDoParallel(cl)
# 
# # Exportar variables necesarias al cluster
# clusterExport(cl, c("Y", "K", "n_sams", "n_burn", "n_skip", "folds"), envir = environment())
# 
# # Cross-validation using parallel processing
# cv <- foreach(l = seq_len(L), .inorder = FALSE, .packages = c("Rcpp")) %dopar% {
#      # Load necessary functions
#      source   ("r_functions.R")
#      source   ("eigen_functions.R")
#      sourceCpp("eigen_functions.cpp")
#      
#      # Identify indices for the current fold
#      na_indices <- folds == l
#      
#      # Extract true values and create masked dataset
#      y_true <- Y[na_indices]
#      Yna <- Y
#      Yna[na_indices] <- NA
#      
#      # Run model with missing values
#      y_ppp <- YPPP(Yna, na_indices, K, n_sams, n_burn, n_skip)
#      
#      # Return results as a named list
#      list(l = l, K = K, y_true = y_true, y_ppp = y_ppp)
# }
# 
# # Stop parallel processing
# stopCluster(cl)
# 
# # save cv
# save(cv, file = "cv_eigen_lazega.RData")

# load cv
load("cv_eigen_lazega.RData")

# Initialize AUCs vector
aucs <- numeric(L)

# Set up ROC plot
pdf(file = "fig_lazega_cv_eigen.pdf", pointsize = 22)
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))

plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "False Positive Rate", ylab = "True Positive Rate", main = "")
abline(a = 0, b = 1, col = "gray", lwd = 2)

# Compute ROC curves and AUCs
for (l in seq_len(L)) {
     # Compute performance metrics
     pred <- prediction(predictions = cv[[l]]$y_ppp, labels = cv[[l]]$y_true)
     perf_roc <- performance(pred, measure = "tpr", x.measure = "fpr")
     
     # Plot ROC curve
     lines(perf_roc@x.values[[1]], perf_roc@y.values[[1]], type = "l", col = 2)
     
     # Compute AUC
     perf_auc <- performance(pred, measure = "auc")
     aucs[l] <- perf_auc@y.values[[1]]
}

legend("bottomright", 
       legend = c(paste0("Mean AUC = ", round(mean(aucs), 3)), paste0("SD AUC = ", round(sd(aucs), 3))), 
       bty = "n")

dev.off()