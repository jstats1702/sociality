# Format a number to `k` decimal places
dec <- function(x, k) formatC(round(x, k), format = "f", digits = k)

# Procrustes alignment: Align Z to Z0
procus <- function(Z, Z0, K) {
     # Center Z by aligning its columns with Z0
     Z <- sweep(Z, 2, colMeans(Z)) + colMeans(Z0)
     
     # Compute alignment matrix
     A <- t(Z) %*% (Z0 %*% t(Z0)) %*% Z
     eA <- eigen(A, symmetric = TRUE)
     Ahalf <- eA$vectors[, 1:K] %*% diag(sqrt(eA$values[1:K])) %*% t(eA$vectors[, 1:K])
     
     # Align Z to Z0 using Procrustes transformation
     t(t(Z0) %*% Z %*% solve(Ahalf) %*% t(Z))
}

# Infer the number of nodes (I) from a binary adjacency matrix (Y)
get_I <- function(Y) (1 + sqrt(1 + 8 * nrow(Y))) / 2

# Convert (i, ii) pair to a flattened index for a matrix of size I
get_k <- function(i, ii, I) I * (I - 1) / 2 - (I - i + 1) * (I - i) / 2 + ii - i

# Convert (k, kk) pair to a flattened index for a diagonal matrix of size K
get_k_diag <- function(k, kk, K) K * (k - 1) + kk - (k - 1) * k / 2

# Generate fold assignments for cross-validation
get_folds <- function(M, J, L) {
     folds <- matrix(NA_integer_, nrow = M, ncol = J)
     
     for (j in seq_len(J)) {
          fold_assignments <- rep(seq_len(L), length.out = M)  # Ensures balanced assignment
          folds[, j] <- sample(fold_assignments, M, replace = FALSE)  # Shuffle assignments
     }
     
     folds
}
