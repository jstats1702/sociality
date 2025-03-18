# compute MC standard errors
mcse <- function(samples) {
     samples <- as.matrix(samples)
     n_eff  <- c(apply(samples, 2, effectiveSize))
     sd_est <- c(apply(samples, 2, sd))
     sd_est/sqrt(n_eff)
}

# Function to compute the log-likelihood
log_likelihood_sociality <- function(y, samples) {
     n <- nrow(y)
     log_lik_samples <- numeric(length(samples$mu))
     for (s in seq_along(samples$mu)) {
          mu <- samples$mu[s]
          delta <- samples$delta[s, ]
          log_lik <- 0
          for (i in 1:(n - 1)) {
               for (j in (i + 1):n) {
                    p_ij <- pnorm(mu + delta[i] + delta[j])
                    log_lik <- log_lik + y[i, j] * log(p_ij + 1e-10) + (1 - y[i, j]) * log(1 - p_ij + 1e-10)
               }
          }
          log_lik_samples[s] <- log_lik
     }
     return(log_lik_samples)
}

# Function to automatically determine the optimal number of clusters
find_optimal_k <- function(delta) {
     max_k <- min(10, length(delta)) # Maximum number of clusters to evaluate
     wss <- numeric(max_k) # Within-cluster sum of squares
     # Calculate WSS for each k
     for (k in 2:max_k) {
          model <- kmeans(delta, centers = k, nstart = 10)
          wss[k] <- model$tot.withinss # Within-cluster sum of squares
     }
     # Elbow method to find the optimal point
     optimal_k <- which(diff(diff(wss)) > 0)[1] + 1
     if (is.na(optimal_k)) optimal_k <- 2 # Default to at least 2 clusters
     return(optimal_k)
}

# Function to calculate co-clustering probabilities
compute_coclustering <- function(samples) {
     n_samples <- length(samples$mu)  # Number of samples
     n <- ncol(samples$delta)         # Number of nodes
     coclustering <- matrix(0, n, n)  # Initialize co-clustering matrix
     # Iterate over each sample
     for (s in 1:n_samples) {
          # Get the deltas for the current sample
          delta <- samples$delta[s, ]
          # Automatically determine the optimal number of clusters
          k <- find_optimal_k(delta)
          # Apply k-means with the optimal k
          clusters <- kmeans(delta, centers = k, nstart = 10)$cluster
          # Update the co-clustering matrix
          for (i in 1:(n - 1)) {
               for (j in (i + 1):n) {
                    if (clusters[i] == clusters[j]) {
                         coclustering[i, j] <- coclustering[i, j] + 1
                         coclustering[j, i] <- coclustering[i, j] # Symmetry
                    }
               }
          }
     }
     # Average over the total number of samples
     coclustering <- coclustering / n_samples
     return(coclustering)
}

# Point estimation of clusters using Mclust
estimate_clusters_mclust <- function(coclustering_probs) {
     # Apply Mclust directly on the probabilities
     mclust_result <- Mclust(coclustering_probs, verbose = F) # Gaussian mixture model
     clusters <- mclust_result$classification # Cluster labels
     return(clusters)
}

## Create binary matrix E for same cluster membership
compute_binary_E <- function(clusters) {
     n <- length(clusters)
     E <- matrix(0, n, n)
     for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
               if (clusters[i] == clusters[j]) {
                    E[i, j] <- 1
                    E[j, i] <- 1 # Symmetry
               }
          }
     }
     return(E)
}

# Function to reorder matrices based on estimated clusters
reorder_matrix <- function(matrix_data, cluster_labels) {
     order <- order(cluster_labels)
     return(matrix_data[order, order])
}

# Function to calculate the matrix of probabilities theta_ij
compute_theta <- function(samples) {
     n_samples <- length(samples$mu)
     n <- ncol(samples$delta)
     theta_avg <- matrix(0, n, n) # Initialize the average matrix
     # Iterate over each sample
     for (s in 1:n_samples) {
          mu <- samples$mu[s]
          delta <- samples$delta[s, ]
          theta <- matrix(0, n, n)
          # Calculate theta_ij for each pair (i, j)
          for (i in 1:(n - 1)) {
               for (j in (i + 1):n) {
                    theta[i, j] <- pnorm(mu + delta[i] + delta[j])
                    theta[j, i] <- theta[i, j]  # Symmetry
               }
          }
          # Average over the samples
          theta_avg <- theta_avg + theta/n_samples
     }
     return(theta_avg)
}

# Function to visualize test statistics with 99% and 95% credible intervals
test_statistic_viz <- function(k, test_statistic, obs_test_statistics, 
                               test_stats_sociality, test_stats_distance, 
                               test_stats_class, test_stats_eigen) {
     
     ## Compute 99% and 95% credible intervals and posterior means
     ci99_sociality <- quantile(test_stats_sociality[, k], probs = c(0.005, 0.995), na.rm = TRUE)
     ci99_distance  <- quantile(test_stats_distance [, k], probs = c(0.005, 0.995), na.rm = TRUE)
     ci99_class     <- quantile(test_stats_class    [, k], probs = c(0.005, 0.995), na.rm = TRUE)
     ci99_eigen     <- quantile(test_stats_eigen    [, k], probs = c(0.005, 0.995), na.rm = TRUE)
     
     ci95_sociality <- quantile(test_stats_sociality[, k], probs = c(0.025, 0.975), na.rm = TRUE)
     ci95_distance  <- quantile(test_stats_distance [, k], probs = c(0.025, 0.975), na.rm = TRUE)
     ci95_class     <- quantile(test_stats_class    [, k], probs = c(0.025, 0.975), na.rm = TRUE)
     ci95_eigen     <- quantile(test_stats_eigen    [, k], probs = c(0.025, 0.975), na.rm = TRUE)
     
     mean_sociality <- mean(test_stats_sociality[, k], na.rm = TRUE)
     mean_distance  <- mean(test_stats_distance [, k], na.rm = TRUE)
     mean_class     <- mean(test_stats_class    [, k], na.rm = TRUE)
     mean_eigen     <- mean(test_stats_eigen    [, k], na.rm = TRUE)
     
     ## Define y positions for the models
     model_names <- c("S", "D", "C", "E")
     y_positions <- c(4, 3, 2, 1)
     x_range <- range(
          obs_test_statistics[k],
          ci99_sociality, 
          ci99_distance,
          ci99_class,
          ci99_eigen
     )
     
     ## Set up an empty plot
     plot(NA, NA, xlim = x_range, ylim = c(0.5, 4.5), xlab = test_statistic, ylab = "", yaxt = "n", main = "")
     
     ## Add 99% credible intervals (thicker lines)
     segments(ci99_sociality[1], y_positions[1], ci99_sociality[2], y_positions[1], lwd = 1)
     segments(ci99_distance[1],  y_positions[2], ci99_distance[2],  y_positions[2], lwd = 1)
     segments(ci99_class[1],     y_positions[3], ci99_class[2],     y_positions[3], lwd = 1)
     segments(ci99_eigen[1],     y_positions[4], ci99_eigen[2],     y_positions[4], lwd = 1)
     
     ## Add 95% credible intervals (thinner lines)
     segments(ci95_sociality[1], y_positions[1], ci95_sociality[2], y_positions[1], lwd = 4)
     segments(ci95_distance[1],  y_positions[2], ci95_distance[2],  y_positions[2], lwd = 4)
     segments(ci95_class[1],     y_positions[3], ci95_class[2],     y_positions[3], lwd = 4)
     segments(ci95_eigen[1],     y_positions[4], ci95_eigen[2],     y_positions[4], lwd = 4)
     
     ## Add posterior means as black dots
     points(mean_sociality, y_positions[1], pch = 16, col = 1, cex = 1.3)
     points(mean_distance,  y_positions[2], pch = 16, col = 1, cex = 1.3)
     points(mean_class,     y_positions[3], pch = 16, col = 1, cex = 1.3)
     points(mean_eigen,     y_positions[4], pch = 16, col = 1, cex = 1.3)
     
     ## Add the observed statistic as a red dashed line
     abline(v = obs_test_statistics[k], lty = 2, col = 2, lwd = 2)
     
     ## Add model labels on the y-axis
     axis(2, at = y_positions, labels = model_names, las = 1)
}

# Function to simulate data based on the model
simulate_data <- function(n, mu_true, tau2_true) {
     # Generate true node effects
     delta_true <- rnorm(n, mean = 0, sd = sqrt(tau2_true))  
     
     # Initialize latent variables and adjacency matrix
     z_true <- matrix(0, n, n)  
     y <- matrix(0, n, n)  
     
     # Generate the adjacency matrix
     for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
               # Sample latent variable from normal distribution
               z_val <- rnorm(1, mean = mu_true + delta_true[i] + delta_true[j], sd = 1)
               z_true[i, j] <- z_true[j, i] <- z_val  # Ensure symmetry
               
               # Apply threshold to determine link presence
               y[i, j] <- y[j, i] <- as.integer(z_val > 0)  # Ensure symmetry
          }
     }
     
     # Return the simulated data
     list(y = y, z_true = z_true, delta_true = delta_true, mu_true = mu_true, tau2_true = tau2_true)
}

remove_outliers <- function(x)
{
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  x[!(x < lower_bound | x > upper_bound)]
}
