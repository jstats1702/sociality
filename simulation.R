# Bayesian Sociality Models: A Scalable and Flexible Alternative for Network Analysis
#
# Section 5.3

# Settings ---------------------------------------------------------------------

## Working directory
setwd("~/Dropbox/PAPERS/projects/sociality")

## Clean global environment
rm(list = ls())

## Required libraries
library(igraph)
library(sand)
library(doParallel)
library(foreach)

## Load R functions
source("MCMC.R")
source("VI.R")
source("helper functions.R")
source("r_functions.R")

# Simulation settings-----------------------------------------------------------

# Number of simulations
N <- 100

# model parameters
mu_true <- -2
tau2_true <- 0.5

# hyperparameters
a_sigma <- 2 
b_sigma <- 1/3
a_tau   <- 2 
b_tau   <- 1/3

# model fit settings
n_iter <- 25000 + 10000
n_burn <- 10000
n_thin <- 1

global_bound <- 3
epsilon <- 1e-06
max_iter <- 1000

# Simulation 1 -----------------------------------------------------------------

n <- 25

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Set random seed for reproducibility
set.seed(420)

# Run simulations in parallel
results <- foreach(i = 1:N, .inorder = F, .combine = "rbind") %dopar% {
     
     # Data generation
     sim_data <- simulate_data(n, mu_true, tau2_true)
     y <- sim_data$y
     
     # Model fitting: Gibbs sampling
     start.time <- Sys.time()
     samples <- gibbs_sampler(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)
     end.time <- Sys.time()
     time_mcmc <- as.numeric(end.time - start.time)
     
     # Model fitting: Variational inference
     start.time <- Sys.time()
     variational <- vi_sociality(y, a_sigma, b_sigma, a_tau, b_tau, global_bound, epsilon, max_iter)
     end.time <- Sys.time()
     time_vi <- as.numeric(end.time - start.time)
     
     # RMSE
     delta_true <- samples$delta
     delta_hat_mcmc <- colMeans(samples$delta)
     delta_hat_vi <- variational$mu_delta
     rmse_mcmc <- sqrt(mean((delta_hat_mcmc - delta_true)^2))
     rmse_vi <- sqrt(mean((delta_hat_vi - delta_true)^2))
     
     # Coverage
     lower_bound_mcmc <- apply(samples$delta, 2, quantile, probs = 0.025)
     upper_bound_mcmc <- apply(samples$delta, 2, quantile, probs = 0.975)
     coverage_mcmc <- mean((delta_true >= lower_bound_mcmc) & (delta_true <= upper_bound_mcmc))
     
     lower_bound_vi <- variational$mu_delta + qnorm(0.025) * sqrt(variational$sigma_delta2)
     upper_bound_vi <- variational$mu_delta + qnorm(0.975) * sqrt(variational$sigma_delta2)
     coverage_vi <- mean((delta_true >= lower_bound_vi) & (delta_true <= upper_bound_vi))
     
     metrics <- c(time_mcmc, time_vi, rmse_mcmc, rmse_vi, coverage_mcmc, coverage_vi)
     
     return(metrics)  # Optionally collect results
}

# Stop parallel backend
stopCluster(cl)

save(results, file = "simulation_sociality_1.RData")

# Simulation 2 -----------------------------------------------------------------

n <- 50

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Set random seed for reproducibility
set.seed(420)

# Run simulations in parallel
results <- foreach(i = 1:N, .inorder = F, .combine = "rbind") %dopar% {
     
     # Data generation
     sim_data <- simulate_data(n, mu_true, tau2_true)
     y <- sim_data$y
     
     # Model fitting: Gibbs sampling
     start.time <- Sys.time()
     samples <- gibbs_sampler(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)
     end.time <- Sys.time()
     time_mcmc <- as.numeric(end.time - start.time)
     
     # Model fitting: Variational inference
     start.time <- Sys.time()
     variational <- vi_sociality(y, a_sigma, b_sigma, a_tau, b_tau, global_bound, epsilon, max_iter)
     end.time <- Sys.time()
     time_vi <- as.numeric(end.time - start.time)
     
     # RMSE
     delta_true <- samples$delta
     delta_hat_mcmc <- colMeans(samples$delta)
     delta_hat_vi <- variational$mu_delta
     rmse_mcmc <- sqrt(mean((delta_hat_mcmc - delta_true)^2))
     rmse_vi <- sqrt(mean((delta_hat_vi - delta_true)^2))
     
     # Coverage
     lower_bound_mcmc <- apply(samples$delta, 2, quantile, probs = 0.025)
     upper_bound_mcmc <- apply(samples$delta, 2, quantile, probs = 0.975)
     coverage_mcmc <- mean((delta_true >= lower_bound_mcmc) & (delta_true <= upper_bound_mcmc))
     
     lower_bound_vi <- variational$mu_delta + qnorm(0.025) * sqrt(variational$sigma_delta2)
     upper_bound_vi <- variational$mu_delta + qnorm(0.975) * sqrt(variational$sigma_delta2)
     coverage_vi <- mean((delta_true >= lower_bound_vi) & (delta_true <= upper_bound_vi))
     
     metrics <- c(time_mcmc, time_vi, rmse_mcmc, rmse_vi, coverage_mcmc, coverage_vi)
     
     return(metrics)  # Optionally collect results
}

# Stop parallel backend
stopCluster(cl)

save(results, file = "simulation_sociality_2.RData")

# Simulation 3 -----------------------------------------------------------------

n <- 100

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Set random seed for reproducibility
set.seed(420)

# Run simulations in parallel
results <- foreach(i = 1:N, .inorder = F, .combine = "rbind") %dopar% {
     
     # Data generation
     sim_data <- simulate_data(n, mu_true, tau2_true)
     y <- sim_data$y
     
     # Model fitting: Gibbs sampling
     start.time <- Sys.time()
     samples <- gibbs_sampler(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)
     end.time <- Sys.time()
     time_mcmc <- as.numeric(end.time - start.time)
     
     # Model fitting: Variational inference
     start.time <- Sys.time()
     variational <- vi_sociality(y, a_sigma, b_sigma, a_tau, b_tau, global_bound, epsilon, max_iter)
     end.time <- Sys.time()
     time_vi <- as.numeric(end.time - start.time)
     
     # RMSE
     delta_true <- samples$delta
     delta_hat_mcmc <- colMeans(samples$delta)
     delta_hat_vi <- variational$mu_delta
     rmse_mcmc <- sqrt(mean((delta_hat_mcmc - delta_true)^2))
     rmse_vi <- sqrt(mean((delta_hat_vi - delta_true)^2))
     
     # Coverage
     lower_bound_mcmc <- apply(samples$delta, 2, quantile, probs = 0.025)
     upper_bound_mcmc <- apply(samples$delta, 2, quantile, probs = 0.975)
     coverage_mcmc <- mean((delta_true >= lower_bound_mcmc) & (delta_true <= upper_bound_mcmc))
     
     lower_bound_vi <- variational$mu_delta + qnorm(0.025) * sqrt(variational$sigma_delta2)
     upper_bound_vi <- variational$mu_delta + qnorm(0.975) * sqrt(variational$sigma_delta2)
     coverage_vi <- mean((delta_true >= lower_bound_vi) & (delta_true <= upper_bound_vi))
     
     metrics <- c(time_mcmc, time_vi, rmse_mcmc, rmse_vi, coverage_mcmc, coverage_vi)
     
     return(metrics)  # Optionally collect results
}

# Stop parallel backend
stopCluster(cl)

save(results, file = "simulation_sociality_3.RData")

# Simulation 4 -----------------------------------------------------------------

n <- 200

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Set random seed for reproducibility
set.seed(420)

# Run simulations in parallel
results <- foreach(i = 1:N, .inorder = F, .combine = "rbind") %dopar% {
     
     # Data generation
     sim_data <- simulate_data(n, mu_true, tau2_true)
     y <- sim_data$y
     
     # Model fitting: Gibbs sampling
     start.time <- Sys.time()
     samples <- gibbs_sampler(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)
     end.time <- Sys.time()
     time_mcmc <- as.numeric(end.time - start.time)
     
     # Model fitting: Variational inference
     start.time <- Sys.time()
     variational <- vi_sociality(y, a_sigma, b_sigma, a_tau, b_tau, global_bound, epsilon, max_iter)
     end.time <- Sys.time()
     time_vi <- as.numeric(end.time - start.time)
     
     # RMSE
     delta_true <- samples$delta
     delta_hat_mcmc <- colMeans(samples$delta)
     delta_hat_vi <- variational$mu_delta
     rmse_mcmc <- sqrt(mean((delta_hat_mcmc - delta_true)^2))
     rmse_vi <- sqrt(mean((delta_hat_vi - delta_true)^2))
     
     # Coverage
     lower_bound_mcmc <- apply(samples$delta, 2, quantile, probs = 0.025)
     upper_bound_mcmc <- apply(samples$delta, 2, quantile, probs = 0.975)
     coverage_mcmc <- mean((delta_true >= lower_bound_mcmc) & (delta_true <= upper_bound_mcmc))
     
     lower_bound_vi <- variational$mu_delta + qnorm(0.025) * sqrt(variational$sigma_delta2)
     upper_bound_vi <- variational$mu_delta + qnorm(0.975) * sqrt(variational$sigma_delta2)
     coverage_vi <- mean((delta_true >= lower_bound_vi) & (delta_true <= upper_bound_vi))
     
     metrics <- c(time_mcmc, time_vi, rmse_mcmc, rmse_vi, coverage_mcmc, coverage_vi)
     
     return(metrics)  # Optionally collect results
}

# Stop parallel backend
stopCluster(cl)

save(results, file = "simulation_sociality_4.RData")

# Results ----------------------------------------------------------------------

#-------------#
#   Table 4   #
#-------------#

out <- NULL

for (i in 1:4) {
  load(paste0("simulation_sociality_", i, ".RData"))
  if (i == 4) {  # results in seconds instead of hours
    results[ , 1] <- 60 * results[ , 1]
    results[ , 2] <- 60 * results[ , 2]
  }
  tmp <- NULL
  for (j in 1:4) {
    tmp[j] <- mean(remove_outliers(results[ , j]))
  }
  out <- rbind(out, tmp) 
}

rownames(out) <- paste0("Scenario ", 1:4)
colnames(out) <- c("Time MCMC", "Time VI", "RMSE MCMC", "RMSE VI")
round(out, 3)

# End --------------------------------------------------------------------------
