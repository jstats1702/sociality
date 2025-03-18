# Full conditional distributions
sample_z <- function(y, mu, delta, z, n) {
     for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
               mean_z <- mu + delta[i] + delta[j]
               if (y[i, j] == 1) {
                    z[i, j] <- truncnorm::rtruncnorm(n = 1, a = 0, b = Inf, mean = mean_z, sd = 1)
               } else {
                    z[i, j] <- truncnorm::rtruncnorm(n = 1, a = -Inf, b = 0, mean = mean_z, sd = 1)
               }
               z[j, i] <- z[i, j]
          }
     }
     return(z)
}
#
sample_mu <- function(z, delta, sigma2) {
     v2_mu <- 1/(1/sigma2 + sum(upper.tri(z)))
     m_mu <- v2_mu * sum(z[upper.tri(z)] - delta[row(z)[upper.tri(z)]] - delta[col(z)[upper.tri(z)]])
     return(rnorm(1, mean = m_mu, sd = sqrt(v2_mu)))
}
#
sample_delta <- function(z, mu, tau2, delta, n) {
     for (i in 1:n) {
          neighbors <- setdiff(1:n, i)
          v2_delta <- 1/(1/tau2 + length(neighbors))
          m_delta <- v2_delta*sum(z[i, neighbors] - mu - delta[neighbors])
          delta[i] <- rnorm(1, mean = m_delta, sd = sqrt(v2_delta))
     }
     delta <- delta - mean(delta)
     return(delta)
}
#
sample_sigma2 <- function(mu) {
     a_sigma_post <- a_sigma + 0.5
     b_sigma_post <- b_sigma + 0.5*mu^2
     return(1 / rgamma(1, shape = a_sigma_post, rate = b_sigma_post))
}
#
sample_tau2 <- function(delta, n) {
     a_tau_post <- a_tau + n/2
     b_tau_post <- b_tau + 0.5*sum(delta^2)
     return(1 / rgamma(1, shape = a_tau_post, rate = b_tau_post))
}
#
sample_Y <- function(Yna_matrix, na_indices_matrix, mu, delta, n) {
     for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
               if (na_indices_matrix[i, j]) {
                    prob <- pnorm(mu + delta[i] + delta[j])
                    Yna_matrix[i, j] <- rbinom(1, size = 1, prob = prob)
                    Yna_matrix[j, i] <- Yna_matrix[i, j]
               }
          }
     }
     
     return(Yna_matrix)
}
# MCMC
gibbs_sampler <- function(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau) {
     # Sizes
     n <- nrow(y)
     
     # Initialize
     mu <- 0
     delta <- rnorm(n, 0, 1)
     sigma2 <- 1
     tau2 <- 1
     z <- matrix(0, n, n)
     
     # Storage
     n_samples <- (n_iter - n_burn) / n_thin
     samples <- list(mu = numeric(n_samples), 
                     delta = matrix(0, nrow = n_samples, ncol = n), 
                     sigma2 = numeric(n_samples), 
                     tau2 = numeric(n_samples))
     
     # Sampling
     cat("Initializing Gibbs Sampler...\n")
     for (t in 1:n_iter) {
          # Sample
          z <- sample_z(y, mu, delta, z, n)
          mu <- sample_mu(z, delta, sigma2)
          delta <- sample_delta(z, mu, tau2, delta, n)
          sigma2 <- sample_sigma2(mu)
          tau2 <- sample_tau2(delta, n)
          
          # Store
          if (t > n_burn && (t - n_burn) %% n_thin == 0) {
               idx <- (t - n_burn) / n_thin
               samples$mu[idx] <- mu
               samples$delta[idx, ] <- delta
               samples$sigma2[idx] <- sigma2
               samples$tau2[idx] <- tau2
          }
          
          # Progress
          if (t %% (n_iter/10) == 0) {
               cat(sprintf("Progress: %d%% completed\n", (t/n_iter)*100))
          }
     }
     cat("Sampling completed.\n")
     return(samples)
}

# MCMC
YPPP <- function(Yna_matrix, na_indices_matrix, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau) {
     # Sizes
     n <- nrow(Yna_matrix)
     n_samples <- (n_iter - n_burn) / n_thin
     
     # Initialize
     mu <- 0
     delta <- rnorm(n, 0, 1)
     sigma2 <- 1
     tau2 <- 1
     z <- matrix(0, n, n)
     
     y_ppp <- numeric(sum(na_indices_matrix)/2)

     # Sampling
     for (t in 1:n_iter) {
          # Sample
          Yna_matrix <- sample_Y(Yna_matrix, na_indices_matrix, mu, delta, n)
          z <- sample_z(Yna_matrix, mu, delta, z, n)
          mu <- sample_mu(z, delta, sigma2)
          delta <- sample_delta(z, mu, tau2, delta, n)
          sigma2 <- sample_sigma2(mu)
          # store
          if (t > n_burn && t %% n_thin == 0) {
               k <- 1
               for (i in 1:(n - 1)) {
                    for (j in (i + 1):n) {
                         if (na_indices_matrix[i, j]) {
                              y_ppp[k] <- y_ppp[k] + Yna_matrix[i, j] / n_samples
                              k <- k + 1
                         }
                    }
               }
          }
     }
     
     return(y_ppp)
}