# Define hyperparameters
get_hyperpars <- function(I, K) {
     # Calculate hyperparameters
     a_ome <- 3
     b_ome <- 2
     a_sig <- 3
     b_sig <- 2  # (a_sig - 1) * (pi^(K / 2) / exp(lgamma(K / 2 + 1)) * I^(2 / K))
     
     # Metropolis-Hastings parameters
     list(
          a_ome = a_ome, b_ome = b_ome,
          a_sig = a_sig, b_sig = b_sig,
          del2_U = 0.1, n_U = 0, n_tun_U = 100,
          del2_zeta = 0.1, n_zeta = 0, n_tun_zeta = 100
     )
}

# Initialize values for MCMC
get_initial_values <- function(I, K, hyps) {
     omesq <- 1 / rgamma(1, shape = hyps$a_ome, rate = hyps$b_ome)
     sigsq <- 1 / rgamma(1, shape = hyps$a_sig, rate = hyps$b_sig)
     zeta <- rnorm(1, mean = 0, sd = sqrt(omesq))
     U <- matrix(rnorm(I * K, mean = 0, sd = sqrt(sigsq)), nrow = I, ncol = K)
     
     list(omesq = omesq, sigsq = sigsq, zeta = zeta, U = U)
}

# Create chains data structure
get_chains_data <- function(I, K, n_sams) {
     list(
          U_chain = matrix(NA, n_sams, I * K),
          zeta_chain = numeric(n_sams),
          sigsq_chain = numeric(n_sams),
          omesq_chain = numeric(n_sams)
     )
}

# Main MCMC function
MCMC <- function(Y, K, n_sams, n_burn, n_skip) {
     I <- get_I(Y)
     THETA <- get_chains_data(I, K, n_sams)
     
     # Initialize hyperparameters and initial values
     hyps <- get_hyperpars(I, K)
     params <- get_initial_values(I, K, hyps)
     U <- params$U
     zeta <- params$zeta
     sigsq <- params$sigsq
     omesq <- params$omesq
     
     # Metropolis-Hastings parameters
     del2_U <- hyps$del2_U
     n_U <- hyps$n_U
     n_tun_U <- hyps$n_tun_U
     del2_zeta <- hyps$del2_zeta
     n_zeta <- hyps$n_zeta
     n_tun_zeta <- hyps$n_tun_zeta
     
     # Total iterations
     B <- n_burn + n_skip * n_sams
     n_disp <- floor(0.1 * B)
     
     for (b in seq_len(B)) {
          # Update U
          tmp <- sample_U(b, n_tun_U, del2_U, n_U, n_burn, I, K, sigsq, zeta, U, Y)
          U <- tmp$U
          del2_U <- tmp$del2_U
          n_U <- tmp$n_U
          n_tun_U <- tmp$n_tun_U
          
          # Update zeta
          tmp <- sample_zeta(b, n_tun_zeta, del2_zeta, n_zeta, n_burn, I, omesq, zeta, U, Y)
          zeta <- tmp$zeta
          del2_zeta <- tmp$del2_zeta
          n_zeta <- tmp$n_zeta
          n_tun_zeta <- tmp$n_tun_zeta
          
          # Update sigsq and omesq
          sigsq <- sample_sigsq(I, K, hyps$a_sig, hyps$b_sig, U)
          omesq <- sample_omesq(hyps$a_ome, hyps$b_ome, zeta)
          
          # Store samples after burn-in
          if (b > n_burn && b %% n_skip == 0) {
               i <- (b - n_burn) / n_skip
               THETA$U_chain[i, ] <- c(U)
               THETA$zeta_chain[i] <- zeta
               THETA$sigsq_chain[i] <- sigsq
               THETA$omesq_chain[i] <- omesq
          }
          
          # Display progress
          if (b %% n_disp == 0) {
               cat(sprintf(
                    "Progress: %.2f%% complete\nU mixing rate = %.2f%%, del2_U = %.5f\nzeta mixing rate = %.2f%%, del2_zeta = %.5f\n\n",
                    100 * b / B, 100 * n_U / (b * I), del2_U, 100 * n_zeta / b, del2_zeta
               ))
          }
     }
     
     THETA
}

# Posterior predictive probabilities
YPPP <- function(Yna, na_indices, K, n_sams, n_burn, n_skip) {
     I <- get_I(Yna)
     hyps <- get_hyperpars(I, K)
     params <- get_initial_values(I, K, hyps)
     U <- params$U
     zeta <- params$zeta
     sigsq <- params$sigsq
     omesq <- params$omesq
     
     del2_U <- hyps$del2_U
     n_U <- hyps$n_U
     del2_zeta <- hyps$del2_zeta
     n_tun_U <- hyps$n_tun_U
     n_zeta <- hyps$n_zeta
     n_tun_zeta <- hyps$n_tun_zeta
     
     y_ppp <- numeric(sum(na_indices))
     B <- n_burn + n_skip * n_sams
     
     for (b in seq_len(B)) {
          Yna <- sample_Y(I, zeta, U, na_indices, Yna)
          tmp <- sample_U(b, n_tun_U, del2_U, n_U, n_burn, I, K, sigsq, zeta, U, Yna)
          U <- tmp$U
          del2_U <- tmp$del2_U
          n_U <- tmp$n_U
          n_tun_U <- tmp$n_tun_U
          
          tmp <- sample_zeta(b, n_tun_zeta, del2_zeta, n_zeta, n_burn, I, omesq, zeta, U, Yna)
          zeta <- tmp$zeta
          del2_zeta <- tmp$del2_zeta
          n_zeta <- tmp$n_zeta
          n_tun_zeta <- tmp$n_tun_zeta
          
          sigsq <- sample_sigsq(I, K, hyps$a_sig, hyps$b_sig, U)
          omesq <- sample_omesq(hyps$a_ome, hyps$b_ome, zeta)
          
          if (b > n_burn && b %% n_skip == 0) {
               y_ppp <- y_ppp + Yna[na_indices] / n_sams
          }
     }
     
     y_ppp
}