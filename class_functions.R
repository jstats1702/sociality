get_hyperpars <- function() {
     list(
          # Hyperparameters
          mu_mu = 0, sig2_mu = 3, a_sig = 3, b_sig = 2, 
          a_alpha = 1, b_alpha = 1,
          
          # Metropolis-Hastings parameters
          del2_Lambda = 1, n_Lambda = 0, n_tun_Lambda = 100,
          del2_alpha = 1, n_alpha = 0, n_tun_alpha = 100
     )
}

get_initial_values <- function(I, K, hyps) {
     alpha  <- 1
     omega  <- rep(alpha / K, K)
     
     Xi <- matrix(c(rep(0, I %% K), rep(seq_len(K) - 1, each = floor(I / K))), ncol = 1)
     
     sigsq  <- 1 / rgamma(1, shape = hyps$a_sig, rate = hyps$b_sig)
     mu     <- rnorm(1, mean = hyps$mu_mu, sd = sqrt(hyps$sig2_mu))
     Lambda <- rnorm(K * (K + 1) / 2, mean = mu, sd = sqrt(sigsq))
     
     list(alpha = alpha, omega = omega, Xi = Xi, sigsq = sigsq, mu = mu, Lambda = Lambda)
}

get_chains_data <- function(I, K, n_sams) {
     list(
          Lambda_chain = matrix(NA_real_, n_sams, K * (K + 1) / 2),
          Xi_chain     = matrix(NA_real_, n_sams, I),
          mu_chain     = matrix(NA_real_, n_sams, 1),
          sigsq_chain  = matrix(NA_real_, n_sams, 1),
          omega_chain  = matrix(NA_real_, n_sams, K),
          alpha_chain  = matrix(NA_real_, n_sams, 1),
          loglik_chain = matrix(NA_real_, n_sams, 1)
     )
}

MCMC <- function(Y, K, n_sams, n_burn, n_skip) {
     # Get data dimensions and initialize storage
     I     <- get_I(Y)
     THETA <- get_chains_data(I, K, n_sams)
     
     # Retrieve hyperparameters and MH parameters
     hyps    <- get_hyperpars()
     mu_mu   <- hyps$mu_mu
     sig2_mu <- hyps$sig2_mu
     a_sig   <- hyps$a_sig
     b_sig   <- hyps$b_sig
     a_alpha <- hyps$a_alpha
     b_alpha <- hyps$b_alpha
     
     del2_Lambda  <- hyps$del2_Lambda
     n_Lambda     <- hyps$n_Lambda
     n_tun_Lambda <- hyps$n_tun_Lambda
     del2_alpha   <- hyps$del2_alpha
     n_alpha      <- hyps$n_alpha
     n_tun_alpha  <- hyps$n_tun_alpha
     
     # Initialize values
     init    <- get_initial_values(I, K, hyps)
     Lambda  <- init$Lambda
     mu      <- init$mu
     sigsq   <- init$sigsq
     Xi      <- init$Xi
     omega   <- init$omega
     alpha   <- init$alpha
     
     # Define total iterations and display frequency
     B      <- n_burn + n_skip * n_sams
     n_disp <- floor(0.1 * B)
     
     # MCMC Sampling
     for (b in 1:B) {
          # Sample parameters
          Lambda_out   <- sample_Lambda(b, n_tun_Lambda, del2_Lambda, n_Lambda, n_burn, I, K, sigsq, mu, Lambda, Xi, Y)
          Lambda       <- Lambda_out$Lambda
          del2_Lambda  <- Lambda_out$del2_Lambda
          n_Lambda     <- Lambda_out$n_Lambda
          n_tun_Lambda <- Lambda_out$n_tun_Lambda
          
          alpha_out    <- sample_alpha(b, n_tun_alpha, del2_alpha, n_alpha, n_burn, K, a_alpha, b_alpha, alpha, omega)
          alpha        <- alpha_out$alpha
          del2_alpha   <- alpha_out$del2_alpha
          n_alpha      <- alpha_out$n_alpha
          n_tun_alpha  <- alpha_out$n_tun_alpha
          
          mu     <- sample_mu(K, mu_mu, sig2_mu, sigsq, Lambda)
          sigsq  <- sample_sigsq(K, a_sig, b_sig, mu, Lambda)
          Xi     <- sample_Xi(I, K, omega, Lambda, Xi, Y)
          omega  <- sample_omega(K, alpha, Xi)
          
          # Store sampled values
          if ((b > n_burn) & (b %% n_skip == 0)) {
               i <- (b - n_burn) / n_skip
               THETA$loglik_chain[i]   <- loglik(I, K, Lambda, Xi, Y)
               THETA$Xi_chain[i, ]     <- Xi
               THETA$Lambda_chain[i, ] <- Lambda
               THETA$alpha_chain[i]    <- alpha
               THETA$mu_chain[i]       <- mu
               THETA$sigsq_chain[i]    <- sigsq
               THETA$omega_chain[i, ]  <- omega
          }
          
          # Display progress
          if (b %% n_disp == 0) {
               cat(
                    "Progress: ", formatC(100 * b / B, digits = 1, format = "f"), "% complete\n",
                    "Lambda mixing rate = ", formatC(100 * n_Lambda / (b * 0.5 * K * (K + 1)), digits = 2, format = "f"), "%, del2_Lambda = ", formatC(del2_Lambda, digits = 5, format = "f"), "\n",
                    "Alpha mixing rate  = ", formatC(100 * n_alpha / b, digits = 2, format = "f"), "%, del2_alpha = ", formatC(del2_alpha, digits = 5, format = "f"), "\n",
                    "---------------------------------------------------\n",
                    sep = ""
               )
          }
     }
     
     return(THETA)
}

YPPP <- function(Yna, na_indices, K, n_sams, n_burn, n_skip) {
     I     <- get_I(Yna)
     hyps  <- get_hyperpars()
     
     # Extract hyperparameters
     mu_mu   <- hyps$mu_mu
     sig2_mu <- hyps$sig2_mu
     a_sig   <- hyps$a_sig
     b_sig   <- hyps$b_sig
     a_alpha <- hyps$a_alpha
     b_alpha <- hyps$b_alpha
     
     # Metropolis-Hastings parameters
     del2_Lambda  <- hyps$del2_Lambda
     n_Lambda     <- hyps$n_Lambda
     n_tun_Lambda <- hyps$n_tun_Lambda
     del2_alpha   <- hyps$del2_alpha
     n_alpha      <- hyps$n_alpha
     n_tun_alpha  <- hyps$n_tun_alpha
     
     # Initial values
     init    <- get_initial_values(I, K, hyps)
     Lambda  <- init$Lambda
     mu      <- init$mu
     sigsq   <- init$sigsq
     Xi      <- init$Xi
     omega   <- init$omega
     alpha   <- init$alpha
     
     # Posterior predictive probabilities
     y_ppp <- rep(0, sum(na_indices))
     
     # Define total iterations
     B <- n_burn + n_skip * n_sams
     
     # MCMC Sampling
     for (b in 1:B) {
          Yna_out     <- sample_Y(I, K, Lambda, Xi, na_indices, Yna)
          Yna         <- Yna_out
          
          Lambda_out  <- sample_Lambda(b, n_tun_Lambda, del2_Lambda, n_Lambda, n_burn, I, K, sigsq, mu, Lambda, Xi, Yna)
          Lambda      <- Lambda_out$Lambda
          del2_Lambda <- Lambda_out$del2_Lambda
          n_Lambda    <- Lambda_out$n_Lambda
          n_tun_Lambda <- Lambda_out$n_tun_Lambda
          
          alpha_out   <- sample_alpha(b, n_tun_alpha, del2_alpha, n_alpha, n_burn, K, a_alpha, b_alpha, alpha, omega)
          alpha       <- alpha_out$alpha
          del2_alpha  <- alpha_out$del2_alpha
          n_alpha     <- alpha_out$n_alpha
          n_tun_alpha <- alpha_out$n_tun_alpha
          
          mu     <- sample_mu(K, mu_mu, sig2_mu, sigsq, Lambda)
          sigsq  <- sample_sigsq(K, a_sig, b_sig, mu, Lambda)
          Xi     <- sample_Xi(I, K, omega, Lambda, Xi, Yna)
          omega  <- sample_omega(K, alpha, Xi)
          
          # Posterior predictive probabilities
          if ((b > n_burn) & (b %% n_skip == 0)) {
               y_ppp <- y_ppp + Yna[na_indices] / n_sams
          }
     }
     
     return(y_ppp)
}

incidence_matrix <- function(THETA) {
     # Extract dimensions
     I <- ncol(THETA$Xi_chain)
     B <- nrow(THETA$Xi_chain)
     K <- (-1 + sqrt(1 + 8 * ncol(THETA$Lambda_chain))) / 2
     
     # Compute incidence matrix vector
     A_vec <- incidence_matrix0(I, K, B, THETA$Xi_chain)
     
     # Initialize symmetric matrix
     A <- matrix(0, I, I)
     
     # Fill upper triangle using nested for loops
     for (i in 1:(I - 1)) {
          for (ii in (i + 1):I) {
               k <- get_k(i, ii, I)
               A[i, ii] <- A_vec[k]
               A[ii, i] <- A[i, ii]
          }
     }
     
     diag(A) <- 1
     
     return(A)
}