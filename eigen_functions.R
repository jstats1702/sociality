get_hyperpars <- function() {
     list(
          # Hyperparameters
          a_ome = 3, b_ome = 2,
          a_sig = 3, b_sig = 2,
          a_kap = 3, b_kap = 2,
          
          # Metropolis-Hastings (MH) parameters
          del2_U = 0.1, n_U = 0, n_tun_U = 100,
          del2_Lambda = 0.1, n_Lambda = 0, n_tun_Lambda = 100,
          del2_zeta = 0.1, n_zeta = 0, n_tun_zeta = 100
     )
}

get_initial_values <- function(I, K, hyps) {
     # Sample variance components
     omesq  <- 1 / rgamma(1, shape = hyps$a_ome, rate = hyps$b_ome)
     kapsq  <- 1 / rgamma(1, shape = hyps$a_kap, rate = hyps$b_kap)
     sigsq  <- 1 / rgamma(1, shape = hyps$a_sig, rate = hyps$b_sig)
     
     # Sample parameters
     zeta   <- rnorm(1, mean = 0, sd = sqrt(omesq))
     Lambda <- rnorm(K, mean = 0, sd = sqrt(kapsq))
     U      <- matrix(rnorm(I * K, mean = 0, sd = sqrt(sigsq)), I, K)
     
     # Return as a structured list
     list(
          omesq = omesq, sigsq = sigsq, kapsq = kapsq,
          zeta = zeta, Lambda = Lambda, U = U
     )
}

get_chains_data <- function(I, K, n_sams) {
     list(
          U_chain      = matrix(NA_real_, n_sams, I * K),
          Lambda_chain = matrix(NA_real_, n_sams, K),
          zeta_chain   = numeric(n_sams),
          sigsq_chain  = numeric(n_sams),
          kapsq_chain  = numeric(n_sams),
          omesq_chain  = numeric(n_sams),
          loglik_chain = numeric(n_sams)
     )
}

MCMC <- function(Y, K, n_sams, n_burn, n_skip) {
     # Number of nodes
     I <- get_I(Y)
     
     # Preallocate storage for MCMC chains
     THETA <- get_chains_data(I, K, n_sams)
     
     # Retrieve hyperparameters
     hyps <- get_hyperpars()
     a_ome <- hyps$a_ome; b_ome <- hyps$b_ome
     a_kap <- hyps$a_kap; b_kap <- hyps$b_kap
     a_sig <- hyps$a_sig; b_sig <- hyps$b_sig
     del2_U <- hyps$del2_U; n_U <- hyps$n_U; n_tun_U <- hyps$n_tun_U
     del2_Lambda <- hyps$del2_Lambda; n_Lambda <- hyps$n_Lambda; n_tun_Lambda <- hyps$n_tun_Lambda
     del2_zeta <- hyps$del2_zeta; n_zeta <- hyps$n_zeta; n_tun_zeta <- hyps$n_tun_zeta
     
     # Initialize values
     init_vals <- get_initial_values(I, K, hyps)
     U <- init_vals$U; Lambda <- init_vals$Lambda; zeta <- init_vals$zeta
     sigsq <- init_vals$sigsq; kapsq <- init_vals$kapsq; omesq <- init_vals$omesq
     
     # Total iterations and display frequency
     B <- n_burn + n_skip * n_sams
     n_disp <- floor(0.1 * B)
     
     # MCMC Sampling
     for (b in 1:B) {
          # Sample U
          U_sample <- sample_U(b, n_tun_U, del2_U, n_U, n_burn, I, K, sigsq, zeta, U, Lambda, Y)
          U <- U_sample$U; del2_U <- U_sample$del2_U; n_U <- U_sample$n_U; n_tun_U <- U_sample$n_tun_U
          
          # Sample Lambda
          Lambda_sample <- sample_Lambda(b, n_tun_Lambda, del2_Lambda, n_Lambda, n_burn, I, K, kapsq, zeta, U, Lambda, Y)
          Lambda <- Lambda_sample$Lambda; del2_Lambda <- Lambda_sample$del2_Lambda
          n_Lambda <- Lambda_sample$n_Lambda; n_tun_Lambda <- Lambda_sample$n_tun_Lambda
          
          # Sample zeta
          zeta_sample <- sample_zeta(b, n_tun_zeta, del2_zeta, n_zeta, n_burn, I, omesq, zeta, U, Lambda, Y)
          zeta <- zeta_sample$zeta; del2_zeta <- zeta_sample$del2_zeta
          n_zeta <- zeta_sample$n_zeta; n_tun_zeta <- zeta_sample$n_tun_zeta
          
          # Sample variance components
          sigsq <- sample_sigsq(I, K, a_sig, b_sig, U)
          kapsq <- sample_kapsq(K, a_kap, b_kap, Lambda)
          omesq <- sample_omesq(a_ome, b_ome, zeta)
          
          # Store samples
          if ((b > n_burn) & (b %% n_skip == 0)) {
               i <- (b - n_burn) / n_skip
               THETA$loglik_chain[i]   <- loglik(I, zeta, U, Lambda, Y)
               THETA$U_chain[i, ]      <- c(U)
               THETA$Lambda_chain[i, ] <- c(Lambda)
               THETA$zeta_chain[i]     <- zeta
               THETA$sigsq_chain[i]    <- sigsq
               THETA$kapsq_chain[i]    <- kapsq
               THETA$omesq_chain[i]    <- omesq
          }
          
          # Progress display
          if (b %% n_disp == 0) {
               cat(
                    "Progress: ", sprintf("%.1f", 100 * b / B), "% complete", "\n",
                    "U      mixing rate = ", sprintf("%.2f", 100 * n_U / (b * I)), "%, del2_U      = ", sprintf("%.5f", del2_U), "\n",
                    "Lambda mixing rate = ", sprintf("%.2f", 100 * n_Lambda / (b * K)), "%, del2_Lambda = ", sprintf("%.5f", del2_Lambda), "\n",
                    "zeta   mixing rate = ", sprintf("%.2f", 100 * n_zeta / b), "%, del2_zeta   = ", sprintf("%.5f", del2_zeta), "\n\n",
                    sep = ""
               )
          }
     }
     
     return(THETA)
}

YPPP <- function(Yna, na_indices, K, n_sams, n_burn, n_skip) {
     # Number of nodes
     I <- get_I(Yna)
     
     # Retrieve hyperparameters
     hyps <- get_hyperpars()
     a_ome <- hyps$a_ome; b_ome <- hyps$b_ome
     a_kap <- hyps$a_kap; b_kap <- hyps$b_kap
     a_sig <- hyps$a_sig; b_sig <- hyps$b_sig
     del2_U <- hyps$del2_U; n_U <- hyps$n_U; n_tun_U <- hyps$n_tun_U
     del2_Lambda <- hyps$del2_Lambda; n_Lambda <- hyps$n_Lambda; n_tun_Lambda <- hyps$n_tun_Lambda
     del2_zeta <- hyps$del2_zeta; n_zeta <- hyps$n_zeta; n_tun_zeta <- hyps$n_tun_zeta
     
     # Initialize values
     init_vals <- get_initial_values(I, K, hyps)
     U <- init_vals$U; Lambda <- init_vals$Lambda; zeta <- init_vals$zeta
     sigsq <- init_vals$sigsq; kapsq <- init_vals$kapsq; omesq <- init_vals$omesq
     
     # Posterior predictive probabilities
     y_ppp <- numeric(sum(na_indices))
     
     # Total iterations
     B <- n_burn + n_skip * n_sams
     
     # MCMC Sampling
     for (b in 1:B) {
          # Update missing values
          Yna <- sample_Y(I, zeta, U, Lambda, na_indices, Yna)
          
          # Sample U
          U_sample <- sample_U(b, n_tun_U, del2_U, n_U, n_burn, I, K, sigsq, zeta, U, Lambda, Yna)
          U <- U_sample$U; del2_U <- U_sample$del2_U; n_U <- U_sample$n_U; n_tun_U <- U_sample$n_tun_U
          
          # Sample Lambda
          Lambda_sample <- sample_Lambda(b, n_tun_Lambda, del2_Lambda, n_Lambda, n_burn, I, K, kapsq, zeta, U, Lambda, Yna)
          Lambda <- Lambda_sample$Lambda; del2_Lambda <- Lambda_sample$del2_Lambda
          n_Lambda <- Lambda_sample$n_Lambda; n_tun_Lambda <- Lambda_sample$n_tun_Lambda
          
          # Sample zeta
          zeta_sample <- sample_zeta(b, n_tun_zeta, del2_zeta, n_zeta, n_burn, I, omesq, zeta, U, Lambda, Yna)
          zeta <- zeta_sample$zeta; del2_zeta <- zeta_sample$del2_zeta
          n_zeta <- zeta_sample$n_zeta; n_tun_zeta <- zeta_sample$n_tun_zeta
          
          # Sample variance components
          sigsq <- sample_sigsq(I, K, a_sig, b_sig, U)
          kapsq <- sample_kapsq(K, a_kap, b_kap, Lambda)
          omesq <- sample_omesq(a_ome, b_ome, zeta)
          
          # Posterior predictive probabilities update
          if ((b > n_burn) & (b %% n_skip == 0)) {
               y_ppp <- y_ppp + Yna[na_indices] / n_sams
          }
     }
     
     return(y_ppp)
}