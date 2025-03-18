WAIC <- function(I, B, Y, mu_chain, delta_chain) {
     
     lppd <- 0.0
     slp <- 0.0
     pWAIC2 <- 0.0
     
     for (i in 1:(I - 1)) {
          for (ii in (i + 1):I) {
               
               m <- get_k(i, ii, I)
               
               tmp <- 0.0
               sum_a_ib_sq <- 0.0
               sum_a_ib <- 0.0
               
               for (b in 1:B) {
                    mu <- mu_chain[b]
                    delta <- delta_chain[b, ]

                    a_ib <- dbinom(Y[m], 1, pnorm(mu + delta[i] + delta[ii]), log = TRUE)
                    
                    # WAIC computations
                    tmp <- tmp + exp(a_ib) / B
                    slp <- slp + a_ib / B
                    sum_a_ib_sq <- sum_a_ib_sq + a_ib^2
                    sum_a_ib <- sum_a_ib + a_ib
               }
               
               lppd <- lppd + log(tmp)
               pWAIC2 <- pWAIC2 + (sum_a_ib_sq - B * (sum_a_ib / B)^2) / (B - 1)
          }
     }
     
     pWAIC1 <- 2.0 * (lppd - slp)
     waic1 <- -2.0 * (lppd - pWAIC1)
     waic2 <- -2.0 * (lppd - pWAIC2)
     
     return(list(lppd = lppd,
                 pWAIC1 = pWAIC1,
                 pWAIC2 = pWAIC2,
                 waic1 = waic1,
                 waic2 = waic2))
}