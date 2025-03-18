# Lazega

# settings ---------------------------------------------------------------------

setwd("C:/Users/Juan Camilo/Dropbox/PAPERS/projects/sociality")
# setwd("C:/Users/User/Dropbox/PAPERS/projects/sociality")

rm(list = ls())

library(igraph)
library(sand)
library(corrplot)

source("MCMC.R")
source("VI.R")
source("helper functions.R")
source("r_functions.R")

# data -------------------------------------------------------------------------

g <- graph_from_data_frame(d = elist.lazega, directed = "F")
y <- as.matrix(as_adjacency_matrix(graph = g, names = F))
n <- nrow(y)

# model fitting using MCMC -----------------------------------------------------

# prior 1: a_sigma = a_tau = 2 and b_sigma = b_tau = 1/3
# prior 2: a_sigma = a_tau = 3 and b_sigma = b_tau = 1/3
# prior 3: a_sigma = a_tau = 2 and b_sigma = 1/2 and b_tau = 1/4
# prior 4: a_sigma = a_tau = 3 and b_sigma = 1/2 and b_tau = 1/4
# prior 5: a_sigma = a_tau = 2 and b_sigma = b_tau = 1
# prior 6: a_sigma = a_tau = 3 and b_sigma = b_tau = 2

# settings
n_iter <- 250000 + 10000
n_burn <- 10000
n_thin <- 10

global_bound <- 3
epsilon <- 1e-06
max_iter <- 1000

# prior 1 ----------------------------------------------------------------------

# hyperparameters
a_sigma <- 2 
b_sigma <- 1/3
a_tau   <- 2 
b_tau   <- 1/3

set.seed(123)
samples <- gibbs_sampler(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

set.seed(123)
variational <- vi_sociality(y, a_sigma, b_sigma, a_tau, b_tau, global_bound, epsilon, max_iter)

save(samples, variational, file = "samples_prior_1_sociality_lazega.RData")

# prior 2 ----------------------------------------------------------------------

# hyperparameters
a_sigma <- 3 
b_sigma <- 1/3
a_tau   <- 3 
b_tau   <- 1/3

set.seed(123)
samples <- gibbs_sampler(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

set.seed(123)
variational <- vi_sociality(y, a_sigma, b_sigma, a_tau, b_tau, global_bound, epsilon, max_iter)

save(samples, variational, file = "samples_prior_2_sociality_lazega.RData")

# prior 3 ----------------------------------------------------------------------

# hyperparameters
a_sigma <- 2 
b_sigma <- 1/2
a_tau   <- 2 
b_tau   <- 1/4

set.seed(123)
samples <- gibbs_sampler(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

set.seed(123)
variational <- vi_sociality(y, a_sigma, b_sigma, a_tau, b_tau, global_bound, epsilon, max_iter)

save(samples, variational, file = "samples_prior_3_sociality_lazega.RData")

# prior 4 ----------------------------------------------------------------------

# hyperparameters
a_sigma <- 3 
b_sigma <- 1/2
a_tau   <- 3 
b_tau   <- 1/4

set.seed(123)
samples <- gibbs_sampler(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

set.seed(123)
variational <- vi_sociality(y, a_sigma, b_sigma, a_tau, b_tau, global_bound, epsilon, max_iter)

save(samples, variational, file = "samples_prior_4_sociality_lazega.RData")

# prior 5 ----------------------------------------------------------------------

# hyperparameters
a_sigma <- 2 
b_sigma <- 1
a_tau   <- 2 
b_tau   <- 1

set.seed(123)
samples <- gibbs_sampler(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

set.seed(123)
variational <- vi_sociality(y, a_sigma, b_sigma, a_tau, b_tau, global_bound, epsilon, max_iter)

save(samples, variational, file = "samples_prior_5_sociality_lazega.RData")

# prior 6 ----------------------------------------------------------------------

# hyperparameters
a_sigma <- 3 
b_sigma <- 2
a_tau   <- 3 
b_tau   <- 2

set.seed(123)
samples <- gibbs_sampler(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

set.seed(123)
variational <- vi_sociality(y, a_sigma, b_sigma, a_tau, b_tau, global_bound, epsilon, max_iter)

save(samples, variational, file = "samples_prior_6_sociality_lazega.RData")

# comparison -------------------------------------------------------------------

load("samples_prior_1_sociality_lazega.RData")

x_1 <- colMeans(samples$delta)
y_1 <- variational$mu_delta

load("samples_prior_2_sociality_lazega.RData")

x_2 <- colMeans(samples$delta)
y_2 <- variational$mu_delta

load("samples_prior_3_sociality_lazega.RData")

x_3 <- colMeans(samples$delta)
y_3 <- variational$mu_delta

load("samples_prior_4_sociality_lazega.RData")

x_4 <- colMeans(samples$delta)
y_4 <- variational$mu_delta

load("samples_prior_5_sociality_lazega.RData")

x_5 <- colMeans(samples$delta)
y_5 <- variational$mu_delta

load("samples_prior_6_sociality_lazega.RData")

x_6 <- colMeans(samples$delta)
y_6 <- variational$mu_delta

# delta data under different priors
D <- cbind(
  x_1, x_2, x_3, x_4, x_5, x_6,
  y_1, y_2, y_3, y_4, y_5, y_6
)

colnames(D) <- c(paste("MCMC",1:6), paste("VI",1:6))

# correlation matrix
R <- round(cor(D), 3)

summary(R[upper.tri(R)])

pdf(file = "fig_lazega_sensitivity.pdf", pointsize = 18)
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
corrplot(R, method = "circle", type = "upper", 
         cl.cex = 0.9, tl.cex = 0.9,
         tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("red", "white", "blue"))(200))
dev.off()