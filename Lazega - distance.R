# Clear the global environment
rm(list = ls(all.names = TRUE))

# Set working directory
setwd("C:/Users/User/Dropbox/PAPERS/projects/sociality")

# Load required libraries and source files
library(Rcpp)
library(igraph)
library(sand)

# Load compiled C++ functions
sourceCpp("dist_functions.cpp")

# Load R functions
source("r_functions.R")
source("dist_functions.R")

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

# Model fitting settings
K      <- 4
n_sams <- 25000
n_burn <- 10000
n_skip <- 10

# Run MCMC
start.time <- Sys.time()
set.seed(123)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)
end.time <- Sys.time()
save(samples, start.time, end.time, file = "samples_distance_lazega.RData")

load("samples_distance_lazega.RData")
end.time - start.time

# test statistics
B <- nrow(samples$U_chain)  # Get number of samples
n_test_stats <- 6
test_stats_distance <- matrix(NA, B, n_test_stats)  # Preallocate matrix

set.seed(42)
for (b in seq_len(B)) {
     # Extract parameters
     zeta <- samples$zeta_chain[b]
     U <- matrix(samples$U_chain[b, ], I, K)

     # Simulate adjacency matrix and convert to graph
     A <- simulate_data(I, zeta, U)
     g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")

     # Compute node degrees once
     degrees <- igraph::degree(g)

     # Compute statistics efficiently
     test_stats_distance[b, ] <- c(
          igraph::edge_density(g),
          igraph::transitivity(g, type = "global"),
          igraph::assortativity_degree(g),
          igraph::mean_distance(g),
          mean(degrees),
          sd(degrees)
     )

     # Display progress every 10%
     if (b %% ceiling(0.1 * B) == 0 || b == B) {
          cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
     }
}

save(test_stats_distance, file = "test_statistics_distance_lazega.RData")

load("test_statistics_distance_lazega.RData")