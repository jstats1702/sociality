# Clear the global environment
rm(list = ls(all.names = TRUE))

# Set working directory
setwd("C:/Users/User/Dropbox/PAPERS/projects/sociality")

# Load required libraries and source files
library(Rcpp)
library(igraph)
library(sand)

# Load compiled C++ functions
sourceCpp("class_functions.cpp")

# Load R functions
source("r_functions.R")
source("class_functions.R")

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
K      <- 10
n_sams <- 25000
n_burn <- 10000
n_skip <- 10

# Run MCMC
start.time <- Sys.time()
set.seed(123)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)
end.time <- Sys.time()
save(samples, start.time, end.time, file = "samples_class_lazega.RData")

load("samples_class_lazega.RData")
end.time - start.time

# test statistics
B <- nrow(samples$Xi_chain)  # Get number of samples
n_test_stats <- 6
test_stats_class <- matrix(NA, B, n_test_stats)  # Preallocate matrix

# set.seed(42)
for (b in seq_len(B)) {
     # Extract parameters
     Lambda <- samples$Lambda_chain[b,]
     Xi     <- samples$Xi_chain[b,]

     # Simulate adjacency matrix and convert to graph
     A <- simulate_data(I, K, Lambda, Xi)
     g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")

     # Compute node degrees once
     degrees <- igraph::degree(g)

     # Compute statistics efficiently
     test_stats_class[b, ] <- c(
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

save(test_stats_class, file = "test_statistics_class_lazega.RData")

load("test_statistics_class_lazega.RData")