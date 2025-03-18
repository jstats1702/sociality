distance_test <- function(Ycube, samples) {
     # settings
     Rcpp::sourceCpp("dist_functions.cpp")
     
     # observed statistics
     g <- igraph::graph_from_adjacency_matrix(Ycube, mode = "undirected")
     dens0  <- igraph::edge_density(g)
     trans0 <- igraph::transitivity(g, type = "global")
     assor0 <- igraph::assortativity_degree(g)
     
     B <- nrow(samples$U_chain)
     n_test_stats <- 3
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
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
          test_stats[b, ] <- c(
               igraph::edge_density(g),
               igraph::transitivity(g, type = "global"),
               igraph::assortativity_degree(g)
          )

          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     return(c(
          mean(test_stats[,1] < dens0 ),
          mean(test_stats[,2] < trans0),
          mean(test_stats[,3] < assor0)
     ))
}

class_test <- function(Ycube, samples) {
     # settings
     Rcpp::sourceCpp("class_functions.cpp")
     
     # observed statistics
     g <- igraph::graph_from_adjacency_matrix(Ycube, mode = "undirected")
     dens0  <- igraph::edge_density(g)
     trans0 <- igraph::transitivity(g, type = "global")
     assor0 <- igraph::assortativity_degree(g)
     
     B <- nrow(samples$Xi_chain)
     n_test_stats <- 3
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
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
          test_stats[b, ] <- c(
               igraph::edge_density(g),
               igraph::transitivity(g, type = "global"),
               igraph::assortativity_degree(g)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     return(c(
          mean(test_stats[,1] < dens0 ),
          mean(test_stats[,2] < trans0),
          mean(test_stats[,3] < assor0)
     ))
}

eigen_test <- function(Ycube, samples) {
     # settings
     Rcpp::sourceCpp("eigen_functions.cpp")
     
     # observed statistics
     g <- igraph::graph_from_adjacency_matrix(Ycube, mode = "undirected")
     dens0  <- igraph::edge_density(g)
     trans0 <- igraph::transitivity(g, type = "global")
     assor0 <- igraph::assortativity_degree(g)
     
     B <- nrow(samples$Lambda_chain)
     n_test_stats <- 3
     test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
     
     for (b in seq_len(B)) {
          # Extract parameters
          zeta   <- samples$zeta_chain[b]
          U      <- matrix(samples$U_chain[b,], I, K)
          Lambda <- samples$Lambda[b,]
          
          # Simulate adjacency matrix and convert to graph
          A <- simulate_data(I, zeta, U, Lambda)
          g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
          
          # Compute node degrees once
          degrees <- igraph::degree(g)
          
          # Compute statistics efficiently
          test_stats[b, ] <- c(
               igraph::edge_density(g),
               igraph::transitivity(g, type = "global"),
               igraph::assortativity_degree(g)
          )
          
          # Display progress every 10%
          if (b %% ceiling(0.1 * B) == 0 || b == B) {
               cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
          }
     }
     
     return(c(
          mean(test_stats[,1] < dens0 ),
          mean(test_stats[,2] < trans0),
          mean(test_stats[,3] < assor0)
     ))
}

sociality_test <- function(Ycube, samples) {
  # settings
  n <- nrow(Ycube)

  # observed statistics
  g <- igraph::graph_from_adjacency_matrix(Ycube, mode = "undirected")
  dens0  <- igraph::edge_density(g)
  trans0 <- igraph::transitivity(g, type = "global")
  assor0 <- igraph::assortativity_degree(g)
  
  B <- nrow(samples$delta)
  n_test_stats <- 3
  test_stats <- matrix(NA, B, n_test_stats)  # Preallocate matrix
  
  for (b in seq_len(B)) {
    # Extract parameters
    mu <- samples$mu[b]
    delta <- samples$delta[b, ]

    # Simulate adjacency matrix and convert to graph
    A <- matrix(0, n, n)
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
             A[i, j] <- rbinom(1, size = 1, prob = pnorm(mu + delta[i] + delta[j]))
             A[j, i] <- A[i, j]
          }
     }
    g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
    
    # Compute node degrees once
    degrees <- igraph::degree(g)
    
    # Compute statistics efficiently
    test_stats[b, ] <- c(
      igraph::edge_density(g),
      igraph::transitivity(g, type = "global"),
      igraph::assortativity_degree(g)
    )
    
    # Display progress every 10%
    if (b %% ceiling(0.1 * B) == 0 || b == B) {
      cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
    }
  }
  
  return(c(
    mean(test_stats[,1] < dens0 ),
    mean(test_stats[,2] < trans0),
    mean(test_stats[,3] < assor0)
  ))
}