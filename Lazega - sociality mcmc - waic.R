# Clear the global environment
rm(list = ls(all.names = TRUE))

# Set working directory
setwd("C:/Users/User/Dropbox/PAPERS/projects/sociality")

# Load required libraries and source files
library(igraph)
library(sand)

# Load R functions
source("WAIC.R")
source("r_functions.R")

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
K <- 4
B <- 25000

load("samples_sociality_lazega.RData")

# WAIC

(waic <- WAIC(I, B, Y, mu_chain = samples$mu, delta_chain = samples$delta))

