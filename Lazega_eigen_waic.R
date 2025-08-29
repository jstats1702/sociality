# Bayesian Sociality Models: A Scalable and Flexible Alternative for Network Analysis
#
# Section 5.1

# Settings ---------------------------------------------------------------------

## Working directory
setwd("~/Dropbox/PAPERS/projects/sociality")

## Clean global environment
rm(list = ls())

## Load required libraries and source files
library(Rcpp)
library(igraph)
library(sand)

## Load compiled C++ functions
sourceCpp("eigen_functions.cpp")

## Load R functions
source("r_functions.R")
source("eigen_functions.R")

# Data -------------------------------------------------------------------------

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

# WAIC computation -------------------------------------------------------------

## Run Lazega_eigen.R

## Model fitting settings
K <- 4
B <- 25000

load("samples_eigen_lazega.RData")

(waic <- WAIC(I, K, B, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain, Lambda_chain = samples$Lambda_chain))

# End --------------------------------------------------------------------------