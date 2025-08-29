# Bayesian Sociality Models: A Scalable and Flexible Alternative for Network Analysis
#
# Section 5.1.6

# Settings ---------------------------------------------------------------------

## Set working directory
setwd("~/Dropbox/PAPERS/projects/sociality")

## Clear global environment
rm(list = ls(all.names = TRUE))

## Load required libraries and source files
library(igraph)
library(sand)

## Load R functions
source("WAIC.R")
source("r_functions.R")

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

## Run Lazega_sociality_mcmc.R

## Posterior draws
B <- 25000

## Load posterior draws
load("samples_sociality_lazega.RData")

(waic <- WAIC(I, B, Y, mu_chain = samples$mu, delta_chain = samples$delta))

# End --------------------------------------------------------------------------
