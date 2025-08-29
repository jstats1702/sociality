# Bayesian Sociality Models: A Scalable and Flexible Alternative for Network Analysis
#
# Section 5.2

# Settings ---------------------------------------------------------------------

## Working directory
setwd("~/Dropbox/PAPERS/projects/sociality")

## Clean global environment
rm(list = ls())

## Required libraries
library(Rcpp)

## Load R functions
source("dist_functions.R")
source("r_functions.R")
source("test_functions.R")
source("cv_functions.R")

## Load C++ functions
sourceCpp("dist_functions.cpp")

## MCMC settings
K      <- 4
n_sams <- 25000
n_burn <- 10000
n_skip <- 10

#-------------------------------------------------------------------------------

dataset <- "zach"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "bktec"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "foot"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "lazega"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "hitech"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "kaptail"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "bkham"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "dol"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "glossgt"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "lesmis"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "salter"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "polbooks"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "adjnoun"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "football"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "nine"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "gen"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "fblog"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "jazz"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "partner"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "indus"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))

#-------------------------------------------------------------------------------

dataset <- "science"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- distance_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, zeta_chain = samples$zeta_chain, U_chain = samples$U_chain)

# cv
set.seed(42)
aucs <- distance_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_distance.RData"))