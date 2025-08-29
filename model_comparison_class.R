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
source("class_functions.R")
source("r_functions.R")
source("test_functions.R")
source("cv_functions.R")

## Load C++ functions
sourceCpp("class_functions.cpp")

## MCMC settings
K      <- 10
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
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "bktec"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "foot"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "lazega"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "hitech"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "kaptail"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "bkham"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "dol"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "glossgt"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "lesmis"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "salter"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "polbooks"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "adjnoun"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "football"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "nine"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "gen"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "fblog"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "jazz"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "partner"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "indus"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))

#-------------------------------------------------------------------------------

dataset <- "science"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- MCMC(Y, K, n_sams, n_burn, n_skip)

# test
set.seed(42)
ppps <- class_test(Ycube, samples)

# waic
waic <- WAIC(I, K, B = n_sams, Y, Lambda_chain = samples$Lambda_chain, Xi_chain = samples$Xi_chain)

# cv
set.seed(42)
aucs <- class_cv(Y, K, n_sams, n_burn, n_skip)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_class.RData"))