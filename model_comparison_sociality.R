# Bayesian Sociality Models: A Scalable and Flexible Alternative for Network Analysis
#
# Section 5.2

# Settings ---------------------------------------------------------------------

## Working directory
setwd("~/Dropbox/PAPERS/projects/sociality")

## Clean global environment
rm(list = ls())

## Load R functions
source("MCMC.R")
source("WAIC.R")
source("r_functions.R")
source("test_functions.R")
source("cv_functions.R")

## MCMC settings
a_sigma <- 2 
b_sigma <- 1/3
a_tau   <- 2 
b_tau   <- 1/3

n_iter <- 250000 + 10000
n_burn <- 10000
n_thin <- 10

#-------------------------------------------------------------------------------

dataset <- "zach"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "bktec"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "foot"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "lazega"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "hitech"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "kaptail"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "bkham"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "dol"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "glossgt"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "lesmis"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "salter"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "polbooks"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "adjnoun"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "football"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "nine"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "gen"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "fblog"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "jazz"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "partner"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "indus"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))

#-------------------------------------------------------------------------------

dataset <- "science"

# load
load(paste0(dataset, "_data.RData"))

# mcmc
set.seed(42)
samples <- gibbs_sampler(Ycube, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)

# test
set.seed(42)
ppps <- sociality_test(Ycube, samples)

# waic
waic <- WAIC(I, B = nrow(samples$delta), Y, mu_chain = samples$mu, delta_chain = samples$delta)

# cv
set.seed(42)
aucs <- sociality_cv(Y, n_iter, n_burn, n_thin)

# results
out <- c(waic$waic1, mean(aucs))

# save
save(samples, ppps, waic, aucs, out, file = paste0(dataset, "_sociality.RData"))