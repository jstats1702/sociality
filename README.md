# Bayesian Sociality Models: A Scalable and Flexible Alternative for Network Analysis

This repository contains all the necessary code to reproduce the results presented in the paper *Bayesian Sociality Models: A Scalable and Flexible Alternative for Network Analysis*.

## Implementation notes

- The implementation is carried out in R software version 4.4.2, using RStudio version 2024.12.0 Build 467 as the IDE, and executed on a standard laptop equipped with an 11th Gen Intel(R) Core(TM) i7-1165G7 processor (2.80 GHz) and 8.00 GB of RAM.

- The user should expect some variation in computation time, particularly when running on different hardware configurations.

- To ensure the software runs properly, the user should place all files in a single folder and update the file path specified at the beginning of each script.

- The user must have the following libraries installed in the global environment to avoid any potential runtime conflicts.

## Section 3.2: Prior elicitation

#### Prior simulation for the sociality model

File: `prior simulation.R`

Required libraries: `ggplot2`, `dplyr`

## Section 4.1: Markov chain Monte Carlo

#### Sociality model: MCMC algorithm

File: `MCMC.R`

Required libraries: `truncnorm`

## Section 4.2: Variational inference

#### Sociality model: VI algorithm

File: `VI.R`

Required libraries: `truncnorm`

## Section 5.1: Collaborative working relationships

#### Sociality model: Exploratory data analysis, model fitting, inference, posterior predictive checks

File: `Lazega_sociality_mcmc_vi.R`

Required libraries: `truncnorm`, `truncnorm`, `igraph`, `sand`, `coda`, `corrplot`, `ggplot2`, `reshape2`, `gridExtra`, `cluster`, `mclust`, `truncnorm`

Depends on: `MCMC.R`, `VI.R`, `helper_functions.R`, `r_functions.R`

#### Sociality model: Sensitivity analysis

File: `Lazega_sociality_mcmc_sensitivity.R`

Required libraries: `truncnorm`,`igraph`, `corrplot`

Depends on: `MCMC.R`, `VI.R`, `helper_functions.R`, `r_functions.R`

#### Sociality model: Model fit (WAIC)

File: `Lazega_sociality_waic.R`

Required libraries: `igraph`, `sand`

Depends on: `WAIC.R`, `r_functions.R`

#### Sociality model: Predictive accuracy (CV)

File: `Lazega_sociality_cv.R`

Required libraries: `igraph`, `sand`, `doParallel`, `ROCR`

Depends on: `MCMC.R`, `r_functions.R`

#### Class model: Model fitting, posterior predictive checks

File: `Lazega_class.R`

Required libraries: `Rcpp`, `igraph`, `sand`

Depends on: `class_functions.cpp`, `r_functions.R`, `class_functions.R`

#### Class model: Model fit (WAIC)

File: `Lazega_class_waic.R`

Required libraries: `Rcpp`, `igraph`, `sand`

Depends on: `class_functions.cpp`, `r_functions.R`, `class_functions.R`

#### Class model: Predictive accuracy (CV)

File: `Lazega_class_cv.R`

Required libraries: `igraph`, `sand`, `doParallel`, `ROCR`, `Rcpp`

Depends on: `class_functions.cpp`, `r_functions.R`, `class_functions.R`

#### Distance model: Model fitting, posterior predictive checks

File: `Lazega_distance.R`

Required libraries: `Rcpp`, `igraph`, `sand`

Depends on: `distance_functions.cpp`, `r_functions.R`, `distance_functions.R`

#### Distance model: Model fit (WAIC)

File: `Lazega_distance_waic.R`

Required libraries: `Rcpp`, `igraph`, `sand`

Depends on: `distance_functions.cpp`, `r_functions.R`, `distance_functions.R`

#### Distance model: Predictive accuracy (CV)

File: `Lazega_distance_cv.R`

Required libraries: `igraph`, `sand`, `doParallel`, `ROCR`, `Rcpp`

Depends on: `distance_functions.cpp`, `r_functions.R`, `distance_functions.R`

#### Eigen model: Model fitting, posterior predictive checks

File: `Lazega_eigen.R`

Required libraries: `Rcpp`, `igraph`, `sand`

Depends on: `eigen_functions.cpp`, `r_functions.R`, `eigen_functions.R`

#### Eigen model: Model fit (WAIC)

File: `Lazega_eigen_waic.R`

Required libraries: `Rcpp`, `igraph`, `sand`

Depends on: `eigen_functions.cpp`, `r_functions.R`, `eigen_functions.R`

#### Eigen model: Predictive accuracy (CV)

File: `Lazega_eigen_cv.R`

Required libraries: `igraph`, `sand`, `doParallel`, `ROCR`, `Rcpp`

Depends on: `eigen_functions.cpp`, `r_functions.R`, `eigen_functions.R`
