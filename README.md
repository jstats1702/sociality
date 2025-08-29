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

Content:

- Settings  
- Sociality model: Prior simulation 1
- Sociality model: Prior simulation 2

## Section 4.1: Markov chain Monte Carlo

#### Sociality model: MCMC algorithm

File: `MCMC.R`

Required libraries: `truncnorm`

Content: 

- Full conditional distributions: MCMC sociality model
- MCMC: Sociality model for data analysis
- MCMC: Sociality model for cross-validation experiments

## Section 4.2: Variational inference

#### Sociality model: VI algorithm

File: `VI.R`

Required libraries: `truncnorm`

Content:

- Function to initialize variational parameters
- Function to compute the ELBO
- VI algorithm

## Section 5.1: Collaborative working relationships

File: `Lazega_sociality_mcmc.R`

Required libraries: `truncnorm`, `truncnorm`, `igraph`, `sand`, `coda`, `ggplot2`, `reshape2`, `gridExtra`, `cluster`, `mclust`, `truncnorm`

Content:

- Settings
- Data
- Exploratory data analysis
- Model fitting using MCMC
- Convergence diagnostics
- Model fitting using VI
- Inference on $\mu$
- Inference on $\sigma^2$
- Inference on $\tau^2$
- Inference on $\mu$, $\sigma^2$ and $\tau^2$
- Inference on delta using MCMC
- Inference on delta using VI
- Clustering sociality effects using MCMC
- Inference on sociality clusters using MCMC
- Clustering sociality effects using VI
- Inference on sociality clusters using VI
- Interaction probabilities
- Posterior predictive checks
- Test statistics visualization
- Degree predictive check

