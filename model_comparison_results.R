# Bayesian Sociality Models: A Scalable and Flexible Alternative for Network Analysis
#
# Section 5.2

## Working directory
setwd("~/Dropbox/PAPERS/projects/sociality")

## Clean global environment
rm(list = ls())

## Run model comparison_sociality.R
## Run model comparison_class.R
## Run model comparison_distance.R
## Run model comparison_eigen.R

#-------------#
#   Table 3   #
#-------------#

datasets <- c("zach", "bktec", "foot", "hitech", "kaptail", "bkham",
              "dol", "glossgt", "lesmis", "salter", "polbooks", "adjnoun",
              "football", "nine", "gen", "fblog", "jazz", "partner", "indus",
              "science")

results <- NULL
for (i in 1:length(datasets)) {
     tmp <- NULL
     load(file = paste0(datasets[i], "_sociality.RData"))
     tmp <- c(tmp, out)
     load(file = paste0(datasets[i], "_distance.RData"))
     tmp <- c(tmp, out)
     load(file = paste0(datasets[i], "_class.RData"))
     tmp <- c(tmp, out)
     load(file = paste0(datasets[i], "_eigen.RData"))
     tmp <- c(tmp, out)
     results <- rbind(results, tmp)
}

results_waic <- round(results[ , c(1,3,5,7)], 3)
results_auc  <- round(results[ , c(2,4,6,8)], 3)

colnames(results_waic) <- c("soci", "dist", "class", "eigen")
colnames(results_auc)  <- c("soci", "dist", "class", "eigen")
rownames(results_waic) <- datasets
rownames(results_auc)  <- datasets

round(results_waic, 1)
round(results_auc,  3)
