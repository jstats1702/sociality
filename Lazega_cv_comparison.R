# Bayesian Sociality Models: A Scalable and Flexible Alternative for Network Analysis
#
# Section 5.1

# Settings ---------------------------------------------------------------------

## Working directory
setwd("~/Dropbox/PAPERS/projects/sociality")

## Clean global environment
rm(list = ls())

## Required libraries
library(ROCR)

# CV comparison ----------------------------------------------------------------

## Run Lazega_sociality_cv.R
## Run Lazega_class_cv.R
## Run Lazega_distance_cv.R
## Run Lazega_eigen_cv.R

## Number of folds
L <- 5

## Models
model <- c("sociality", "distance", "class", "eigen")

model_legend <- c("sociality ", 
                  "distance", 
                  "class     ", 
                  "eigen    ")

#---------------#
#   Figure 11   #
#---------------#

## Plots
for (l in 1:L) {

     # Set up ROC plot
     pdf(file = paste0("fig_lazega_cv_comparison_fold_", l, ".pdf"), pointsize = 25)
     par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
     
     plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), 
          xlab = "False Positive Rate", ylab = "True Positive Rate", main = "")
     abline(a = 0, b = 1, col = "gray", lwd = 2)
     
     # Store AUCs
     aucs <- numeric(4)
     
     for (k in 1:4) {
          # Load cv
          load(paste0("cv_", model[k],"_lazega.RData"))
          
          # Compute performance metrics
          pred <- prediction(predictions = cv[[l]]$y_ppp, labels = cv[[l]]$y_true)
          perf_roc <- performance(pred, measure = "tpr", x.measure = "fpr")
               
          # Plot ROC curve
          lines(perf_roc@x.values[[1]], perf_roc@y.values[[1]], type = "l", 
                lwd = 5, col = adjustcolor(k, 0.8))
               
          # Compute AUC
          perf_auc <- performance(pred, measure = "auc")
          aucs[k] <- perf_auc@y.values[[1]]
     }
     
     legend("bottomright", fill = 1:4, border = 1:4,
            legend = paste0(model_legend, " = ", round(aucs, 3)), 
            bty = "n")
     
     dev.off()
}

# End --------------------------------------------------------------------------