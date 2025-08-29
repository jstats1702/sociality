# Bayesian Sociality Models: A Scalable and Flexible Alternative for Network Analysis
#
# Section 3.2

# Settings ---------------------------------------------------------------------

# Clean global environment
rm(list = ls())

# Working directory
setwd("~/Dropbox/PAPERS/projects/sociality")

# Load required libraries
library(ggplot2)
library(dplyr)

# Sociality model: Prior simulation 1 ------------------------------------------

## Generate data
set.seed(123)
mu <- rnorm(100000, 0, sqrt(1/3))
deltai <- rnorm(100000, 0, sqrt(1/3))
deltaj <- rnorm(100000, 0, sqrt(1/3))
x <- mu + deltai + deltaj
pnorm_x <- pnorm(x)

## Create a data frame for ggplot2
data <- data.frame(pnorm_x = pnorm_x)

#-----------------------#
#   Figure 2, panel a   #
#-----------------------#

## Plot the histogram using ggplot2
pdf(file = "fig_prior_simulation_1.pdf")
ggplot(data, aes(x = pnorm_x)) +
     geom_histogram(aes(y = ..density..), bins = 20, fill = "gray", color = "gray", alpha = 0.9) +
     labs(
          title = "",
          x = "x",
          y = "Density"
     ) +
     theme_minimal(base_size = 25) +
     theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold")
     )
dev.off()

# Sociality model: Prior simulation 1 ------------------------------------------

# Generate data
set.seed(123)
mu <- rnorm(100000, 0, sqrt(1))
deltai <- rnorm(100000, 0, sqrt(1))
deltaj <- rnorm(100000, 0, sqrt(1))
x <- mu + deltai + deltaj
pnorm_x <- pnorm(x)

## Create a data frame for ggplot2
data <- data.frame(pnorm_x = pnorm_x)

#-----------------------#
#   Figure 2, panel b   #
#-----------------------#

## Plot the histogram using ggplot2
pdf(file = "fig_prior_simulation_2.pdf")
ggplot(data, aes(x = pnorm_x)) +
     geom_histogram(aes(y = ..density..), bins = 20, fill = "gray", color = "gray", alpha = 0.9) +
     labs(
          title = "",
          x = "x",
          y = "Density"
     ) +
     theme_minimal(base_size = 25) +
     theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title = element_text(face = "bold")
     )
dev.off()

# End --------------------------------------------------------------------------