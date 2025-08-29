# Lazega

# settings ---------------------------------------------------------------------
# setwd("C:/Users/User/Dropbox/PAPERS/projects/sociality")
setwd("~/Dropbox/PAPERS/projects/sociality")

rm(list = ls())

suppressMessages(suppressWarnings(library(igraph)))
suppressMessages(suppressWarnings(library(sand)))
suppressMessages(suppressWarnings(library(coda)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(gridExtra)))
suppressMessages(suppressWarnings(library(cluster)))
suppressMessages(suppressWarnings(library(mclust)))

source("MCMC.R")
source("VI.R")
source("helper functions.R")
source("r_functions.R")

# data -------------------------------------------------------------------------

g <- graph_from_data_frame(d = elist.lazega, directed = "F")
y <- as.matrix(as_adjacency_matrix(graph = g, names = F))
n <- nrow(y)

# observed test statistics
obs_test_statistics <- c(
     edge_density(g),
     transitivity(g, type = "global"),
     assortativity_degree(g, directed = FALSE),
     mean_distance(g, directed = FALSE),
     mean(degree(g)),
     sd(degree(g))
)

obs_degree <- degree(g)

## adjacency matrix visualization
pdf(file = "fig_lazega_adjacency.pdf", pointsize = 15)
corrplot::corrplot(corr = y, 
                   col.lim = c(0,1), 
                   method = "color", 
                   tl.col = "black", 
                   addgrid.col = "white", 
                   cl.pos = "n", 
                   col = colorRampPalette(c("red4", "white", "black"))(200))
rect(xleft = 0.5, ybottom = 0.5, xright = n + 0.5, ytop = n + 0.5, border = "black", lwd = 1)
dev.off()

## graph visualization
pdf(file = "fig_lazega_graph.pdf", pointsize = 18)
par(mar = c(0, 0, 0, 0))
set.seed(123)
plot(g, 
     layout = layout_with_kk, 
     vertex.label = 1:vcount(g), 
     vertex.label.cex = 0.8,  # Scale down the label size
     vertex.size = 9, 
     vertex.frame.color = "black", 
     vertex.color = 0, 
     vertex.label.color = "black", 
     edge.color = adjustcolor("black", 0.5))
dev.off()

# exploratory data analysis ----------------------------------------------------

## Global Properties
(num_nodes <- vcount(g))
(num_edges <- ecount(g))
(dens <- edge_density(g))
(diam <- diameter(g))
(avg_dist <- mean_distance(g))
(trans <- transitivity(g, type = "global"))
(assor <- assortativity_degree(g))

## Local Properties
(deg <- degree(g))
(avg_degree <- mean(deg))
(degree_dist <- degree_distribution(g))

## Centrality Metrics
(betweenness_centrality <- betweenness(g))
(closeness_centrality   <- closeness(g))
(eigenvector_centrality <- eigen_centrality(g)$vector)

## Community Detection
set.seed(123)
communities <- cluster_fast_greedy(g)
(num_communities <- length(communities))
(modularity_value <- modularity(communities))

## Plot the network with communities
pdf(file = "fig_lazega_graph_communities_fast_greedy.pdf", pointsize = 18)
par(mar = c(0, 0, 0, 0))
community_colors <- membership(communities)
palette <- rainbow(max(community_colors))
set.seed(123)
plot(g, 
     layout = layout_with_kk,
     vertex.size = 9,
     vertex.label = 1:vcount(g),
     vertex.label.cex = 0.8,
     vertex.label.color = "black",
     vertex.color = adjustcolor(palette[community_colors], 0.5),
     vertex.frame.color = adjustcolor(palette[community_colors], 0.5),
     edge.color = adjustcolor("black", 0.5),
     main = "")
dev.off()

# 

# Cluster memberships
clusters <- membership(communities)

## Reorder the matrices based on the reordered clusters
y_reordered <- reorder_matrix(y, clusters)

diag(y_reordered) <- 0

rownames(y_reordered) <- (1:n)[order(clusters)]
colnames(y_reordered) <- (1:n)[order(clusters)]

## adjacency matrix visualization (ordered)
pdf(file = "fig_lazega_adjacency_ordered.pdf", pointsize = 15)
corrplot::corrplot(corr = y_reordered, 
                   col.lim = c(0,1), 
                   method = "color", 
                   tl.col = "black", 
                   addgrid.col = "white", 
                   cl.pos = "n", 
                   col = colorRampPalette(c("red4", "white", "black"))(200))
rect(xleft = 0.5, ybottom = 0.5, xright = n + 0.5, ytop = n + 0.5, border = "black", lwd = 1)
dev.off()

# model fitting using MCMC -----------------------------------------------------

# hyperparameters
a_sigma <- 2 
b_sigma <- 1/3
a_tau   <- 2 
b_tau   <- 1/3

n_iter <- 250000 + 10000
n_burn <- 10000
n_thin <- 10

start.time <- Sys.time()
set.seed(123)
samples <- gibbs_sampler(y, n_iter, n_burn, n_thin, a_sigma, b_sigma, a_tau, b_tau)
end.time <- Sys.time()
save(samples, start.time, end.time, file = "samples_sociality_lazega.RData")

load("samples_sociality_lazega.RData")
end.time - start.time

# convergence ------------------------------------------------------------------

## compute effective sample sizes
round(effectiveSize(samples$mu), 0)
round(effectiveSize(samples$sigma2), 0)
round(effectiveSize(samples$tau2), 0)
round(summary(effectiveSize(samples$delta)), 0)

## compute MC standard errors
round(mcse(samples$mu), 4)
round(mcse(samples$sigma2), 4)
round(mcse(samples$tau2), 4)
round(summary(mcse(samples$delta)), 4)

## compute sampler's log-likelihood in each iteration
# log_lik <- log_likelihood_sociality(y, samples)

## viz log-likelihood trace plot
# log_lik_df <- data.frame(
#      Iteration = seq_along(log_lik),
#      LogLikelihood = log_lik
# )

# ggplot(log_lik_df, aes(x = Iteration, y = LogLikelihood)) +
#      geom_point(size = 0.5, color = "black", alpha = 0.5) +
#      theme_minimal() +
#      labs(
#           title = "", 
#           x = "Iteration", 
#           y = "Log-likelihood"
#      ) +
#      theme(
#           plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
#           axis.title = element_text(size = 12),
#           axis.text  = element_text(size = 10)
#      )

# model fitting using VI ---------------------------------------------------

global_bound <- 3
epsilon <- 1e-06
max_iter <- 1000

start.time <- Sys.time()
variational <- vi_sociality(y, a_sigma, b_sigma, a_tau, b_tau, global_bound, epsilon, max_iter)
end.time <- Sys.time()
end.time - start.time

# inference on mu ------------------------------------------------

## Extract posterior samples for μ
posterior_mu <- samples$mu

## Extract variational approximation parameters for μ
mu_mean_vi <- variational$mu_mu
mu_sd_vi <- sqrt(variational$sigma_mu2)

## Compute kernel density estimate for posterior samples
kde <- density(posterior_mu, from = -1.2, to = -0.6)

## Create a data frame for ggplot
kde_df <- data.frame(x = kde$x, y = kde$y, Type = "MCMC")

## Compute variational density
x_vals <- seq(-1.2, -0.6, length.out = 1000)
vi_density <- dnorm(x_vals, mean = mu_mean_vi, sd = mu_sd_vi)
vi_df <- data.frame(x = x_vals, y = vi_density, Type = "VI")

## Combine both distributions
combined_df <- rbind(kde_df, vi_df)

## Plot using ggplot2
pdf(file = "fig_lazega_posterior_mu.pdf")
ggplot(combined_df, aes(x = x, y = y, color = Type, linetype = Type)) +
     geom_line(size = 0.5) +
     scale_color_manual(values = c("blue", "red")) +
     scale_linetype_manual(values = c("solid", "solid")) +
     xlim(-1.2, -0.6) +
     labs(
          title = "",
          x = expression(mu),
          y = "Density"
     ) +
     theme_minimal(base_size = 30) +
     theme(
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.title = element_blank(),
          legend.background = element_blank()
     )
dev.off()

# inference on sigma^2 ---------------------------------------------------------

## Extract posterior samples for σ²
posterior_sigma2 <- samples$sigma2

## Compute kernel density estimate for posterior samples
kde <- density(posterior_sigma2, from = 0, to = 2)

## Create a data frame for ggplot
kde_df <- data.frame(x = kde$x, y = kde$y, Type = "MCMC")

## Compute variational density (Inverse Gamma approximation)
x_vals <- seq(0, 2, length.out = 1000)
vi_density <- invgamma::dinvgamma(x_vals, shape = variational$alpha_sigma, rate = variational$beta_sigma)
vi_df <- data.frame(x = x_vals, y = vi_density, Type = "VI")

## Combine both distributions
combined_df <- rbind(kde_df, vi_df)

## Plot using ggplot2
pdf(file = "fig_lazega_posterior_sigma2.pdf")
ggplot(combined_df, aes(x = x, y = y, color = Type, linetype = Type)) +
     geom_line(size = 0.5) +
     scale_color_manual(values = c("blue", "red")) +
     scale_linetype_manual(values = c("solid", "solid")) +
     xlim(0, 2) +
     labs(
          title = "",
          x = expression(sigma^2),
          y = "Density"
     ) +
     theme_minimal(base_size = 30) +
     theme(
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.title = element_blank(),
          legend.background = element_blank()
     )
dev.off()
# inference on tau^2 -----------------------------------------------------------

## Extract posterior samples for τ²
posterior_tau2 <- samples$tau2

## Compute kernel density estimate for posterior samples
kde <- density(posterior_tau2, from = 0, to = 0.5)

## Create a data frame for ggplot
kde_df <- data.frame(x = kde$x, y = kde$y, Type = "MCMC")

## Compute variational density (Inverse Gamma approximation)
x_vals <- seq(0, 0.5, length.out = 1000)
vi_density <- invgamma::dinvgamma(x_vals, shape = variational$alpha_tau, rate = variational$beta_tau)
vi_df <- data.frame(x = x_vals, y = vi_density, Type = "VI")

## Combine both distributions
combined_df <- rbind(kde_df, vi_df)

## Plot using ggplot2
pdf(file = "fig_lazega_posterior_tau2.pdf")
ggplot(combined_df, aes(x = x, y = y, color = Type, linetype = Type)) +
     geom_line(size = 0.5) +
     scale_color_manual(values = c("blue", "red")) +
     scale_linetype_manual(values = c("solid", "solid")) +
     xlim(0, 0.5) +
     labs(
          title = "",
          x = expression(tau^2),
          y = "Density"
     ) +
     theme_minimal(base_size = 30) +
     theme(
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.title = element_blank(),
          legend.background = element_blank()
     )
dev.off()

# inference on mu, sigma^2 and tau^2 -------------------------------------------

## MCMC Estimates
mcmc_mu_mean <- mean(samples$mu)
mcmc_mu_sd <- sd(samples$mu)
mcmc_mu_ci <- quantile(samples$mu, probs = c(0.025, 0.975))

mcmc_sigma2_mean <- mean(samples$sigma2)
mcmc_sigma2_sd <- sd(samples$sigma2)
mcmc_sigma2_ci <- quantile(samples$sigma2, probs = c(0.025, 0.975))

mcmc_tau2_mean <- mean(samples$tau2)
mcmc_tau2_sd <- sd(samples$tau2)
mcmc_tau2_ci <- quantile(samples$tau2, probs = c(0.025, 0.975))

## Create MCMC Table
mcmc_table <- data.frame(
     Parameter = c("mu", "sigma^2", "tau^2"),
     Mean = c(mcmc_mu_mean, mcmc_sigma2_mean, mcmc_tau2_mean),
     SD = c(mcmc_mu_sd, mcmc_sigma2_sd, mcmc_tau2_sd),
     CI95_Lower = c(mcmc_mu_ci[1], mcmc_sigma2_ci[1], mcmc_tau2_ci[1]),
     CI95_Upper = c(mcmc_mu_ci[2], mcmc_sigma2_ci[2], mcmc_tau2_ci[2])
)

## VI Estimates
vi_mu_mean <- variational$mu_mu
vi_mu_sd <- sqrt(variational$sigma_mu2)
vi_mu_ci <- qnorm(c(0.025, 0.975), mean = vi_mu_mean, sd = vi_mu_sd)

vi_sigma2_mean <- variational$beta_sigma / (variational$alpha_sigma - 1) # E[X] = beta / (alpha - 1) for alpha > 1
vi_sigma2_sd <- sqrt(variational$beta_sigma^2 / ((variational$alpha_sigma - 1)^2 * (variational$alpha_sigma - 2))) # Var[X] = beta^2 / ((alpha - 1)^2 * (alpha - 2)), for alpha > 2
vi_sigma2_ci <- invgamma::qinvgamma(c(0.025, 0.975), shape = variational$alpha_sigma, rate = variational$beta_sigma)

vi_tau2_mean <- variational$beta_tau / (variational$alpha_tau - 1) # E[X] = beta / (alpha - 1) for alpha > 1
vi_tau2_sd <- sqrt(variational$beta_tau^2 / ((variational$alpha_tau - 1)^2 * (variational$alpha_tau - 2))) # Var[X] = beta^2 / ((alpha - 1)^2 * (alpha - 2)), for alpha > 2
vi_tau2_ci <- invgamma::qinvgamma(c(0.025, 0.975), shape = variational$alpha_tau, rate = variational$beta_tau)

## Create VI Table
vi_table <- data.frame(
     Parameter = c("mu", "sigma^2", "tau^2"),
     Mean = c(vi_mu_mean, vi_sigma2_mean, vi_tau2_mean),
     SD = c(vi_mu_sd, vi_sigma2_sd, vi_tau2_sd),
     CI95_Lower = c(vi_mu_ci[1], vi_sigma2_ci[1], vi_tau2_ci[1]),
     CI95_Upper = c(vi_mu_ci[2], vi_sigma2_ci[2], vi_tau2_ci[2])
)

# Display the Tables
cat("MCMC Table:\n")
print(mcmc_table)
     
cat("\nVI Table:\n")
print(vi_table)

# inference on delta MCMC ------------------------------------------------------

## Create the DataFrame with estimates
delta_mean <- colMeans(samples$delta)  # Posterior mean of delta
delta_ci95 <- apply(samples$delta, 2, quantile, probs = c(0.025, 0.975))

## Create the DataFrame ordered by the posterior mean
delta_df <- data.frame(
     Node = 1:n,  # Original node
     Delta_Est = delta_mean,
     CI95_Lower = delta_ci95[1, ],
     CI95_Upper = delta_ci95[2, ]
)

## Order by the posterior mean
delta_df <- delta_df[order(delta_df$Delta_Est), ]
delta_df$Order <- 1:n 

## Identify intervals that contain 0
delta_df$IntervalType <- ifelse(
     delta_df$CI95_Lower > 0, "Above 0",
     ifelse(delta_df$CI95_Upper < 0, "Below 0", "Contains 0")
)

## Plot using ggplot2
pdf(file = "fig_lazega_posterior_delta_mcmc.pdf")
ggplot(delta_df, aes(x = Order)) +
     # 95% Credible Intervals (thin lines)
     geom_segment(aes(
          x = Order, xend = Order,
          y = CI95_Lower, yend = CI95_Upper,
          color = IntervalType
     ), linewidth = 0.8) +
     
     # Small horizontal lines at the ends
     geom_segment(aes(
          x = Order - 0.2, xend = Order + 0.2, y = CI95_Lower, yend = CI95_Lower,
          color = IntervalType
     ), linewidth = 0.8) +
     
     geom_segment(aes(
          x = Order - 0.2, xend = Order + 0.2, y = CI95_Upper, yend = CI95_Upper,
          color = IntervalType
     ), linewidth = 0.8) +
     
     # Point estimates (colored based on interval type)
     geom_point(aes(y = Delta_Est, color = IntervalType), size = 2) +
     
     # Add numbers above each interval (original node)
     geom_text(aes(y = CI95_Upper + 0.07, label = Node), size = 5, hjust = 0.5, angle = 90) +
     
     # Customization
     scale_color_manual(values = c("Above 0" = "green", "Below 0" = "red", "Contains 0" = "gray70")) +
     labs(
          title = "",
          x = NULL,
          y = expression(delta)
     ) +
     ylim(-1.25, 1.25) +
     theme_minimal(base_size = 20) +
     theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
     )
dev.off()

# inference on delta VI --------------------------------------------------------

## Extract VI estimates for delta
delta_mean_vi <- variational$mu_delta  # Variational posterior mean of delta
delta_sd_vi <- sqrt(variational$sigma_delta2)  # Variational posterior standard deviations
delta_ci95_vi <- t(apply(cbind(delta_mean_vi, delta_sd_vi), 1, function(x) {
     c(Lower = x[1] - 1.96 * x[2], Upper = x[1] + 1.96 * x[2])
}))

## Create the DataFrame with VI estimates
delta_df_vi <- data.frame(
     Node = 1:length(delta_mean_vi),  # Original node
     Delta_Est = delta_mean_vi,
     CI95_Lower = delta_ci95_vi[, 1],
     CI95_Upper = delta_ci95_vi[, 2]
)

## Order by the posterior mean
delta_df_vi <- delta_df_vi[order(delta_df_vi$Delta_Est), ]
delta_df_vi$Order <- 1:length(delta_mean_vi)

## Identify intervals that contain 0
delta_df_vi$IntervalType <- ifelse(
     delta_df_vi$CI95_Lower > 0, "Above 0",
     ifelse(delta_df_vi$CI95_Upper < 0, "Below 0", "Contains 0")
)

## Plot using ggplot2
pdf(file = "fig_lazega_posterior_delta_vi.pdf")
ggplot(delta_df_vi, aes(x = Order)) +
     # 95% Credible Intervals (thin lines)
     geom_segment(aes(
          x = Order, xend = Order,
          y = CI95_Lower, yend = CI95_Upper,
          color = IntervalType
     ), linewidth = 0.8) +
     
     # Small horizontal lines at the ends
     geom_segment(aes(
          x = Order - 0.2, xend = Order + 0.2, y = CI95_Lower, yend = CI95_Lower,
          color = IntervalType
     ), linewidth = 0.8) +
     
     geom_segment(aes(
          x = Order - 0.2, xend = Order + 0.2, y = CI95_Upper, yend = CI95_Upper,
          color = IntervalType
     ), linewidth = 0.8) +
     
     # Point estimates (colored based on interval type)
     geom_point(aes(y = Delta_Est, color = IntervalType), size = 2) +
     
     # Add numbers above each interval (original node)
     geom_text(aes(y = CI95_Upper + 0.07, label = Node), size = 5, hjust = 0.5, angle = 90) +
     
     # Customization
     scale_color_manual(values = c("Above 0" = "green", "Below 0" = "red", "Contains 0" = "gray70")) +
     labs(
          title = "",
          x = NULL,
          y = expression(delta)
     ) +
     ylim(-1.25, 1.25) +
     theme_minimal(base_size = 14) +
     theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
     )
dev.off()

# clustering MCMC --------------------------------------------------------------

## Calculate co-clustering probabilities
coclustering_probs <- compute_coclustering(samples)

## Estimate clusters using Mclust
clusters <- estimate_clusters_mclust(coclustering_probs)

## Create the binary matrix E
E <- compute_binary_E(clusters)

## Reorder the matrices based on the reordered clusters
coclustering_reordered <- reorder_matrix(coclustering_probs, clusters)
E_reordered <- reorder_matrix(E, clusters)

diag(coclustering_reordered) <- 1
diag(E_reordered) <- 1

rownames(coclustering_reordered) <- (1:n)[order(clusters)]
colnames(coclustering_reordered) <- (1:n)[order(clusters)]

rownames(E_reordered) <- (1:n)[order(clusters)]
colnames(E_reordered) <- (1:n)[order(clusters)]

## Visualization of coclustering probabilities
pdf(file = "fig_lazega_coclustering_probabilities.pdf", pointsize = 15)
corrplot::corrplot(corr = coclustering_reordered,
                   col.lim = c(0, 1), 
                   method = "color", 
                   tl.col = "black", 
                   addgrid.col = "white", 
                   cl.pos = "n", 
                   col = colorRampPalette(c("red", "white", "red"))(200))

## Add border around the matrix
rect(xleft = 0.5, ybottom = 0.5, xright = n + 0.5, ytop = n + 0.5, border = "black", lwd = 1)
dev.off()

## Visualization of partition
pdf(file = "fig_lazega_partition_mcmc.pdf", pointsize = 15)
corrplot::corrplot(corr = E_reordered,
                   col.lim = c(0, 1), 
                   method = "color", 
                   tl.col = "black", 
                   addgrid.col = "white", 
                   cl.pos = "n", 
                   col = colorRampPalette(c("red", "white", "red"))(200))

## Add border around the matrix
rect(xleft = 0.5, ybottom = 0.5, xright = n + 0.5, ytop = n + 0.5, border = "black", lwd = 1)
dev.off()

# inference on delta clusters MCMC --------------------------------- 

clusters_mcmc <- clusters

## Create the DataFrame with estimates
delta_mean <- colMeans(samples$delta)  # Posterior mean of delta
delta_ci95 <- apply(samples$delta, 2, quantile, probs = c(0.025, 0.975))

## Create the DataFrame with cluster assignments
delta_df <- data.frame(
     Node = 1:n,  # Original node
     Delta_Est = delta_mean,
     CI95_Lower = delta_ci95[1, ],
     CI95_Upper = delta_ci95[2, ],
     Cluster = clusters_mcmc  # Add cluster assignments
)

## Order by the posterior mean
delta_df <- delta_df[order(delta_df$Delta_Est), ]
delta_df$Order <- 1:n 

## Assign unique colors to clusters
cluster_colors <- scales::hue_pal()(length(unique(delta_df$Cluster)))

## Plot using ggplot2
pdf(file = "fig_lazega_posterior_delta_mcmc_clustering.pdf")
ggplot(delta_df, aes(x = Order)) +
     # 95% Credible Intervals (thin lines)
     geom_segment(aes(
          x = Order, xend = Order,
          y = CI95_Lower, yend = CI95_Upper,
          color = as.factor(Cluster)  # Color by cluster
     ), linewidth = 0.8) +
     
     # Small horizontal lines at the ends
     geom_segment(aes(
          x = Order - 0.2, xend = Order + 0.2, y = CI95_Lower, yend = CI95_Lower,
          color = as.factor(Cluster)  # Color by cluster
     ), linewidth = 0.8) +
     
     geom_segment(aes(
          x = Order - 0.2, xend = Order + 0.2, y = CI95_Upper, yend = CI95_Upper,
          color = as.factor(Cluster)  # Color by cluster
     ), linewidth = 0.8) +
     
     # Point estimates (colored by cluster)
     geom_point(aes(y = Delta_Est, color = as.factor(Cluster)), size = 2) +
     
     # Add numbers above each interval (original node)
     geom_text(aes(y = CI95_Upper + 0.07, label = Node), size = 5, hjust = 0.5, angle = 90) +
     
     # Customization
     scale_color_manual(values = cluster_colors) +
     labs(
          title = "",
          x = NULL,
          y = expression(delta)
     ) +
     ylim(-1.25, 1.25) +
     theme_minimal(base_size = 20) +
     theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
     )
dev.off()

# visualizing graph clusters MCMC ----------------------------------------------

## Map the same cluster colors from the previous plot
cluster_colors <- scales::hue_pal()(length(unique(clusters_mcmc)))

## Assign colors to each node based on clusters_mcmc
vertex_colors <- cluster_colors[clusters_mcmc]

## Plot the network with communities
pdf(file = "fig_lazega_graph_communities_mcmc.pdf", pointsize = 18)
par(mar = c(0, 0, 0, 0))  # Remove margins
set.seed(123)  # Ensure layout consistency
plot(g, 
     layout = layout_with_kk,  # Kamada-Kawai layout
     vertex.size = 9, 
     vertex.label = 1:vcount(g),  # Node labels
     vertex.label.cex = 0.8,  # Label size
     vertex.label.color = "black",  # Label color
     vertex.color = adjustcolor(vertex_colors, 0.5),  # Node fill color based on cluster
     vertex.frame.color = adjustcolor(vertex_colors, 0.5),  # Node border color matches fill
     edge.color = adjustcolor("black", 0.5),  # Edge color
     main = "")  # No title
dev.off()

# clustering VI ----------------------------------------------------------------

## Extract the values to be clustered
mu_delta <- variational$mu_delta

## Determine the optimal number of clusters using the elbow method
wss <- numeric(10)
for (k in 1:10) {
     kmeans_result <- kmeans(mu_delta, centers = k, nstart = 25)
     wss[k] <- kmeans_result$tot.withinss
}

## Choose the optimal number of clusters
optimal_k <- which(diff(diff(wss)) == min(diff(diff(wss)))) + 1

## Perform k-means clustering with the optimal number of clusters
kmeans_final <- kmeans(mu_delta, centers = optimal_k, nstart = 25)

# Cluster assignments
clusters <- kmeans_final$cluster

# inference on delta clusters VI ----------------------------------------------- 

clusters_vi <- clusters

## Extract VI estimates for delta
delta_mean_vi <- variational$mu_delta  # Variational posterior mean of delta
delta_sd_vi <- sqrt(variational$sigma_delta2)  # Variational posterior standard deviations
delta_ci95_vi <- t(apply(cbind(delta_mean_vi, delta_sd_vi), 1, function(x) {
     c(Lower = x[1] - 1.96 * x[2], Upper = x[1] + 1.96 * x[2])
}))

## Create the DataFrame with VI estimates
delta_df_vi <- data.frame(
     Node = 1:length(delta_mean_vi),  # Original node
     Delta_Est = delta_mean_vi,
     CI95_Lower = delta_ci95_vi[, 1],
     CI95_Upper = delta_ci95_vi[, 2],
     Cluster = clusters_vi  # Add cluster assignments
)

## Order by the posterior mean
delta_df_vi <- delta_df_vi[order(delta_df_vi$Delta_Est), ]
delta_df_vi$Order <- 1:length(delta_mean_vi)

## Assign unique colors to clusters
cluster_colors <- scales::hue_pal()(length(unique(delta_df_vi$Cluster)))

## Plot using ggplot2
pdf(file = "fig_lazega_posterior_delta_vi_clustering.pdf")
ggplot(delta_df_vi, aes(x = Order)) +
     # 95% Credible Intervals (thin lines)
     geom_segment(aes(
          x = Order, xend = Order,
          y = CI95_Lower, yend = CI95_Upper,
          color = as.factor(Cluster)  # Color by cluster
     ), linewidth = 0.8) +
     
     # Small horizontal lines at the ends
     geom_segment(aes(
          x = Order - 0.2, xend = Order + 0.2, y = CI95_Lower, yend = CI95_Lower,
          color = as.factor(Cluster)  # Color by cluster
     ), linewidth = 0.8) +
     
     geom_segment(aes(
          x = Order - 0.2, xend = Order + 0.2, y = CI95_Upper, yend = CI95_Upper,
          color = as.factor(Cluster)  # Color by cluster
     ), linewidth = 0.8) +
     
     # Point estimates (colored by cluster)
     geom_point(aes(y = Delta_Est, color = as.factor(Cluster)), size = 2) +
     
     # Add numbers above each interval (original node)
     geom_text(aes(y = CI95_Upper + 0.07, label = Node), size = 5, hjust = 0.5, angle = 90) +
     
     # Customization
     scale_color_manual(values = cluster_colors) +
     labs(
          title = "",
          x = NULL,
          y = expression(delta)
     ) +
     ylim(-1.25, 1.25) +
     theme_minimal(base_size = 14) +
     theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
     )
dev.off()

# visualizing graph clusters VI ------------------------------------------------

## Map the same cluster colors from the previous plot
cluster_colors <- scales::hue_pal()(length(unique(clusters_vi)))

## Assign colors to each node based on clusters_mcmc
vertex_colors <- cluster_colors[clusters_vi]

## Plot the network with communities
pdf(file = "fig_lazega_graph_communities_vi.pdf", pointsize = 18)
par(mar = c(0, 0, 0, 0))  # Remove margins
set.seed(123)  # Ensure layout consistency
plot(g, 
     layout = layout_with_kk,  # Kamada-Kawai layout
     vertex.size = 9, 
     vertex.label = 1:vcount(g),  # Node labels
     vertex.label.cex = 0.8,  # Label size
     vertex.label.color = "black",  # Label color
     vertex.color = adjustcolor(vertex_colors, 0.5),  # Node fill color based on cluster
     vertex.frame.color = adjustcolor(vertex_colors, 0.5),  # Node border color matches fill
     edge.color = adjustcolor("black", 0.5),  # Edge color
     main = "")  # No title
dev.off()

# interaction probabilities ----------------------------------------------------

# Number of posterior draws and number of nodes
B <- nrow(samples$delta)
n <- ncol(samples$delta)

# Initialize interaction probabilities matrix
Theta <- matrix(0, nrow = n, ncol = n)

# Loop over posterior draws
for (b in seq_len(B)) {
     
     # Extract current draw
     mu    <- samples$mu[b]
     delta <- samples$delta[b, ]
     
     # Fill only the upper triangle (i < j) to avoid double work
     for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
               Theta[i, j] <- Theta[i, j] + pnorm(mu + delta[i] + delta[j])
          }
     }
}

# Average over posterior draws
Theta <- Theta / B

# Reflect to make the matrix symmetric (since we only filled i < j)
Theta <- Theta + t(Theta)

# No self-ties: set diagonal to NA (use 0 instead if preferred)
diag(Theta) <- 0

# 

# Cluster memberships
clusters <- membership(communities)

## Reorder the matrices based on the reordered clusters
Theta_reordered <- reorder_matrix(Theta, clusters)

diag(Theta_reordered) <- 0

rownames(Theta_reordered) <- (1:n)[order(clusters)]
colnames(Theta_reordered) <- (1:n)[order(clusters)]

## Visualization of coclustering probabilities (ordered)
pdf(file = "fig_lazega_interaction_probabilities_ordered.pdf", pointsize = 15)
corrplot::corrplot(corr = Theta_reordered,
                   col.lim = c(0, 1), 
                   method = "color", 
                   tl.col = "black", 
                   addgrid.col = "white", 
                   cl.pos = "n", 
                   col = colorRampPalette(c("red", "white", "red"))(200))

## Add border around the matrix
rect(xleft = 0.5, ybottom = 0.5, xright = n + 0.5, ytop = n + 0.5, border = "black", lwd = 1)
dev.off()

# test statistics --------------------------------------------------------------

# T_1 : density
# T_2 : transitivity
# T_3 : degree assortativity
# T_4 : mean geodecic distance
# T_5 : mean degree
# T_6 : sd degree 

## settings
B <- nrow(samples$delta)  # Get number of samples
n_test_stats <- 6
test_stats_sociality <- matrix(NA, B, n_test_stats)  # Preallocate matrix

set.seed(42)
for (b in seq_len(B)) {
     # Extract parameters
     mu <- samples$mu[b]
     delta <- samples$delta[b, ]

     # Simulate adjacency matrix and convert to graph
     A <- matrix(0, n, n)
     for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
               A[i, j] <- rbinom(1, size = 1, prob = pnorm(mu + delta[i] + delta[j]))
               A[j, i] <- A[i, j]
          }
     }
     g <- graph_from_adjacency_matrix(A, mode = "undirected")

     # Compute node degrees once
     degrees <- igraph::degree(g)

     # Compute statistics efficiently
     test_stats_sociality[b, ] <- c(
          igraph::edge_density(g),
          igraph::transitivity(g, type = "global"),
          igraph::assortativity_degree(g),
          igraph::mean_distance(g),
          mean(degrees),
          sd(degrees)
     )

     # Display progress every 10%
     if (b %% ceiling(0.1 * B) == 0 || b == B) {
          cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
     }
}

save(test_stats_sociality, obs_test_statistics, file = "test_statistics_sociality_lazega.RData")

# load tests statistics data
load("test_statistics_sociality_lazega.RData")
load("test_statistics_distance_lazega.RData")
load("test_statistics_class_lazega.RData")
load("test_statistics_eigen_lazega.RData")

# test statistics visualization ------------------------------------------------

pdf(file = "fig_lazega_test_statistics_density.pdf", pointsize = 22)
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
test_statistic_viz(1, "Density", obs_test_statistics, test_stats_sociality, 
                   test_stats_distance, test_stats_class, test_stats_eigen)
dev.off()

pdf(file = "fig_lazega_test_statistics_transitivity.pdf", pointsize = 22)
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
test_statistic_viz(2, "Transitivity", obs_test_statistics, test_stats_sociality, 
                   test_stats_distance, test_stats_class, test_stats_eigen)
dev.off()

pdf(file = "fig_lazega_test_statistics_assortativity.pdf", pointsize = 22)
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
test_statistic_viz(3, "Assortativity", obs_test_statistics, test_stats_sociality, 
                   test_stats_distance, test_stats_class, test_stats_eigen)
dev.off()

pdf(file = "fig_lazega_test_statistics_mean_distance.pdf", pointsize = 22)
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
test_statistic_viz(4, "Mean distance", obs_test_statistics, test_stats_sociality, 
                   test_stats_distance, test_stats_class, test_stats_eigen)
dev.off()

pdf(file = "fig_lazega_test_statistics_mean_degree.pdf", pointsize = 22)
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
test_statistic_viz(5, "Mean degree", obs_test_statistics, test_stats_sociality, 
                   test_stats_distance, test_stats_class, test_stats_eigen)
dev.off()

pdf(file = "fig_lazega_test_statistics_sd_degree.pdf", pointsize = 22)
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
test_statistic_viz(6, "SD degree", obs_test_statistics, test_stats_sociality, 
                   test_stats_distance, test_stats_class, test_stats_eigen)
dev.off()

# degree predictive check ------------------------------------------------------

## settings
B <- nrow(samples$delta)  # Get number of samples
test_degree_sociality <- matrix(NA, B, n)  # Preallocate matrix

set.seed(42)
for (b in seq_len(B)) {
     # Extract parameters
     mu <- samples$mu[b]
     delta <- samples$delta[b, ]
     
     # Simulate adjacency matrix and convert to graph
     A <- matrix(0, n, n)
     for (i in 1:(n - 1)) {
          for (j in (i + 1):n) {
               A[i, j] <- rbinom(1, size = 1, prob = pnorm(mu + delta[i] + delta[j]))
               A[j, i] <- A[i, j]
          }
     }
     g <- graph_from_adjacency_matrix(A, mode = "undirected")
     
     # Compute statistics efficiently
     test_degree_sociality[b, ] <- igraph::degree(g)
     
     # Display progress every 10%
     if (b %% ceiling(0.1 * B) == 0 || b == B) {
          cat("Simulation progress:", round(100 * b / B, 1), "% completed\n")
     }
}

save(test_degree_sociality, obs_degree, file = "test_degree_sociality_lazega.RData")

## 

model <- c("sociality", "distance", "class", "eigen")

ord <- order(obs_degree, decreasing = F)
obs_ord <- obs_degree[ord]

## x positions and labels (original node indices)
x <- seq_len(n)
x_labels <- ord

## Plot
pdf(file = "fig_lazega_degree_statistics.pdf", pointsize = 14)
par(mar = c(2.75,2.75,0.5,0.5), mgp = c(1.7,0.7,0))
plot(x, mean_ord, type = "n",
     xlab = "Node",
     ylab = "Degree",
     ylim = c(0, 20),
     xaxt = "n")  # custom x-axis below

axis(side = 1, at = x, labels = ord, cex.axis = 0.9, las = 2)

for (j in 1:4) {

     load(paste0("test_degree_", model[j], "_lazega.RData"))
     
     test_degree <- get(paste0("test_degree_", model[j]))
     
     ## Posterior predictive summaries for degrees by node
     post_mean <- colMeans(test_degree, na.rm = TRUE)
     ci95_lo <- apply(X = test_degree, MARGIN = 2, FUN = quantile, prob = 0.025)
     ci95_hi <- apply(X = test_degree, MARGIN = 2, FUN = quantile, prob = 0.975)
     
     ## Reorder all per-node vectors
     mean_ord    <- post_mean [ord]
     ci95_lo_ord <- ci95_lo   [ord]
     ci95_hi_ord <- ci95_hi   [ord]
     
     ## posterior mean and 95% CI
     lines(x, ci95_lo_ord, type = "l", col = adjustcolor(j, 0.5))
     lines(x, ci95_hi_ord, type = "l", col = adjustcolor(j, 0.5))
     lines(x, mean_ord, type = "l", col = adjustcolor(j, 0.5))
}

lines(x, obs_ord, type = "l", col = "gray", lwd = 8)

legend("topleft", legend = model, fill = 1:4, border = 1:4, bty = "n")

dev.off()


