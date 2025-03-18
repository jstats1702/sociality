# Clear the global environment
rm(list = ls(all.names = TRUE))

# Set working directory
setwd("C:/Users/User/Dropbox/PAPERS/projects/sociality")
# setwd("C:/Users/Juan Camilo/Dropbox/PAPERS/projects/sociality")

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

results1 <- NULL  # density
results2 <- NULL  # transitivity
results3 <- NULL  # assortativity
for (i in 1:length(datasets)) {
     tmp1 <- NULL
     tmp2 <- NULL
     tmp3 <- NULL
     load(file = paste0(datasets[i], "_sociality.RData"))
     tmp1 <- c(tmp1, ppps[1])
     tmp2 <- c(tmp2, ppps[2])
     tmp3 <- c(tmp3, ppps[3])
     load(file = paste0(datasets[i], "_distance.RData"))
     tmp1 <- c(tmp1, ppps[1])
     tmp2 <- c(tmp2, ppps[2])
     tmp3 <- c(tmp3, ppps[3])
     load(file = paste0(datasets[i], "_class.RData"))
     tmp1 <- c(tmp1, ppps[1])
     tmp2 <- c(tmp2, ppps[2])
     tmp3 <- c(tmp3, ppps[3])
     load(file = paste0(datasets[i], "_eigen.RData"))
     tmp1 <- c(tmp1, ppps[1])
     tmp2 <- c(tmp2, ppps[2])
     tmp3 <- c(tmp3, ppps[3])
     results1 <- rbind(results1, tmp1)
     results2 <- rbind(results2, tmp2)
     results3 <- rbind(results3, tmp3)
}

colnames(results1) <- c("soci", "dist", "class", "eigen")
colnames(results2) <- c("soci", "dist", "class", "eigen")
colnames(results3) <- c("soci", "dist", "class", "eigen")
rownames(results1) <- datasets
rownames(results2) <- datasets
rownames(results3) <- datasets

round(results1, 3)
round(results2, 3)
round(results3, 3)

# density
par(mfrow = c(1,1), mar = c(3,3,1.5,1.5), mgp = c(1.75,0.75,0))
plot(NA, NA, xlim = c(0.5,4.5), ylim = c(0, 1), xlab = "Model", ylab = "ppp", xaxt = "n")
axis(side = 1, at = 1:4, labels = c("soc", "dist", "class", "eigen"))
abline(h = 0.5, col = "darkgrey", lty = 2, lwd = 2)
for (i in 1:length(datasets)) {
     load(file = paste0(datasets[i], "_sociality.RData"))
     points(x = 1, y = ppps[1], pch = 16, cex = 1.5, col = adjustcolor(1, 0.5))
     load(file = paste0(datasets[i], "_distance.RData"))
     points(x = 2, y = ppps[1], pch = 16, cex = 1.5, col = adjustcolor(2, 0.5))
     load(file = paste0(datasets[i], "_class.RData"))
     points(x = 3, y = ppps[1], pch = 16, cex = 1.5, col = adjustcolor(3, 0.5))
     load(file = paste0(datasets[i], "_eigen.RData"))
     points(x = 4, y = ppps[1], pch = 16, cex = 1.5, col = adjustcolor(4, 0.5))
}

# transitivity
par(mfrow = c(1,1), mar = c(3,3,1.5,1.5), mgp = c(1.75,0.75,0))
plot(NA, NA, xlim = c(0.5,4.5), ylim = c(0, 1), xlab = "Model", ylab = "ppp", xaxt = "n")
axis(side = 1, at = 1:4, labels = c("soc", "dist", "class", "eigen"))
abline(h = 0.5, col = "darkgrey", lty = 2, lwd = 2)
for (i in 1:length(datasets)) {
     load(file = paste0(datasets[i], "_sociality.RData"))
     points(x = 1, y = ppps[2], pch = 16, cex = 1.5, col = adjustcolor(1, 0.5))
     load(file = paste0(datasets[i], "_distance.RData"))
     points(x = 2, y = ppps[2], pch = 16, cex = 1.5, col = adjustcolor(2, 0.5))
     load(file = paste0(datasets[i], "_class.RData"))
     points(x = 3, y = ppps[2], pch = 16, cex = 1.5, col = adjustcolor(3, 0.5))
     load(file = paste0(datasets[i], "_eigen.RData"))
     points(x = 4, y = ppps[2], pch = 16, cex = 1.5, col = adjustcolor(4, 0.5))
}

# assortativity
par(mfrow = c(1,1), mar = c(3,3,1.5,1.5), mgp = c(1.75,0.75,0))
plot(NA, NA, xlim = c(0.5,4.5), ylim = c(0, 1), xlab = "Model", ylab = "ppp", xaxt = "n")
axis(side = 1, at = 1:4, labels = c("soc", "dist", "class", "eigen"))
abline(h = 0.5, col = "darkgrey", lty = 2, lwd = 2)
for (i in 1:length(datasets)) {
     load(file = paste0(datasets[i], "_sociality.RData"))
     points(x = 1, y = ppps[3], pch = 16, cex = 1.5, col = adjustcolor(1, 0.5))
     load(file = paste0(datasets[i], "_distance.RData"))
     points(x = 2, y = ppps[3], pch = 16, cex = 1.5, col = adjustcolor(2, 0.5))
     load(file = paste0(datasets[i], "_class.RData"))
     points(x = 3, y = ppps[3], pch = 16, cex = 1.5, col = adjustcolor(3, 0.5))
     load(file = paste0(datasets[i], "_eigen.RData"))
     points(x = 4, y = ppps[3], pch = 16, cex = 1.5, col = adjustcolor(4, 0.5))
}