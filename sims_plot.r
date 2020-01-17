###############################################################################
# Make plot
library(ggplot2)
source("data_gen_sims.r")


###############################################################################
# Plot the simulation results

# Read in the log likelihood results
list_llks_gca <- list_llks_tca <- list_llks_ica <- list_llk_stats <- vector("list", length(vec_dep))

for (dep in 1:length(vec_dep)) {
  
  mat_llks_gca <- mat_llks_tca <- mat_llks_ica <- matrix(rep(NA, nsims*length(vec_n)), nrow = length(vec_n), ncol = nsims)
  llk_sd_gca <- llk_sd_tca <- llk_sd_ica <- rep(NA, length(vec_n))
  
  for (n in 1:length(vec_n)) {
    
    for (sim in 1:nsims) {
      mat_llks_gca[n, sim] <- read.csv(paste0("plots_gca/dep_", dep, "_n_", n, "_sim_", sim, "_llk_gca.csv"), header = TRUE)[, 2]
      mat_llks_tca[n, sim] <- read.csv(paste0("plots_tca/dep_", dep, "_n_", n, "_sim_", sim, "_llk_tca.csv"), header = TRUE)[, 2]
      mat_llks_ica[n, sim] <- read.csv(paste0("plots_ica/dep_", dep, "_n_", n, "_sim_", sim, "_llk_ica.csv"), header = TRUE)[, 2]
    }
    
    llk_sd_gca[n] <- sd(mat_llks_gca[n, ])
    llk_sd_tca[n] <- sd(mat_llks_tca[n, ])
    llk_sd_ica[n] <- sd(mat_llks_ica[n, ])
  }
  
  llk_mean_gca <- rowMeans(mat_llks_gca, na.rm=TRUE)
  list_llks_gca[[dep]] <- mat_llks_gca
  
  llk_mean_tca <- rowMeans(mat_llks_tca, na.rm=TRUE)
  list_llks_tca[[dep]] <- mat_llks_tca
  
  llk_mean_ica <- rowMeans(mat_llks_ica, na.rm=TRUE)
  list_llks_ica[[dep]] <- mat_llks_ica
  
  list_llk_stats[[dep]] <- as.data.frame(cbind(vec_n, llk_mean_gca, llk_sd_gca,
                                               llk_mean_tca, llk_sd_tca,
                                               llk_mean_ica, llk_sd_ica))
}

# Basic plot of heldout log likelihood (no sd, no ggplot2) 
pdf("plot_llk_mean_basic.pdf")
par(mfrow = c(1, 3))
for (i in 1:length(vec_dep)) {
  plot(vec_n, list_llk_stats[[i]]$llk_mean_gca, col = "red",
       ylim=c(min(list_llk_stats[[i]]$llk_mean_gca, list_llk_stats[[i]]$llk_mean_tca, list_llk_stats[[i]]$llk_mean_ica) - 100,
              max(list_llk_stats[[i]]$llk_mean_gca, list_llk_stats[[i]]$llk_mean_tca, list_llk_stats[[i]]$llk_mean_ica) + 100), 
       xlab = "Sample size", ylab = "Heldout log likelihood")
  points(vec_n, list_llk_stats[[i]]$llk_mean_tca, col = "green")
  points(vec_n, list_llk_stats[[i]]$llk_mean_ica, col = "blue")
  legend("bottomright", legend=c("GCA", "TCA", "ICA"), pch=1, col=c("red", "green", "blue"))
}
dev.off()



# Read in the ROC results
list_rocs_gca <- list_rocs_tca <- vector("list", length(vec_dep))
for (dep in 1:length(vec_dep)) {
  mat_rocs_gca <- mat_rocs_tca <- matrix(rep(NA, 2*nsims*length(vec_n)), nrow = nsims*length(vec_n), ncol = 2)
  for (n in 1:length(vec_n)) {
    for (sim in 1:nsims) {
      mat_rocs_gca[(n-1)*length(vec_n) + sim, ] <- read.csv(paste0("plots_gca/dep_", dep, "_n_", n, "_sim_", sim, "_roc_gca.csv"), header = TRUE)[, 2]
      mat_rocs_tca[(n-1)*length(vec_n) + sim, ] <- read.csv(paste0("plots_tca/dep_", dep, "_n_", n, "_sim_", sim, "_roc_tca.csv"), header = TRUE)[, 2]
    }
  }
  
  list_rocs_gca[[dep]] <- mat_rocs_gca
  list_rocs_tca[[dep]] <- mat_rocs_tca
}





# ROC plot
pdf("plot_roc.pdf")
par(mfrow = c(1, 2))
for (i in 2:length(vec_dep)) {
  plot(list_rocs_gca[[i]], col = "red", 
       xlab = "FPR", ylab = "TPR", 
       xlim = c(0, 1), ylim = c(0, 1))
  points(list_rocs_tca[[i]], col = "green")
  legend("bottomright", legend=c("GCA", "TCA"), pch=1, col=c("red", "green"))
}
dev.off()


# Save the ROC information
save(list_rocs_tca[[2]], file="list_rocs_tca_2.RData")
save(list_rocs_tca[[3]], file="list_rocs_tca_3.RData")
save(list_rocs_gca[[2]], file="list_rocs_gca_2.RData")
save(list_rocs_gca[[3]], file="list_rocs_gca_3.RData")


# Plot llk_mean using ggplot2
#p <- lapply(1:length(vec_dep), function(dep) ggplot(list_llk_stats[[dep]], aes(x=vec_n, y=llk_mean_gca)) +
#              geom_errorbar(aes(ymin=llk_mean_gca-llk_sd_gca, ymax=llk_mean_gca+llk_sd_gca), width=.1) +
#              geom_errorbar(aes(ymin=llk_mean_ica-llk_sd_ica, ymax=llk_mean_ica+llk_sd_ica), width=.1) +
#              geom_line() + geom_point() +
#              labs(x="Sample size", y="Heldout log likelihood", title = paste0("Dependence level:", vec_dep[dep])))

#pdf(paste0("plot_llk_mean.pdf"))
#multiplot(p[[1]], p[[2]], p[[3]], cols=length(vec_dep))
#dev.off()


