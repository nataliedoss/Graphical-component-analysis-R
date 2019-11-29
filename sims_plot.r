###############################################################################
# Make plot
library(ggplot2)



###############################################################################
# Plot the simulation results

nsims = 10
vec_n <- seq(500, 5000, length=10)
vec_dep <- c(0.0, 0.3, 0.3)


list_llks_gca <- list_llks_tca <- list_llks_ica <- list_llk_stats <- vector("list", length(vec_dep))

for (dep in 1:length(vec_dep)) {
  
  mat_llks_gca <- mat_llks_tca <- mat_llks_ica <- matrix(rep(NA, nsims*length(vec_n)), nrow = length(vec_n), ncol = nsims)
  llk_sd_gca <- llk_sd_tca <- llk_sd_ica <- rep(NA, length(vec_n))
  
  for (n in 1:length(vec_n)) {
    
    for (sim in 1:nsims) {
      mat_llks_gca[n, sim] <- read.csv(paste0("plots_gca/dep_", dep, "_n_", n, "_sim_", sim, "_llk_gca.csv"), header = TRUE)[, 2]
      mat_llks_tca[n, sim] <- read.csv(paste0("plots_ica/dep_", dep, "_n_", n, "_sim_", sim, "_llk_ica.csv"), header = TRUE)[, 2]
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



# Instead do basic plotting:
pdf("plot_llk_mean_basic.pdf")
par(mfrow = c(1, 3))
for (i in 1:length(vec_dep)) {
  plot(list_llk_stats[[i]]$llk_mean_gca, col = "red",
       ylim=c(min(list_llk_stats[[i]]$llk_mean_gca, list_llk_stats[[i]]$llk_mean_tca, list_llk_stats[[i]]$llk_mean_ica) - 100,
              max(list_llk_stats[[i]]$llk_mean_gca, list_llk_stats[[i]]$llk_mean_tca, list_llk_stats[[i]]$llk_mean_ica) + 100))
  points(list_llk_stats[[i]]$llk_mean_tca, col = "green")
  points(list_llk_stats[[i]]$llk_mean_ica, col = "blue")
  legend("topleft", legend=c("GCA", "TCA", "ICA"), pch=1, col=c("red", "green", "blue"))
}
dev.off()



# # 
# # # Plot llk_mean using ggplot2
# # #p <- lapply(1:length(vec_dep), function(dep) ggplot(list_llk_stats[[dep]], aes(x=vec_n, y=llk_mean_gca)) + 
# # #              geom_errorbar(aes(ymin=llk_mean_gca-llk_sd_gca, ymax=llk_mean_gca+llk_sd_gca), width=.1) + 
# # #              geom_errorbar(aes(ymin=llk_mean_ica-llk_sd_ica, ymax=llk_mean_ica+llk_sd_ica), width=.1) + 
# # #              geom_line() + geom_point() + 
# # #              labs(x="Sample size", y="Held out log likelihood", title = paste0("Dependence level:", vec_dep[dep])))
# # 
# # #pdf(paste0("plot_llk_mean.pdf"))
# # #multiplot(p[[1]], p[[2]], p[[3]], cols=length(vec_dep))
# # #dev.off()
