###############################################################################
# Test GCA and run sims
library(ggplot2)
source("data_gen_sims.r")
source("gca.r")
source("eval.r")


###############################################################################

# score-matching parameters
m1 = 4
m2 = 1
rho = 1
cores = 1

len_lam = 20
niter = 30
tol = 0.001
N = 100 # for density estimation in computing perm_vec
max_component_size = 3
init = "random"
seed = 3 # Not used right now but is necessary as placeholder so don't comment out


# Obtain the arguments from the command line
args = commandArgs(trailingOnly=TRUE)
ind <- as.integer(args)

# Extract the indices of the vectors of dependence levels and sample sizes
#for (ind in 1:(length(vec_dep)*length(vec_n)*nsims)) {
dep = ceiling(ind/ (length(vec_n)*nsims))
n = ceiling((ind %% (length(vec_n)*nsims)) / nsims)
if (n == 0) n = length(vec_n)
sim = (ind %% nsims)
if (sim == 0) sim = nsims
print(c(dep, n, sim))


# GCA
llk_gca <- gca_eval(list_x_train[[dep]][[n]][[sim]], list_x_test[[dep]][[n]][[sim]], 
           list_W_true[[dep]][[n]][[sim]], 
           list_s_true_train[[dep]][[n]][[sim]], list_s_true_test[[dep]][[n]][[sim]], 
           list_L_true[[dep]][[n]][[sim]], paste0("plots_gca/dep_", dep, "_n_", n, "_sim_", sim), 
           m1, m2, rho, cores, len_lam, niter, tol, init, seed, N, max_component_size, TRUE)
write.csv(llk_gca, file=paste0("plots_gca/dep_", dep, "_n_", n, "_sim_", sim, "_llk_gca.csv"))

# ICA
llk_ica <- ica_eval(list_x_train[[dep]][[n]][[sim]], list_x_test[[dep]][[n]][[sim]],
                    list_W_true[[dep]][[n]][[sim]],
                    list_s_true_train[[dep]][[n]][[sim]], list_s_true_test[[dep]][[n]][[sim]], 
                    list_L_true[[dep]][[n]][[sim]], paste0("plots_ica/dep_", dep, "_n_", n, "_sim_", sim), 
                    N, TRUE)
write.csv(llk_ica, file=paste0("plots_ica/dep_", dep, "_n_", n, "_sim_", sim, "_llk_ica.csv"))

# TCA
# W_est = as.matrix(read.table(paste0("tca/Store_estimates/W_est_dep_", dep, "_n_", n, "_sim_", sim, ".csv"), sep=","))
# tree_est = read.table(paste0("tca/Store_estimates/tree_est_dep_", dep, "_n_", n, "_sim_", sim, ".txt"), 
#                       sep=",", header = TRUE)
# llk_tca <- tca_eval(list_x_train[[dep]][[n]][[sim]], list_x_test[[dep]][[n]][[sim]],
#                     list_W_true[[dep]][[n]][[sim]],
#                     list_s_true_train[[dep]][[n]][[sim]], list_s_true_test[[dep]][[n]][[sim]], 
#                     list_L_true[[dep]][[n]][[sim]], vec_names_tca[sim],
#                     N, FALSE, W_est, tree_est)
# write.csv(llk_tca, file=paste0("plots_tca/dep_", dep, "_n_", n, "_sim_", sim, "_llk_tca.csv"))
# 




###############################################################################
# Plot the simulation results

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




