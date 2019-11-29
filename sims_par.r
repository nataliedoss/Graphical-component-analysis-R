###############################################################################
# Test GCA and run sims
library(ggplot2)
library(foreach)
library(doParallel)
source("data_gen_sims.r")
source("gca.r")
source("eval.r")


###############################################################################

# score-matching parameters
m1 = 4
m2 = 1
rho = 1
cores = 1 # keep at one for admm

len_lam = 10
niter = 10
tol = 0.01
N = 100 # for density estimation in computing perm_vec
max_component_size = 3
init = "identity"
seed = 3 # Not really used now but is necessary as placeholder so don't comment out


# Obtain the arguments from the command line
args = commandArgs(trailingOnly=TRUE)
ind <- as.integer(args)

# Extract the indices of the vectors of dependence levels and sample sizes
#for (ind in 1:(length(vec_dep)*length(vec_n))) {
dep = ceiling(ind/length(vec_n))
n = ind %% length(vec_n)
if (n == 0) n = length(vec_n)

print(c(n, dep))


# Create names vector
vec_names_gca <- paste0("plots_gca/dep_", dep, "_n_", n, "_sims_", 1:nsims)
vec_names_tca <- paste0("plots_tca/dep_", dep, "_n_", n, "_sims_", 1:nsims)
vec_names_ica <- paste0("plots_ica/dep_", dep, "_n_", n, "_sims_", 1:nsims)

# Run GCA in parallel
registerDoParallel(cores=nsims)
print(getDoParWorkers())
list_llk_gca <- foreach(i=1:nsims) %dopar% 
  gca_eval(list_x_train[[dep]][[n]][[i]], list_x_test[[dep]][[n]][[i]], 
           list_W_true[[dep]][[n]][[i]], 
           list_s_true_train[[dep]][[n]][[i]], list_s_true_test[[dep]][[n]][[i]], 
           list_L_true[[dep]][[n]][[i]], vec_names_gca[i], 
           m1, m2, rho, cores, len_lam, niter, tol, init, seed, N, max_component_size, FALSE)
write.csv(unlist(list_llk_gca), file=paste0("plots_gca/dep_", dep, "_n_", n, "_llk_gca.csv"))


# 
# # Run ICA in loop for now
# list_llk_ica <- vector(mode = "list", length = nsims)
# for (i in 1:nsims) {
#   list_llk_ica[[i]] <- ica_eval(list_x_train[[dep]][[n]][[i]], list_x_test[[dep]][[n]][[i]],
#                                 list_W_true[[dep]][[n]][[i]], 
#                                 list_s_true_train[[dep]][[n]][[i]], list_s_true_test[[dep]][[n]][[i]], 
#                                 list_L_true[[dep]][[n]][[i]], vec_names_ica[i], 
#                                 N, FALSE)
# }
# write.csv(unlist(list_llk_ica), file=paste0("plots_ica/dep_", dep, "_n_", n, "_llk_ica.csv"))
# 
#  
# # Run TCA in loop for now
# list_llk_tca <- vector(mode = "list", length = nsims)
# for (sim in 1:nsims) {
#   W_est = as.matrix(read.table(paste0("tca/Estimates/W_est_dep_", dep, "_n_", n, "_sim_", sim, ".csv"), sep=","))
#   tree_est = read.table(paste0("tca/Estimates/tree_est_dep_", dep, "_n_", n, "_sim_", sim, ".txt"), 
#                         sep=",", header = TRUE)
#   
#   list_llk_tca[[sim]] <- tca_eval(list_x_train[[dep]][[n]][[sim]], list_x_test[[dep]][[n]][[sim]], 
#                                   list_W_true[[dep]][[n]][[sim]],
#                                   list_s_true_train[[dep]][[n]][[sim]], list_s_true_test[[dep]][[n]][[sim]], 
#                                   list_L_true[[dep]][[n]][[sim]], vec_names_tca[sim],
#                                   N, FALSE, W_est, tree_est)
# }
# write.csv(unlist(list_llk_tca), file=paste0("plots_tca/dep_", dep, "_n_", n, "_llk_tca.csv"))

#}
# 
# 
# # Run the simulations
# #list_llk_gca <- mcmapply(gca_eval, list_x_train[[dep]][[n]], list_x_test[[dep]][[n]], 
# #                         list_W_true[[dep]][[n]], 
# #                         list_s_true_train[[dep]][[n]], list_s_true_test[[dep]][[n]], 
# #                         list_L_true[[dep]][[n]], vec_names_gca, 
# #                         MoreArgs = list(m1, m2, rho, cores, len_lam, niter, tol, 
# #                                         init, seed, N, max_component_size, TRUE),
# #                         mc.preschedule = TRUE, mc.cores = nsims)
# #write.csv(unlist(list_llk_gca), file=paste0("plots_gca/dep_", dep, "_n_", n, "_llk_gca.csv"))
# 
# 
# #list_llk_ica <- mcmapply(ica_eval, list_x_train[[dep]][[n]], list_x_test[[dep]][[n]], 
# #                         list_W_true[[dep]][[n]], 
# #                         list_s_true_train[[dep]][[n]], list_s_true_test[[dep]][[n]], 
# #                         list_L_true[[dep]][[n]], vec_names_ica, 
# #                         MoreArgs = list(N, TRUE), mc.preschedule = TRUE, mc.cores = nsims)
# #write.csv(unlist(list_llk_ica), file=paste0("plots_ica/dep_", dep, "_n_", n, "_llk_ica.csv"))
# 
# #}
# 
# 
# 
# 
# ###############################################################################
# # Plot the simulation results
# 
# list_llks_gca <- list_llks_tca <- list_llks_ica <- list_llk_stats <- vector("list", length(vec_dep))
# 
# for (dep in 1:length(vec_dep)) {
# 
#   mat_llks_gca <- matrix(rep(NA, nsims*length(vec_n)), nrow = length(vec_n), ncol = nsims)
#   llk_sd_gca <- rep(NA, length(vec_n))
# 
#   mat_llks_tca <- matrix(rep(NA, nsims*length(vec_n)), nrow = length(vec_n), ncol = nsims)
#   llk_sd_tca <- rep(NA, length(vec_n))
# 
#   mat_llks_ica <- matrix(rep(NA, nsims*length(vec_n)), nrow = length(vec_n), ncol = nsims)
#   llk_sd_ica <- rep(NA, length(vec_n))
# 
#   for (n in 1:length(vec_n)) {
#     mat_llks_gca[n, ] <- read.csv(paste0("plots_gca/dep_", dep, "_n_", n, "_llk_gca.csv"), header = TRUE)[, 2]
#     llk_sd_gca[n] <- sd(mat_llks_gca[n, ])
# 
#     mat_llks_tca[n, ] <- read.csv(paste0("plots_tca/dep_", dep, "_n_", n, "_llk_tca.csv"), header = TRUE)[, 2]
#     llk_sd_tca[n] <- sd(mat_llks_tca[n, ])
# 
#     mat_llks_ica[n, ] <- read.csv(paste0("plots_ica/dep_", dep, "_n_", n, "_llk_ica.csv"), header = TRUE)[, 2]
#     llk_sd_ica[n] <- sd(mat_llks_ica[n, ])
#   }
#   llk_mean_gca <- rowMeans(mat_llks_gca)
#   list_llks_gca[[dep]] <- mat_llks_gca
# 
#   llk_mean_tca <- rowMeans(mat_llks_tca)
#   list_llks_tca[[dep]] <- mat_llks_tca
# 
#   llk_mean_ica <- rowMeans(mat_llks_ica)
#   list_llks_ica[[dep]] <- mat_llks_ica
# 
#   list_llk_stats[[dep]] <- as.data.frame(cbind(vec_n, llk_mean_gca, llk_sd_gca,
#                                                llk_mean_tca, llk_sd_tca,
#                                                llk_mean_ica, llk_sd_ica))
# }
# 
# 
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
# # 
# # 
# # Instead do basic plotting:
# pdf("plot_llk_mean_basic.pdf")
# par(mfrow = c(1, 3))
# for (i in 1:length(vec_dep)) {
#   plot(list_llk_stats[[i]]$llk_mean_gca, col = "red",
#        ylim=c(min(list_llk_stats[[i]]$llk_mean_gca, list_llk_stats[[i]]$llk_mean_tca, list_llk_stats[[i]]$llk_mean_ica) - 100,
#               max(list_llk_stats[[i]]$llk_mean_gca, list_llk_stats[[i]]$llk_mean_tca, list_llk_stats[[i]]$llk_mean_ica) + 100))
#   points(list_llk_stats[[i]]$llk_mean_tca, col = "green")
#   points(list_llk_stats[[i]]$llk_mean_ica, col = "blue")
#   legend("topleft", legend=c("GCA", "TCA", "ICA"), pch=1, col=c("red", "green", "blue"))
# }
# dev.off()

# 
# 
# 
# 
