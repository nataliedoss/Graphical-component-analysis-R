###############################################################################
# Test GCA and run sims
library(ggplot2)
source("data_gen_sims.r")
source("eval.r")


###############################################################################

# score-matching parameters
m1 = 4
m2 = 1
rho = 1
cores = 1
len_lam = 10
niter = 30
tol = 0.001
init = "random"
seed = 3 # Not used but leave in. 
N = 100 # for density estimation 
max_component_size = 3
edge_threshold = 1.0

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
sim_gca <- gca_eval(list_x_train[[dep]][[n]][[sim]], list_x_test[[dep]][[n]][[sim]], 
           list_W_true[[dep]][[n]][[sim]], 
           list_s_true_train[[dep]][[n]][[sim]], list_s_true_test[[dep]][[n]][[sim]], 
           list_L_true[[dep]][[n]][[sim]], paste0("plots_gca/dep_", dep, "_n_", n, "_sim_", sim), 
           m1, m2, rho, cores, len_lam, niter, tol, init, seed, N, 
           max_component_size, edge_threshold, TRUE)
write.csv(sim_gca[1], file=paste0("plots_gca/dep_", dep, "_n_", n, "_sim_", sim, "_llk_gca.csv"))
write.csv(sim_gca[c(2, 3)], file=paste0("plots_gca/dep_", dep, "_n_", n, "_sim_", sim, "_roc_gca.csv"))

# ICA
sim_ica <- ica_eval(list_x_train[[dep]][[n]][[sim]], list_x_test[[dep]][[n]][[sim]],
                    list_W_true[[dep]][[n]][[sim]],
                    list_s_true_train[[dep]][[n]][[sim]], list_s_true_test[[dep]][[n]][[sim]], 
                    list_L_true[[dep]][[n]][[sim]], paste0("plots_ica/dep_", dep, "_n_", n, "_sim_", sim), 
                    N, TRUE)
write.csv(sim_ica, file=paste0("plots_ica/dep_", dep, "_n_", n, "_sim_", sim, "_llk_ica.csv"))

# TCA
W_est = as.matrix(read.table(paste0("tca/Store_estimates/W_est_dep_", dep, "_n_", n, "_sim_", sim, ".csv"), sep=","))
tree_est = read.table(paste0("tca/Store_estimates/tree_est_dep_", dep, "_n_", n, "_sim_", sim, ".txt"),
                       sep=",", header = TRUE)
sim_tca <- tca_eval(list_x_train[[dep]][[n]][[sim]], list_x_test[[dep]][[n]][[sim]],
                     list_W_true[[dep]][[n]][[sim]],
                     list_s_true_train[[dep]][[n]][[sim]], list_s_true_test[[dep]][[n]][[sim]],
                     list_L_true[[dep]][[n]][[sim]], paste0("plots_tca/dep_", dep, "_n_", n, "_sim_", sim),
                     N, TRUE, W_est, tree_est)
write.csv(sim_tca[1], file=paste0("plots_tca/dep_", dep, "_n_", n, "_sim_", sim, "_llk_tca.csv"))
write.csv(sim_tca[c(2, 3)], file=paste0("plots_tca/dep_", dep, "_n_", n, "_sim_", sim, "_roc_tca.csv"))


#}





