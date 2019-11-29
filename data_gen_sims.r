###############################################################################
# Generate data for GCA tests
library(huge)
source("data_create.r")




###############################################################################
# Generate the datasets we will use for simulations
set.seed(20)
nsims = 10
vec_n <- seq(500, 5000, length=10)
vec_dep <- c(0.0, 0.5, 0.5)
vec_type = c("indep", "tree", "cycles")
test_size = 500
d <- 15
a <- rep(c(1, 2, 3), d) 

dat = list_x_train = list_x_test = list_W_true = list_s_true_train = list_s_true_test = list_L_true = vector("list", length(vec_dep))
for (i in 1:length(vec_dep)) {
  for (j in 1:length(vec_n)) {
    dat[[i]][[j]] = gen_dat_nsims(nsims, vec_n[j], test_size, d, vec_type[i], vec_dep[i], a)
    list_x_train[[i]][[j]] = dat[[i]][[j]][1, ]
    list_x_test[[i]][[j]] = dat[[i]][[j]][2, ]
    list_W_true[[i]][[j]] = dat[[i]][[j]][3, ]
    list_s_true_train[[i]][[j]] = dat[[i]][[j]][4, ]
    list_s_true_test[[i]][[j]] = dat[[i]][[j]][5, ]
    list_L_true[[i]][[j]] = dat[[i]][[j]][6, ]
  }
}

# Save the x_train data as csv's so TCA can run on them
# dir.create("tca/sim_data")
# for (i in 1:length(vec_dep)) {
#  for (j in 1:length(vec_n)) {
#    for (k in 1:nsims) {
#      write.table(t(list_x_train[[i]][[j]][[k]]),
#                  paste0("tca/Store_sim_data/x_train_dep_", i, "_n_", j, "_sim_", k, ".csv"),
#                  row.names=FALSE, col.names=FALSE, sep=",")
#    }
#  }
# }



# Double check whether dependence is there in s_true
#library(ggplot2)
#s1 = list_s_true_train[[1]][[1]][[1]]
#s2 = list_s_true_train[[2]][[1]][[1]]
#p <- ggplot(as.data.frame(s1), aes(x = s1[, 1], y = s1[, 3])) + geom_density_2d()
#plot(p)



