###############################################################################
# A small experiment
source("data_create.r")
source("gca.r")
source("eval.r")



###############################################################################
set.seed(11)
n = 2000
test_size = 500
d <- 3
type = "tree"
a <- rep(c(1, 2, 3), d)
dep = 0.3
dat = gen_dat_mixed(n, test_size, d, type, dep, a)

x_train = dat[[1]]
x_test = dat[[2]]
W_true = dat[[3]]
s_true_train = dat[[4]]
s_true_test = dat[[5]]
L_true = dat[[6]]






###############################################################################
# Test gca, tca, ica
# score-matching parameters
m1 = 4
m2 = 1
rho = 1
cores = 1

len_lam = 10
niter = 30
tol = 0.01
N = 100 # for density estimation 
max_component_size = 3
init = "random"
seed = 3 # Not used but leave in. 
# It's a placeholder parameter in case you want to run gca on the same seed. 
# It's commented out in that code now.


# Run gca
gca_cycles = gca_eval(x_train, x_test, W_true, s_true_train, s_true_test, L_true, "plots_small_gca", 
                      m1, m2, rho, cores, len_lam, niter, tol, init, seed, 
                      N, max_component_size, TRUE)


# To run TCA, save x_train, run the TCA Matlab code on it
# Then read in the W_est and tree_est from tca file
write.table(t(x_train), "tca/x_train.csv", row.names=FALSE, col.names=FALSE, sep=",")
W_est = as.matrix(read.table("tca/W.csv", sep=","))
tree_est <- read.table("tca/tree_est.txt", sep=",", header = TRUE)

# Run tca
tca_cycles = tca_eval(x_train, x_test, W_true, s_true_train, s_true_test, L_true, "plots_small_tca", 
                     N, TRUE, W_est, tree_est) 

# Run ica
ica_cycles = ica_eval(x_train, x_test, W_true, s_true_train, s_true_test, L_true, "plots_small_ica", 
                      N, TRUE) 








###############################################################################
# Compare GCA, TCA, and ICA via plots

load("plots_small_gca/s_est_train_perm.RData")
s_est_train_gca = s_est_train_perm
load("plots_small_tca/s_est_train_perm.RData")
s_est_train_tca = s_est_train_perm
load("plots_small_ica/s_est_train_perm.RData")
s_est_train_ica = s_est_train_perm

p1 <- lapply(1:d, function(i) ggplot(as.data.frame(s_true_train), aes(x=s_true_train[, i])) + 
               geom_histogram(aes(y=..density..), colour="black", fill="white", bins=20) +
               geom_density(alpha=.2, fill="blue") + labs(x="", y="Density") + ggtitle("True source"))

p2 <- lapply(1:d, function(i) ggplot(as.data.frame(s_est_train_gca), aes(x=s_est_train_gca[, i])) + 
               geom_histogram(aes(y=..density..), colour="black", fill="white", bins=20) +
               geom_density(alpha=.2, fill="blue") + labs(x="", y="Density") + ggtitle("GCA"))

p3 <- lapply(1:d, function(i) ggplot(as.data.frame(s_est_train_tca), aes(x=s_est_train_tca[, i])) + 
               geom_histogram(aes(y=..density..), colour="black", fill="white", bins=20) +
               geom_density(alpha=.2, fill="blue") + labs(x="", y="Density") + ggtitle("TCA"))

p4 <- lapply(1:d, function(i) ggplot(as.data.frame(s_est_train_ica), aes(x=s_est_train_ica[, i])) + 
               geom_histogram(aes(y=..density..), colour="black", fill="white", bins=20) +
               geom_density(alpha=.2, fill="blue") + labs(x="", y="Density") + ggtitle("ICA"))  

pdf("hist_marginals_gca_tca_ica.pdf")
multiplot(p1[[1]], p2[[1]], p3[[1]], p4[[1]],
          p1[[2]], p2[[2]], p3[[2]], p4[[2]],
          p1[[3]], p2[[3]], p3[[3]], p4[[3]], cols = 4)
dev.off()






###############################################################################
# Testing what the dependence parameter v does
# Default for u is 0.1
#library(huge)
#n = 1000
#d = 3
#deps <- seq(0.0, 2.0, by=0.1)
#tmp <- rep(NA, length(deps))
#for (i in 1:length(deps)) {
#  L = huge.generator(n = n, d = d, graph = "cluster", v = deps[i], g=d/3)
#  tmp[i] <- L$sigma[1,2]
#}
#plot(tmp)


