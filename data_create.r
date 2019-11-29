###############################################################################
# Functions to generate data for GCA testing
library(huge)
library(feather)


###############################################################################
# A transformation function
f <- function(x, alpha) {
  return(sign(x) * abs(x)^alpha)
}

# A function to generate a dataset according to our model where x = As
gen_dat_mixed <- function(n, test_size, d, type, dep, a) {

  N = n + test_size # train + test set size
  
  if (type == "cycles") {
    L_true = huge.generator(n = N, d = d, graph = "cluster", v = dep, g = d/3)
    s_true_full = L_true$data
  }
  else if (type == "tree") {
    L_true = huge.generator(n = N, d = d, graph="scale-free", v = dep)
    s_true_full = L_true$data
  }
  else {
    # for independent data, L_true is just a placeholder that isn't used
    L_true = huge.generator(n = N, d = d, graph="scale-free", v = dep) 
    s_true_full = matrix(rnorm(N*d), nrow = N, ncol = d)
  }

  # Regardless of which type of data was chosen, now transform it via the npn:
  for(i in 1:d) {
    s_true_full[,i] <- f(s_true_full[,i], a[i]*0.2)
  }
  
  # Create a mixing matrix at random
  A_true = svd(matrix(rnorm(d^2), d, d))$u
  W_true = solve(A_true)
  
  # Mix the data via the mixing matrix
  ind <- sample(nrow(s_true_full), test_size)
  s_true_train = s_true_full[-ind, ]
  s_true_test = s_true_full[ind, ]
  x_train = s_true_train %*% t(A_true)
  x_test = s_true_test %*% t(A_true)
  
  return(list(x_train, x_test, W_true, s_true_train, s_true_test, L_true)) 
}


# Returns nsims datasets of mixed data (for a fixed n)
gen_dat_nsims <- function(nsims, n, test_size, d, type, dep, a) {
  return(replicate(nsims, gen_dat_mixed(n, test_size, d, type, dep, a)))
}

# Returns length(vec_n) collections 
# Each collection has nsims datasets of mixed data
gen_dat_n_nsims <- function(nsims, vec_n, test_size, d, type, dep, a) {
  return(lapply(vec_n, gen_dat_nsims, nsims=nsims, test_size=test_size, d=d, type=type, dep=dep, a=a))
}




