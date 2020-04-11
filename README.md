# Graphical Component Analysis

This is an R implementation of the algorithm described in <a href="https://github.com/nataliedoss/Thesis/main.pdf" download>Chapter 6: Graphical component analysis for latent signal detection</a>. A sample script to test the algorithm is below. 

## External dependencies

[orthopolynom](https://cran.r-project.org/web/packages/orthopolynom/index.html)

[igraph](https://igraph.org/r/)

[huge](https://cran.r-project.org/web/packages/huge/index.html)

## Example
```
source("data_create.r")
source("gca.r")
source("eval.r")

n = 2000
test_size = 500
d <- 3
type = "cycles"
a <- rep(c(1, 2, 3), d)
dep = 0.3
dat = gen_dat_mixed(n, test_size, d, type, dep, a)

x_train = dat[[1]]
x_test = dat[[2]]
W_true = dat[[3]]
s_true_train = dat[[4]]
s_true_test = dat[[5]]
L_true = dat[[6]]

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
seed = 3 

# Run gca
gca_cycles = gca_eval(x_train, x_test, W_true, s_true_train, s_true_test, L_true, "plots_small_gca", 
                      m1, m2, rho, cores, len_lam, niter, tol, init, seed, 
                      N, max_component_size, TRUE)

```
