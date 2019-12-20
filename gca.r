library(orthopolynom)
library(igraph)
library(parallel)
library(ggplot2)
library(plot3D)
library(grid)
library(fastICA)
source("suff_stats_unbounded.r")
source("admm.r")




###############################################################################
# This is the loss in W, i.e., in eta, which is all we need.

loss_W = function(x, W, theta, m1, m2, rho, cores) {
  n = dim(x)[1]
  d = dim(x)[2]
  s = x %*% t(W)
  
  K = Kvec2(s, m1, m2, cores)
  gr = Gamma.ret2(s, m1, m2, cores)
  Gam = mclapply(1:d,function(x) gr[[x]]$v %*% (((1/n)*gr[[x]]$d^2)*t(gr[[x]]$v)), mc.cores=cores)

  objective_list = mclapply(1:d, function(x) {0.5 * t(theta[[x]]) %*% Gam[[x]] %*% theta[[x]] + t(theta[[x]]) %*% K[[x]]}, 
                            mc.cores=cores)
  objective = Reduce('+', objective_list)
  
  return(objective)
}

penalty_group_sparsity <- function(theta, m1, m2) {
  d = length(theta)
  penalty_group = vector("list", d)
  for (i in 1:d) {
    for (j in 1:d) {
      if (i == j) {
        penalty_group[[i]][j] = sqrt(sum((theta[[i]][((i-1)*m2^2+1):((i-1)*m2^2+m1)])^2))
      }
      else if (i < j) {
        penalty_group[[i]][j] = sqrt(sum((theta[[i]][(m1+(j-2)*m2^2+1):(m1+(j-2)*m2^2+m2^2)])^2))
      }
      else {
        penalty_group[[i]][j] = sqrt(sum((theta[[i]][((j-1)*m2^2+1):(j*m2^2)])^2))
      }
    }
  }
  return(sum(Reduce("+", penalty_group)))
}


loss = function(x, W, theta, m1, m2, rho, cores) {
  
  loss_sm = loss_W(x, W, theta, m1, m2, rho, cores)
  penalty = penalty_group_sparsity(theta, m1, m2)
  
  return(loss_sm + penalty)
}


# Make a function to return the G_eta
Givens_rotate = function(eta, u, v) {
  G_eta = diag(d)
  G_eta[u, u] = cos(eta)
  G_eta[u, v] = -sin(eta)
  G_eta[v, u] = sin(eta)
  G_eta[v, v] = cos(eta)
  return(G_eta)
}

# Make the objective as a function of eta.
loss_eta = function(eta, u, v, x, W, theta, m1, m2, rho, cores) {
  G_eta = Givens_rotate(eta, u, v)
  objective = loss_W(x, W %*% G_eta, theta, m1, m2, rho, cores)
  return(objective)
}


# Optimize in theta
opt_theta = function(x, W, m1, m2, rho, cores, len_lam) {
  n = dim(x)[1]
  d = dim(x)[2]
  
  # obtain training sample
  ind <- sample(c(T, F), n, replace = TRUE, prob = c(0.6, 0.4))
  x_train <- x[ind, ]
  x_test <- x[!ind, ]
  s_train = x_train %*% t(W)
  
  lammax=max(abs(unlist(Kvec2(s_train,m1,m2,cores))))
  lamseq=exp(seq(from=log(.001),to=log(lammax),length.out=len_lam))
  admm=admmSel(s_train,lamseq,m1,m2,rho,cores)
  
  loss_theta <- lapply(1:length(lamseq), function(i) {loss(x_test, W, admm[[i]]$zz, m1, m2, rho, cores)})
  lam_opt <- lamseq[which.min(loss_theta)]
  
  # Use optimal lambda to estimate theta on whole dataset
  s = x %*% t(W)
  admm_opt=admmSel(s,c(lam_opt),m1,m2,rho,cores)
  theta_opt = admm_opt[[1]]$zz

  return(list(theta_opt, lam_opt))
}

# Optimize one rotation of W using grid search. 
opt_W_gs <- function(u, v, x, W, theta, m1, m2, rho, cores, grid_size, iter) {
  eta_grid <- seq(-pi, pi, length = grid_size)
  loss_eta_grid <- rep(NA, length = grid_size)
  loss_eta_grid <- unlist(lapply(1:length(eta_grid), function(j) 
    loss(x, W %*% Givens_rotate(eta_grid[j], u, v), theta, m1, m2, rho, cores)))
  eta_opt = eta_grid[which.min(loss_eta_grid)]
  #pdf(paste0(iter, "_loss_", u, "_", v, ".pdf"))
  #plot(eta_grid, loss_eta_grid)
  #dev.off()
  W_new = W %*% Givens_rotate(eta_opt, u, v)
  return(list(W_new, eta_opt, min(loss_eta_grid), loss_eta_grid, 0, 0))
}

# Optimize one rotation of W using quasi-Newton method
opt_W = function(u, v, x, W, theta, m1, m2, rho, cores) {
  opt = optim(0.0, loss_eta, u = u, v = v, x = x, W = W, theta = theta,
              m1 = m1, m2 = m2, rho = rho, cores = cores, 
              lower = -pi, upper = pi, method = "L-BFGS-B")
  W_new = W %*% Givens_rotate(opt$par, u, v)
  return(list(W_new, opt))
}



###############################################################################
# The GCA algorithm
gca <- function(x, m1, m2, rho, cores, len_lam, niter, tol, init, seed) {
  
  d <- dim(x)[2]
  
  # Initialize W.
  if (init == "ica") {
    W_init <- solve(fastICA(x, n.comp=d)$A)
  }  else if (init == "random") {
    #set.seed(seed)
    W_init = svd(matrix(rnorm(d^2), d, d))$u
  } else W_init = diag(d)
  
  #
  nsubiter <- niter*d*(d-1)/2
  t = 1
  counter_theta = 1
  counter_eta = 1
  W = W_init # Will be updated
  loss_change = 1.0
  
  time_vec_theta = data.frame(matrix(rep(NA, 5*niter), nrow = niter, ncol = 5))
  time_vec_eta = data.frame(matrix(rep(NA, 5*nsubiter), nrow = nsubiter, ncol = 5))
  
  loss_vec = c(.01, rep(NA, niter-1))
  loss_vec_theta = c(.01, rep(NA, niter-1))
  loss_vec_eta = c(.01, rep(NA, nsubiter-1))
  
  opt_vec_theta = rep(NA, niter)
  opt_mat_eta = data.frame(matrix(rep(NA, 6*niter), nrow = niter, ncol = 6)) 
  
  list_pairs <- NULL
  i = 1
  for (u in 1:(d-1)) {
    for (v in (u+1):d) {
      list_pairs[[i]] <- c(u, v)
      i = i + 1
    }
  }
  
  while (loss_change > tol && t <= niter) {
    
    # Update theta
    counter_theta = counter_theta + 1
    start_time_theta <- proc.time()

    opt_theta_fn <- opt_theta(x, W, m1, m2, rho, cores, len_lam) # UPDATE THETA
    
    time_vec_theta[counter_theta, ] <- proc.time() - start_time_theta
    theta = opt_theta_fn[[1]]
    loss_vec_theta[counter_theta] = loss(x, W, theta, m1, m2, rho, cores)
    opt_vec_theta[counter_theta] = c(opt_theta_fn[[2]])
    print(paste0("time for theta opt: ", time_vec_theta[counter_theta, 3], 
                 " lam_opt: ", opt_theta_fn[[2]]))
    
    # Update W
    for (j in 1:length(list_pairs)) {
      counter_eta = counter_eta + 1
      pair <- list_pairs[[j]]
      u <- pair[1]; v <- pair[2]
      
      start_time_eta <- proc.time()
      opt_W_fn = opt_W(u, v, x, W, theta, m1, m2, rho, cores) # UPDATE W (ONE ANGLE)
      end_time_eta <- proc.time() - start_time_eta
      time_vec_eta[counter_eta, ] <- end_time_eta
      
      # Or update optimize the angle via a grid search
      #start_time_eta_gs <- proc.time()
      #opt_W_fn_gs = opt_W_gs(u, v, x, W, theta, m1, m2, rho, cores, 1000, t)
      #end_time_eta_gs <- proc.time() - start_time_eta_gs
      
      W = opt_W_fn[[1]]
      loss_vec_eta[counter_eta] = loss(x, W, theta, m1, m2, rho, cores)
      opt_mat_eta[counter_eta, ] = c(unlist(opt_W_fn[2:6])) 
      print(paste0("pair: ", u, " ", v, " time for eta opt from optim: ", time_vec_eta[counter_eta, 3], 
                   " opt_eta: ", opt_W_fn[[2]][[1]]))
                   #" time for eta opt from gs: ", end_time_eta_gs[3], " opt_eta_gs: ", opt_W_fn_gs[[2]]))
      
    }
    t = t + 1
    loss_vec[t] = loss(x, W, theta, m1, m2, rho, cores)
    loss_change = abs((loss_vec[t] - loss_vec[t-1])/loss_vec_eta[t-1])
  }
  
  # Return all
  return(list(W, theta, W_init, 
              time_vec_theta, time_vec_eta, 
              loss_vec, loss_vec_theta, loss_vec_eta, 
              opt_vec_theta, opt_mat_eta))
}









