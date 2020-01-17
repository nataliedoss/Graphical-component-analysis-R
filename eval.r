library(ggplot2)
library(grid)
library(fastICA)
source("gca.r")
source("density_kde.r")






###############################################################################
# A function to find the best permutation of the estimated source s_est given the true s
# This function sends the label of the estimated node to the label of the true node it's closest to
find_perm <- function(s_est, s_true, N) {
  d <- ncol(s_est)
  dist_mat <- matrix(rep(NA, d^2), nrow = d, ncol = d)
  perm_vec <- rep(NA, d)
  for (i in 1:d) {
    for (j in 1:d) {
      p_s_est <- density(s_est[, i], n = N)$y
      p_s_true <- density(s_true[, j], n = N)$y
      dist_mat[i, j] <- sum((p_s_est - p_s_true)^2) # Can do more efficiently but not bothering now
    }
  }
  perm_vec[1] <- which.min(dist_mat[, 1])
  for (i in 2:d) {
    dist_mat_new <- dist_mat
    dist_mat_new[perm_vec[1:i-1], i] <- 10000  # hacky!!
    perm_vec[i] <- which.min(dist_mat_new[, i])
  }
  
  # Return the permutation
  return(perm_vec)
}



# A function to find the best permutation given a bunch of triangles
# Each triangle has the same set of marginals in it
# This assumes that the triangle_true list is in order 1:d
find_perm_triangles <- function(s_est, s_true, N, list_triangles, list_triangles_true) {
  d <- ncol(s)
  perm_vec <- rep(NA, d)
  for (i in 1:length(list_triangles)) {
    triangle <- list_triangles[[i]]
    triangle_true <- list_triangles_true[[i]]
    s_triangle_est <- s_est[, triangle]
    s_triangle_true <- s_true[, triangle_true] 
    perm_vec_triangle <- find_perm(s_triangle_est, s_triangle_true, N)
    perm_vec[triangle_true] <- triangle[perm_vec_triangle]
  }
  return(perm_vec)
}

###############################################################################
# Functions for graph recovery

# A function to create a graph matrix given a vector theta as from GCA
create_graph_mat_gca <- function(theta, m1, m2, perm_vec, edge_threshold) {
  graph_est_mat <- diag(d)
  for (i in 1:(d-1)) {
    graph_est_mat[i, (i+1):d] = theta[[i]][((i-1)*m2 + 1 + m1):length(theta[[i]])]
  }
  # Symmetrize it - find better way!
  for (j in 1:d) {
    for (i in j:d) {
      graph_est_mat[i, j] = graph_est_mat[j, i]
    }
  }
  graph_est_mat <- graph_est_mat[perm_vec, perm_vec] # Permute
  graph_est_mat[abs(graph_est_mat) < edge_threshold] = 0 # Threshold edges
  graph_est_mat[abs(graph_est_mat) >= edge_threshold] = 1 # Threshold edges
  return(graph_est_mat)
}

# A function to create a graph matrix given a list of edges as from TCA
create_graph_mat_tca <- function(list_edges, perm_vec) {
  graph_est_mat <- diag(d)
  for (i in 1:nrow(list_edges)) {
    graph_est_mat[list_edges[i, 1], list_edges[i, 2]] <- 1
    graph_est_mat[list_edges[i, 2], list_edges[i, 1]] <- 1
  }
  return(graph_est_mat[perm_vec, perm_vec])
}

compute_roc_vals <- function(graph_est_mat, graph_true_mat) {
  tp = tn = fp = fn = 0
  d = nrow(graph_true_mat)
  # Only examine the upper triangle
  for (i in 1:(d-1)) {
    for (j in (i+1):d) {
      if (graph_est_mat[i, j] == 1 & graph_true_mat[i, j] == 1) {
        tp = tp + 1
      }
      if (graph_est_mat[i, j] == 0 & graph_true_mat[i, j] == 0) {
        tn = tn + 1
      }
      if (graph_est_mat[i, j] == 1 & graph_true_mat[i, j] == 0) {
        fp = fp + 1
      }
      if (graph_est_mat[i, j] == 0 & graph_true_mat[i, j] == 1) {
        fn = fn + 1
      }
    }
  }
  fpr = fp/(fp+tn)
  tpr = tp/(tp+fn)
  if (is.nan(fpr)) fpr = 0
  if (is.nan(tpr)) tpr = 0

  return(c(fpr, tpr))
}





###############################################################################
# A function to run gca, do the permutation, and evaluate performance
# This function takes as input some of the truth for comparison
# GCA function itself only observes data of course
gca_eval <- function(x_train, x_test, W_true, s_true_train, s_true_test, L_true, name, 
                     m1, m2, rho, cores, len_lam, niter, tol, init, seed,
                     N, max_component_size, edge_threshold, record) {
  
  start = Sys.time()
  gca_run <- gca(x_train, m1, m2, rho, cores, len_lam, niter, tol, init, seed)
  end = Sys.time() - start
  
  # Extract things from gca
  n = nrow(x_train)
  W_est = gca_run[[1]]
  theta_est = gca_run[[2]]
  W_init = gca_run[[3]]
  
  # Compute the estimated sources
  s_est_train = x_train %*% t(W_est)
  s_est_test = x_test %*% t(W_est)
  
  # Find the best permutation for the training data for plots.
  perm_vec_train = find_perm(s_est_train, s_true_train, N)
  
  # Permute training data for later comparison
  s_est_train_perm = s_est_train[, perm_vec_train]
  
  # Compute sources from initial W just to see how well one step did
  s_init = x_train %*% t(W_init)
  perm_vec_init = find_perm(s_init, s_true_train, N)
  s_init_perm = s_init[, perm_vec_init]
  
  # Find the components (just for observation purposes)
  list_components <- find_components(theta_est, m1, m2, max_component_size, edge_threshold)
  
  # Find estimated and true graph adjacency matrix
  graph_est_mat <- create_graph_mat_gca(theta_est, m1, m2, perm_vec_train, edge_threshold)
  graph_true_mat <- as.matrix(L_true$theta)
  diag(graph_true_mat) <- 1 # not necessary but makes visualizing graphs easier

  if (record) {
    dir.create(name)
    
    save(W_est, file=paste0(name, "/W_est.RData")) 
    save(theta_est, file=paste0(name, "/theta_est.RData"))
    
    save(W_init, file=paste0(name, "/W_init.RData"))
    save(W_true, file=paste0(name, "/W_true.RData"))
    
    save(perm_vec_train, file=paste0(name, "/perm_vec_train.RData"))
    
    # Easy to recompute but store anyway to save one step later
    save(s_est_train, file=paste0(name, "/s_est_train.RData"))
    save(s_est_train_perm, file=paste0(name, "/s_est_train_perm.RData"))
    
    # Save the detected components and both estimated and true graphs
    save(list_components, file=paste0(name, "/list_components.RData"))
    save(graph_est_mat, file=paste0(name, "/graph_est_mat.RData"))
    save(graph_true_mat, file=paste0(name, "/graph_true_mat.RData"))
    
    list_params = list(gca_run[[4]], gca_run[[5]], gca_run[[6]], 
                       gca_run[[7]], gca_run[[8]], gca_run[[9]], gca_run[[10]])
    cat(capture.output(print(list_params), file=paste0(name, "/list_tracking.txt")))
    
    write.csv(paste0(" n = ", n, " d = ", d, " m1 = ", m1, " m2 = ", m2,
                     " rho = ", rho, " cores = ", cores, " length of lam = ", len_lam, 
                     " niter = ", niter, " tol = ", tol, 
                     " complete algorithm time ", end),
              file=paste0(name, "/params.txt"))
    
    # Plots
    loss_vec = gca_run[[6]]
    pdf(paste0(name, "/loss.pdf"))
    plot(loss_vec)
    dev.off()
    pdf(paste0(name, "/graph_vis.pdf"))
    plot(L_true)
    dev.off()
    plot_hist_marginal(s_est_train_perm, s_true_train, x_train, d, paste0(name, "/est_marginal_gca"))
    plot_hist_bvt(s_est_train_perm, s_true_train, x_train, d, paste0(name, "/est_bvt_gca"))
    
    plot_hist_marginal(s_init_perm, s_true_train, x_train, d, paste0(name, "/init_marginal_gca"))
    plot_hist_bvt(s_init_perm, s_true_train, x_train, d, paste0(name, "/init_bvt_gca"))
  }
  
  # Return the heldout log likelihood and roc values
  llhood <- log_kde_dep(theta_est, s_est_train, s_est_test, m1, m2, max_component_size, edge_threshold) + 
    log(abs(det(W_est)))
  roc_vals <- compute_roc_vals(graph_est_mat, graph_true_mat)
  
  return(c(llhood, roc_vals))
}


tca_eval <- function(x_train, x_test, W_true, s_true_train, s_true_test, L_true, name, N, record, 
                     W_est, tree_est) {
  
  # Compute the estimated sources
  s_est_train = x_train %*% t(W_est)
  s_est_test = x_test %*% t(W_est)
  
  # Find the best permutation. 
  perm_vec_train = find_perm(s_est_train, s_true_train, N)
  s_est_train_perm = s_est_train[, perm_vec_train]
  
  # Find estimated and true graph adjacency matrix
  graph_est_mat <- create_graph_mat_tca(tree_est, perm_vec_train)
  graph_true_mat <- as.matrix(L_true$theta)
  diag(graph_true_mat) <- 1 # not necessary but makes visualizing graphs easier
  
  if (record) {
    dir.create(name)
    
    save(W_est, file=paste0(name, "/W_est.RData")) 
    save(W_true, file=paste0(name, "/W_true.RData"))
    
    save(perm_vec_train, file=paste0(name, "/perm_vec_train.RData"))
    
    # Optional if you want to store these because easy to recompute
    save(s_est_train, file=paste0(name, "/s_est_train.RData"))
    save(s_est_train_perm, file=paste0(name, "/s_est_train_perm.RData"))
    
    # Save both estimated and true graphs
    save(tree_est, file=paste0(name, "/list_edges.RData"))
    save(graph_est_mat, file=paste0(name, "/graph_est_mat.RData"))
    save(graph_true_mat, file=paste0(name, "/graph_true_mat.RData"))
    
    # Plots
    pdf(paste0(name, "/graph_vis.pdf"))
    plot(L_true)
    dev.off()
    
    plot_hist_marginal(s_est_train_perm, s_true_train, x_train, d, paste0(name, "/est_marginal_tca"))
    plot_hist_bvt(s_est_train_perm, s_true_train, x_train, d, paste0(name, "/est_bvt_tca"))
  }
  
  llhood <- log_kde_tree(tree_est, s_est_train, s_est_test) + log(abs(det(W_est)))
  roc_vals <- compute_roc_vals(graph_est_mat, graph_true_mat)
  
  return(c(llhood, roc_vals))
}

ica_eval <- function(x_train, x_test, W_true, s_true_train, s_true_test, L_true, name, 
                     N, record) {
  d <- ncol(x_train)
  ica_run <- fastICA(x_train, n.comp=d, row.norm=TRUE)
  W_est <- solve(ica_run$A)
  s_est_train <- x_train %*% W_est
  s_est_test <- x_test %*% W_est
  
  # Find the best permutation. 
  perm_vec_train = find_perm(s_est_train, s_true_train, N)
  s_est_train_perm = s_est_train[, perm_vec_train]
  
  if (record) {
    dir.create(name)
    
    save(W_est, file=paste0(name, "/W_est.RData")) 
    save(W_true, file=paste0(name, "/W_true.RData"))
    
    save(perm_vec_train, file=paste0(name, "/perm_vec_train.RData"))
    
    # Optional if you want to store these because easy to recompute
    save(s_est_train, file=paste0(name, "/s_est_train.RData"))
    save(s_est_train_perm, file=paste0(name, "/s_est_train_perm.RData"))
    
    # Plots
    pdf(paste0(name, "/graph_vis.pdf"))
    plot(L_true)
    dev.off()
    
    plot_hist_marginal(s_est_train_perm, s_true_train, x_train, d, paste0(name, "/est_marginal_ica"))
    plot_hist_bvt(s_est_train_perm, s_true_train, x_train, d, paste0(name, "/est_bvt_ica"))
  }
  
  return(log_kde_indep(s_est_train, s_est_test) + 
           log(abs(det(W_est))))
}


###############################################################################
# Plotting functions

multiplot <- function(..., plotlist=NULL, cols) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}


plot_hist_bvt <- function(s_est, s_true, x, d, name) {
  pdf(paste0(name, ".pdf"))
  for (u in 1:(d-1)) {
    for (v in (u+1):d) {
      p1 <- ggplot(as.data.frame(s_true), aes(x = s_true[, u], y = s_true[, v])) + geom_density_2d()
      p2 <- ggplot(as.data.frame(s_est), aes(x = s_est[, u], y = s_est[, v])) + geom_density_2d()
      p3 <- ggplot(as.data.frame(x), aes(x = x[, u], y = x[, v])) + geom_density_2d()
      multiplot(p1, p2, p3, cols = 3)
    }
  }
  dev.off()
}

plot_hist_marginal <- function(s_est, s_true, x, d, name) {
  pdf(paste0(name, ".pdf"))
  par(mfrow = c(1, 3))
  for (j in 1:d) {
    hist(s_true[, j], freq=FALSE, breaks=20, xlab = "", ylab = "Density", main="True source")
    lines(density(s_true[, j]), col="blue", lwd=2)
    hist(s_est[, j], freq=FALSE, breaks=20, xlab="", ylab = "Density", main="Estimated source")
    lines(density(s_est[, j]), col="blue", lwd=2)
    hist(x[, j], breaks = 20, freq=FALSE, xlab="", main="Observed data")
  }
  dev.off()
}



plot_hist_marginal_gg <- function(s_est, s_true, d, name) {
  pdf(paste0(name, ".pdf"))
  p1 <- vector("list", d)
  p2 <- vector("list", d)
  for (i in 1:d) {
    p1[[i]] <- ggplot(as.data.frame(s_true[, i]), aes(x=s_true[, i])) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white", bins=20) +
      geom_density(alpha=.2, fill="blue") + labs(x="", y="Density") + ggtitle("True source")
    p2[[i]] <- ggplot(as.data.frame(s_est[, i]), aes(x=s_est[, i])) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white", bins=20) +
      geom_density(alpha=.2, fill="blue") + labs(x="", y="") + ggtitle("Estimated source")
  }
  multiplot(p1[[1]], p2[[1]], p1[[2]], p2[[2]], p1[[3]], p2[[3]], cols = 2)
  dev.off()
}


