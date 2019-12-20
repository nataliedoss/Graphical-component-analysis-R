library(orthopolynom)
library(igraph)
library(parallel)
library(ks)


###############################################################################
# Function to find clusters.
find_components <- function(theta, m1, m2, max_component_size) {
  d = length(theta)
  list_components <- NULL 
  i = 1
  while (i <= d) {
    i_in_list_components = sum(unlist(lapply(list_components, function(x) i %in% x))) > 0
    if (!i_in_list_components) {
      nb_vec = rep(NA, d)
      nb_vec[1:(i-1)] = theta[[i]][1:(i-1)]
      nb_vec[i] = 0
      nb_vec[(i+1):d] = theta[[i]][(i+m1):length(theta[[i]])]
      ind <- rep(NA, d)
      for (j in 1:d) {
        ind[j] <- sum(unlist(lapply(list_components, function(x) j %in% x))) > 0
      }
      nb_vec[ind] = 0
      nbs_ordered = order(abs(nb_vec), decreasing = T)
      nbs = nbs_ordered[nbs_ordered != i & (abs(nb_vec[nbs_ordered]) > .99)][1:max_component_size-1]
      list_components = append(list_components, list(sort(c(i, nbs))))
      i = i + 1
      
    }
    else {
      i = i + 1
    }
  }
  return(list_components)
}

log_kde_component <- function(component, s_train, s_test) {
  density <- kde(x=s_train[, component], eval.points=s_test[, component])
  log_kde = log(10^(-20) + density$estimate)
  return(sum(log_kde))
}


log_kde_indep <- function(s_train, s_test) {
  d <- ncol(s_train)
  list_kde <- lapply(1:d, function(i) log_kde_component(i, s_train, s_test))
  return(Reduce("+", list_kde))
}

log_kde_dep <- function(theta, s_train, s_test, m1, m2, max_component_size) {
  list_components = find_components(theta, m1, m2, max_component_size)
  list_kde <- lapply(list_components, 
                     function(component) log_kde_component(component, s_train, s_test))
  return(Reduce("+", list_kde))
}

log_kde_tree <- function(tree, s_train, s_test) {
  list_kde_edges <- lapply(1:nrow(tree_est), function(i) {
    log_kde_component(as.matrix(tree_est)[i, ], s_train, s_test)})
  kde_edges <- Reduce("+", list_kde_edges)
  
  nodes_tree <- unlist(tree_est)
  list_kde_nodes <- lapply(nodes_tree, function(i) log_kde_component(i, s_train, s_test))
  kde_nodes <- Reduce("+", list_kde_nodes)
  
  return(kde_edges - kde_nodes + log_kde_indep(s_train, s_test))
}






