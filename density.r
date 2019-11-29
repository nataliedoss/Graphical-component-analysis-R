library(orthopolynom)
library(igraph)
library(parallel)
library(ks)
source("suff_stats_unbounded.r")


ss.node<-function(y,m){
  d = length(y)
  out=lapply(1:d,function(i){
    ld2 = t(legendre_Pl_array(m,y[i])[-1])
  })
  return(matrix(unlist(out), nrow = d, ncol = m, byrow = T))
}

ss.edge<-function(y,m,cores,edgelist=NULL){
  d = length(y)
  g=graph.full(d)
  if(is.null(edgelist)) edgelist=get.edgelist(g)
  e=dim(edgelist)[1]
  ret=mclapply(1:e,function(i){
    s=edgelist[i,1]
    t=edgelist[i,2]
    out1 = outer( legendre_Pl_array(m,y[s])[-1], legendre_Pl_array(m,y[t])[-1] )/2
    out2 = outer( legendre_Pl_array(m,y[t])[-1], legendre_Pl_array(m,y[s])[-1] )/2
    return(list(out1 = out1, out2 = out2))
  },mc.cores=cores)
  return(list(array(unlist(lapply(ret,function(x)x$out1)),dim=c(m,m,e)),
              array(unlist(lapply(ret,function(x)x$out2)),dim=c(m,m,e))))
}

ss <- function(y,m1,m2,cores,elist=NULL){
  d = length(y)
	if(is.null(elist)) elist=get.edgelist(graph.full(length(y)))
	g=graph.edgelist(elist,directed=F)
	kvn=ss.node(y,m1)
	kve=ss.edge(y,m2,cores,elist)
	YY=list()
	for(i in 1:d){
		XX=numeric()
		for(j in sort(c(i,neighbors(g,i)))){
			if(i==j){
				XX=c(XX,kvn[i,])
			}else if(i<j){
				e=which(elist[,1]%in% c(i,j) & elist[,2]%in%c(i,j))
				XX= c(XX, kve[[1]][,,e])
			}else{
				e=which(elist[,1]%in% c(i,j) & elist[,2]%in%c(i,j))
				XX= c(XX, kve[[2]][,,e])
			}
		}
		YY[[i]]=XX
	}
	return(YY)
}


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
  return(sum(log(density$estimate)))
}

log_kde_dep <- function(theta, s_train, s_test, m1, m2, max_component_size) {
  list_components = find_components(theta, m1, m2, max_component_size)
  list_kde <- lapply(list_components, 
                     function(component) log_kde_component(component, s_train, s_test))
  return(Reduce("+", list_kde))
}

log_kde_tree <- function(tree, s_train, s_test) {
  list_kde <- lapply(1:nrow(tree), function(i) {
    log_kde_component(as.matrix(tree)[i, ], s_train, s_test)})
  return(Reduce("+", list_kde))
}

log_kde_indep <- function(s_train, s_test) {
  d <- ncol(s_train)
  list_kde <- lapply(1:d, function(i) log_kde_component(i, s_train, s_test))
  return(Reduce("+", list_kde))
}





