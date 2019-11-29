library(orthopolynom)
library(igraph)
library(parallel)
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

# Function to find triangles. 
# This currently assumes m2 = 1.
# It also assumes d is divisible by 3.
find_triangles <- function(theta, d, m1, m2) {
  vec_triangles <- NULL
  i = 1
  while (i <= d) {
    if (!(i %in% vec_triangles)) {
      nb_vec = rep(NA, d)
      nb_vec[1:(i-1)] = theta[[i]][1:(i-1)]
      nb_vec[i] = 0
      nb_vec[(i+1):d] = theta[[i]][(i+m1):length(theta[[i]])]
      nbs_ordered = order(abs(nb_vec), decreasing = T)
      nbs = nbs_ordered[nbs_ordered != i & !(nbs_ordered %in% vec_triangles)][1:2]
      vec_triangles = append(vec_triangles, sort(c(i, nbs)))
      i = i + 1
    }
    else 
      i = i + 1
  }
  list_triangles <- vector("list", d/3)
  for (i in 1:(d/3)) {
    list_triangles[[i]] = vec_triangles[((3*(i-1))+1):(3*i)]
  }
  return(list_triangles)
}

# Create a theta for a single triangle
# Output has dimension 2*m2 + m1. Input theta has dimension d*((d-1)*m2 + m1)
create_theta_condensed <- function(theta, triangle, m1, m2) {
  theta_triangle = list(theta[[triangle[1]]], theta[[triangle[2]]], theta[[triangle[3]]])
  theta_condensed = vector("list", 3)
  for (i in 1:3) {
    for (j in 1:3) {
      node = triangle[i]
      nb = triangle[j]
      if (i == j) {
        theta_condensed[[i]][((i-1)*m2+1):((i-1)*m2+m1)] = theta_triangle[[i]][((node-1)*m2+1):((node-1)*m2+m1)]
      }
      else if (i < j) {
        theta_condensed[[i]][(m1+(j-2)*m2+1):(m1+(j-2)+m2)] = theta_triangle[[i]][(m1+(nb-2)*m2+1):(m1+(nb-2)+m2)]
      }
      else {
        theta_condensed[[i]][((j-1)*m2+1):(j*m2)] = theta_triangle[[i]][((nb-1)*m2+1):(nb*m2)]
      }
    }
  }
  return(theta_condensed)
}


# Estimate an unnormalized log likelihood for all triangles in a list
est_llk_unnorm_triangles <- function(ss_val, theta, list_triangles, m1, m2) {
  theta_condensed <- lapply(list_triangles, 
                            function(triangle) create_theta_condensed(theta, triangle, m1, m2))
  
  return(lapply(1:length(list_triangles), function(i) ss_val %*% unlist(theta_condensed[[i]])))
}

# Estimate unnormalized log likelihood
est_llk_unnorm <- function(s, theta, m1, m2, cores) {
  d = ncol(s)
  s_ss = vector("list", length = nrow(s)) 
  llk = rep(NA, nrow(s)) 
  for (i in 1:nrow(s)) {
    s_ss[[i]] = ss(s[i, ], m1, m2, cores)
    llk[i] = Reduce('+', lapply(1:d, function(j) sum(s_ss[[i]][[j]] * theta[[j]])))
  }
  return(llk)
}

# Estimate normalized log likelihood
est_llk <- function(s, ss_val, theta, m1, m2, cores) {
  
  # Obtain density on grid so can obtain normalizer
  list_triangles = find_triangles(theta, d, m1, m2)
  llk_unnorm_triangles = est_llk_unnorm_triangles(ss_val, theta, list_triangles, m1, m2)
  
  # Compute normalizer 
  normalizer_vec = unlist(lapply(1:length(list_triangles), 
                                 function(i) log(mean(exp(llk_unnorm_triangles[[i]])))))
  
  # Compute likelihood at s
  llk_unnorm = est_llk_unnorm(s, theta, m1, m2, cores)
  
  # Normalize the log likelihood
  llk = llk_unnorm - sum(normalizer_vec)
  
  return(sum(llk))
}

# Compute marginals and plot
est_density_marginals <- function(ss_val, theta, d, m1, m2, l, indices, name) {
  
  list_triangles = find_triangles(theta, d, m1, m2)
  llk_unnorm_triangles = est_llk_unnorm_triangles(ss_val, theta, list_triangles, m1, m2)
  
  # Compute normalizer 
  normalizer_vec = unlist(lapply(1:length(list_triangles), 
                                 function(i) log(mean(exp(llk_unnorm_triangles[[i]])))))
  
  llk = lapply(1:length(list_triangles), 
                   function(i) llk_unnorm_triangles[[i]] - sum(normalizer_vec[i]))
  
  # Check integrates to 1
  for (i in 1:length(list_triangles)) {
    print(mean(exp(llk[[i]])))
  }
  
  marginals = vector("list", length(list_triangles))
  for (k in 1:length(list_triangles)) {
    marg = matrix(rep(NA, l*3), nrow = 3, ncol = l)
    lle = cbind(llk[[k]], indices)
    for (i in 1:3) {
      for (j in 1:l) {
        lle = lle[order(lle[, i+1]), ] # don't really need for first one
        marg[i, j] = mean(exp(lle[(((j-1)*l^2)+1):(j*l^2), 1]))
      }
    }
    marginals[[k]] = marg
  }
  
  # Plot them
  pdf(paste0(name, "_density_ests_marginals.pdf"))
  for (i in 1:length(list_triangles)) {
    for (j in 1:3) {
      plot(grid, marginals[[i]][j, ])
    }
  }
  dev.off()
  
  return(marginals)
}


###############################################################################
# Create the indices
l = 100
t = 1
grid = seq(-2, 2, length = l)
indices <- matrix(rep(NA, 3*l^3), nrow = l^3, ncol = 3)
for (i in 1:l) {
  for (j in 1:l) {
    for (k in 1:l) {
      indices[t, ] <- c(i, j, k)
      t = t + 1
    }
  }
}

# Compute the sufficient statistics on a grid on [-2,2]^3 and save
# The parameters MUST MATCH what you have in data_gen.r! BAD. REDO LATER.
# COMMENT THIS OUT 
#grid = seq(-2, 2, length = l)
#m1 = 4
#m2 = 1
#cores = 1
#length = 3*(m1 + 2*m2)) # correct full suff stat vector length for triangle
# Compute
#t = 1
#ss_val = matrix(rep(NA, length*(l^3)), nrow = l^3, ncol = length)
#for (i in 1:l) {
#  for (j in 1:l) {
#    for (k in 1:l) {
#      ss_val[t, ] = unlist(ss(c(grid[i], grid[j], grid[k]), m1, m2, cores))
#      t = t + 1
#    }
#  }
#}

#save(ss_val, file="ss_val_smaller.RData")





