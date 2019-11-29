thresh<-function(x,lam){
	if(lam==0) return(x)
	max(1-lam/sum(x^2)^(1/2) , 0)*x
}

admmSel<-function(dat,lam,m1,m2,rho,cores,elist=NULL){
	n=dim(dat)[1]
	print("Calculating Statistics:")
	d=dim(dat)[2]
	if(is.null(elist)) elist = get.edgelist(graph.full(d))
	g=graph.edgelist(elist,directed=F)
	K = Kvec2(dat,m1,m2,cores,elist)
	G = Gamma.ret2(dat,m1,m2,cores,ge=elist)
	dd=dim(G[[1]]$v)[1]

	Ginv=lapply(1:d,function(x) G[[x]]$v %*% ((1/((1/n)*G[[x]]$d^2+rho)) * t(G[[x]]$v)))
	print("Done.")
	znode = matrix(rnorm(d*m1)*1e-10,d,m1)
	zedge = array(rnorm(m2*m2*dim(elist)[1])*1e-10,dim=c(m2,m2,dim(elist)[1]))
	tol=Inf
	x=xold=y=zz=lapply(1:d,function(x)rep(0,dim(G[[x]]$v)[1]))
	admmout=list()
	lam=rev(sort(lam))
	for(lamind in 1:length(lam)){
		print(paste("Run ",lamind," Begun"))
		if(lamind>1){
			x=admmout[[lamind-1]]$x
			y=admmout[[lamind-1]]$y
			zz=admmout[[lamind-1]]$zz
		}
		tol=Inf
		iter=0
		while(tol>1e-3 & iter<500){
		iter=iter+1
		# x step
		xold=x
		zedge_old=zedge
		znode_old=znode
		x=lapply(1:d,function(i){
			Ginv[[i]]%*%(-K[[i]] - y[[i]] + rho*zz[[i]]) 
		})
		#z-y step
		for(i in 1:d){
			ii=sum(neighbors(g,i)<i)
			c= (ii)*m2^2+1
			ind=c: (c+m1-1)
			znode[i,] = (x[[i]][ind] + y[[i]][ind]/rho )
			znode[i,] = thresh(znode[i,],.01*lam[lamind]/rho)
			y[[i]][ind] = y[[i]][ind] + rho*( x[[i]][ind] - znode[i,])
			zz[[i]][ind] = znode[i,]
			for(j in neighbors(g,i)){
		
						jj = which(neighbors(g,i)==j)
						ii = which(neighbors(g,j)==i)

					if(i<j){
						c1 = m1+(jj-1)*m2^2+1
						c2 = (ii-1)*m2^2+1
					}else{
						c2 = m1+(ii-1)*m2^2+1
						c1 = (jj-1)*m2^2+1
					}

					e= which(elist[,1]%in%c(i,j) & elist[,2]%in%c(i,j))
					ind1 =c1:(c1+m2^2-1)
					ind2 =c2:(c2+m2^2-1)
					zedge[,,e] = .5* (x[[i]][ind1] + x[[j]][ind2] +
									y[[i]][ind1]/rho + y[[j]][ind2]/rho) 
					zedge[,,e]=matrix( thresh(zedge[,,e],lam[lamind]/rho) , m2 , m2)

					y[[i]][ind1] = y[[i]][ind1] + rho*(x[[i]][ind1] - c(zedge[,,e]))
					y[[j]][ind2] = y[[j]][ind2] + rho*(x[[j]][ind2] - c(zedge[,,e]))
					zz[[i]][ind1] = c(zedge[,,e])
					zz[[j]][ind2] = c(zedge[,,e])
				}
			}
	
		tol = ( sum(abs(zedge-zedge_old)) + sum(abs(znode-znode_old)) ) 
			print(paste("tol:",tol))
		}
		admmout[[lamind]]= list(x=x,y=y,zz=zz,znode=znode,zedge=zedge)
	}
	return(admmout)
}


