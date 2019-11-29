library(orthopolynom)
library(igraph)
library(parallel)



legendre_Pl_array<-function(m,y){
  pf<-polynomial.functions(hermite.h.polynomials(m,normalized=T))
  lv<-sapply(X=1:(m+1),FUN=function(x){
    pf[[x]](y) 
  })
  return(t(lv))
}

len_deriv<-function(m,y){
	out=matrix(0,m+1,length(y))
	out[1,]=0
	lpa<-legendre_Pl_array(m,y)[,,drop=F]
	for(i in 2:(m+1)){
		out[i,]<- 2*(i-1)*lpa[i-1,] 
	}
	return(out)
} 

len_deriv2<-function(m,y){
	out=matrix(0,m+1,length(y))
	out[1,]=0
	out[2,]=0
	if(m==1) return(out)
	lpa<-legendre_Pl_array(m,y)[,,drop=F]
	for(i in 3:(m+1)){
		k=i-1
		out[i,] = 4*(i-1)*(i-2)*lpa[i-2,]
	}
	return(out)
}


Kvec.node<-function(dat,m){
	d=dim(dat)[2]
	out=matrix(0,d,m)
	out=lapply(1:d,function(i){
		ld2 = t(len_deriv2(m,dat[,i])[-1,,drop=F])
		rowMeans(t(ld2))
	})
	t(simplify2array(out))
}

Kvec.edge<-function(dat,m,cores,edgelist=NULL){
	d=dim(dat)[2]
	g=graph.full(d)
	if(is.null(edgelist)) edgelist=get.edgelist(g)
	e=dim(edgelist)[1]
	ret=mclapply(1:e,function(i){
		s=edgelist[i,1]
		t=edgelist[i,2]
		out1 = crossprod( t(len_deriv2(m,dat[,s])[-1,,drop=F]), t(legendre_Pl_array(m,dat[,t])[-1,,drop=F]) )
		out2 = legendre_Pl_array(m,dat[,s])[-1,,drop=F] %*% (t(len_deriv2(m,dat[,t])[-1,,drop=F]) )
		return(list(out1=out1/dim(dat)[1],out2=out2/dim(dat)[1]))
	},mc.cores=cores)
	return(list(array(unlist(lapply(ret,function(x)x$out1)),dim=c(m,m,e)),
		array(unlist(lapply(ret,function(x)x$out2)),dim=c(m,m,e))))
}


Gamma.ret2<-function(dat,m1,m2,cores=1,ge=NULL){
	d=dim(dat)[2]
	if(is.null(ge))ge=get.edgelist(graph.full(d))	
	g=graph.edgelist(ge,directed=F)
	nedge=dim(ge)[1]
	out=mclapply(1:d,function(i){
		XX=numeric()
		for(j in sort(c(i,neighbors(g,i)))){	
				if(i==j){
					ld=len_deriv(m1,dat[,i])[-1,,drop=F]
					XX=cbind(XX, t(ld) )
				}else if(i<j){
					e=which(ge[,1]%in%c(i,j) & ge[,2]%in%c(i,j))
					ld=len_deriv(m2,dat[,i])[-1,,drop=F]
					lpa=legendre_Pl_array(m2,dat[,j])[-1,,drop=F]
					for(l in 1:m2){
						XX=cbind(XX, t(ld) * lpa[l,])
					}
				}else{
					e=which(ge[,1]%in%c(i,j) & ge[,2]%in%c(i,j))
					ld=len_deriv(m2,dat[,i])[-1,,drop=F]
					lpa=legendre_Pl_array(m2,dat[,j])[-1,,drop=F]
					for(l in 1:m2){
						XX=cbind(XX, t(lpa) * ld[l,])
					}
				}
		}		
		return(svd(XX,nu=0))
 	},mc.cores=cores)
	return(out)

}

Kvec2<-function(dat,m1,m2,cores,elist=NULL){
	if(is.null(elist)) elist=get.edgelist(graph.full(dim(dat)[2]))
	g=graph.edgelist(elist,directed=F)
	kvn=Kvec.node(dat,m1)
	kve=Kvec.edge(dat,m2,cores,elist)
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







