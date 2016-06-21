classds<-function(param,knn,y,K){
  N<-length(y)
  M=max(y)
  L<-rep(0,N)
  mk <- rbind(matrix(0,M,N),rep(1,N))
  is<-t(knn$nn.index)
  ds<-t(knn$nn.dist)

  for(k in 1:K){
    Is <- is[k,]
    Is = y[Is]
    Tk <- matrix(0,M,N)
    for(j in 1:M){
      pos <- which(Is==j)
      if(length(pos) != 0) Tk[j,pos] <- rep(1,length(pos))
    }
    G <- matrix(param$gamm^2,M,N) * Tk
    gam <- apply(G,2,max)
    s <- param$alpha*exp(-gam *ds[k,])
    m <- rbind(Tk*matrix(s,M,N,byrow=TRUE),1-s)
    mk <- rbind( mk[1:M,]*(m[1:M,]+matrix(m[M+1,],M,N,byrow=TRUE))+
                   m[1:M,]*matrix(mk[M+1,],M,N,byrow=TRUE),mk[M+1,]*m[M+1,])
    Kn <- colSums(mk)
    mk <- mk/ matrix(Kn,M+1,N,byrow=TRUE)
  }
  mk<-t(mk)
  L<-max.col(mk[,1:M])
  err=length(which(L != y))/N
  return(list(m=mk,ypred=L,err=err))
}





