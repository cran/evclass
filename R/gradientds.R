# Gradient and cost function for EkNN, used by optimds.R
gradientds<-function(y,knn,P,alpha,K,lambda){
  N<-length(y)
  M <- max(y)
  Ident <- diag(M)
  T <- Ident[,y]

  gama <- P
  Dgama = rep(0,M)
  Ds <- rep(0,N)
  Dsgama <- matrix(0,M,N)
  mk<-rbind(matrix(0,M,N),rep(1,N))
  mm <- mk

  s <- matrix(0,K,N)

  is<-t(knn$nn.index)
  ds<-t(knn$nn.dist)

  logK=0
  for(k in 1:K){
    Is <- is[k,]
    Is = y[Is]
    Tk <- matrix(0,M,N)
    for(j in 1:M){
      pos <- which(Is==j)
      if(length(pos) != 0) Tk[j,pos] <- rep(1,length(pos))
    }
    G <- matrix(gama^2,M,N) * Tk
    gam <- apply(G,2,max)
    s[k,] <- alpha*exp(-gam *ds[k,])
    m <- rbind(Tk*matrix(s[k,],M,N,byrow=TRUE),1-s[k,])
    mk <- rbind( mk[1:M,]*(m[1:M,]+matrix(m[M+1,],M,N,byrow=TRUE))+
          m[1:M,]*matrix(mk[M+1,],M,N,byrow=TRUE),mk[M+1,]*m[M+1,])
  }

  # Normalization
  Kn <- colSums(mk)
  mkn <- mk/ matrix(Kn,M+1,N,byrow=TRUE)
  Q <- mkn[1:M,]+lambda*matrix(mkn[M+1,],M,N,byrow=TRUE) - T

  ERR <- 0.5*mean(colSums(Q^2))

  #		gradient

  Dsm <- matrix(0,M+1,N)

  for(k in 1:K){
    Is <- is[k,]
    Is <- y[Is]
    Tk <- matrix(0,M,N)
    for(j in 1:M){
      pos <- which(Is==j)
      if(length(pos) != 0) Tk[j,pos] <- rep(1,length(pos))
    }
    m <- rbind(Tk*matrix(s[k,],M,N,byrow=TRUE),1-s[k,])
    if(length(which(m[M+1,] == 0)) > 1){
      Dsgama <- matrix(1e10,M,N)
      k <- K+1
    } else{
      mm[M+1,] <- mk[M+1,]/m[M+1,]
      H = matrix(mm[M+1,],M,N,byrow=TRUE)
      mm[1:M,] <- (mk[1:M,] - H*m[1:M,])/(m[1:M,]+ matrix(m[M+1,],M,N,byrow=TRUE))
      v<-  (mm[1:M,] + matrix(mm[M+1,],M,N,byrow=TRUE)) * Tk - mm[1:M,]
      DsK <- colSums(v) - mm[M+1,]
      Dsm[1:M,] <- (matrix(Kn,M,N,byrow=TRUE)*v -
        mk[1:M,]*matrix(DsK,M,N,byrow=TRUE) ) / matrix(Kn^2,M,N,byrow=TRUE)
      Dsm[M+1,] <- (-Kn*mm[M+1,] - mk[M+1,]*DsK)/(Kn^2)
      Ds <- colSums(Q*(Dsm[1:M,]+lambda*matrix(Dsm[M+1,],M,N,byrow=TRUE)))
      Dsgama <- Dsgama -2*gama%*%t(Ds*ds[k,]*s[k,]) *Tk
    }
  }

  return(list(cost=ERR,grad=rowMeans(Dsgama)))

}



