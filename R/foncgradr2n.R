# Gradient  and cost calculation for the evidential neural network classifier. Called by harris.
foncgradr2n<- function(w,X,t,lambda,mu,optimProto){
  l <- length(w)
  p<-nrow(X)
  N<-ncol(X)
  M<- nrow(t)
  n <- l / (p+2+M)

  W<-matrix(w[1:(n*p)],n,p)
  Gamma <- w[(n*p+1):(n*p+n)]
  BETA <- matrix(w[(n*(p+1)+1):(n*(p+1+M))],n,M)
  alpha <- w[(n*(p+1+M)+1):l]
  BETA2 <- BETA^2
  beta2 <- rowSums(BETA2)
  U <- BETA2 / matrix(beta2,n,M,byrow=TRUE)

  d<-matrix(0,n,N)
  s<-d
  expo<-s
  dU <- matrix(0,n,M)
  Dbeta <- dU
  Ds <- matrix(0,n,N)
  Dgamma <- Ds
  Dalpha <- Ds
  DW<-matrix(0,n,p)
  mk <- rbind(matrix(0,M,N),rep(1,N))
  mm <- mk
  alphap <- 0.99 / (1 + exp(-alpha))

  for(k in 1:n){
    # propagation
    d[k,] <- 0.5*colSums((X - matrix(W[k,],p,N))^2)
    expo[k,] <- exp(- Gamma[k]^2 * d[k,])
    s[k,] <- alphap[k]*expo[k,]
    m = rbind(U[k,] %o% s[k,],1-s[k,])
    mk = rbind( mk[1:M,] * (m[1:M,] + matrix(m[M+1,],M,N,byrow=TRUE)) +
             m[1:M,]*matrix(mk[M+1,],M,N,byrow=TRUE), mk[M+1,]*m[M+1,])
  }
  # normalisation
  K <- colSums(mk)
  K2<-K^2
  mkn<-mk/matrix(K,M+1,N,byrow=TRUE)
  Q <- mkn[1:M,] + lambda * matrix(mkn[M+1,],M,N,byrow=TRUE) - t
  E <- 0.5 * mean(colSums(Q^2)) + mu * sum(alphap)

  dEdm <- matrix(0,M+1,N)
  I <- diag(M)
  for(k in 1:M){
    dEdm[k,] <- colSums( Q * (I[,k]%o%K - mk[1:M,] -
                             lambda*matrix(mk[M+1,],M,N,byrow=TRUE))) / K2
  }
  dEdm[M+1,] <- colSums(Q * (-mk[1:M,] +
                          lambda*matrix(K-mk[M+1,],M,N,byrow=TRUE))) / K2

  mm<-matrix(0,M+1,N)
  for(k in 1:n){
    #  gradient calculation
    m <- rbind(U[k,] %o% s[k,] , 1-s[k,])
    mm[M+1,] <- mk[M+1,] / m[M+1,]
    L <- matrix(mm[M+1,],M,N,byrow=TRUE)
    mm[1:M,] <- (mk[1:M,] - L * m[1:M,]) /
                 (m[1:M,] + matrix(m[M+1,],M,N,byrow=TRUE))
    R <- mm[1:M,] + L
    A <- R * matrix(s[k,],M,N,byrow=TRUE)
    B <- matrix(U[k,],M,N) * R - mm[1:M,]
    dU[k,] <- rowMeans(A * dEdm[1:M,])
    Ds[k,] <- colSums(dEdm[1:M,] * B) - dEdm[M+1,]*mm[M+1,]
    DW[k,] <- - (Ds[k,] * Gamma[k]^2  * s[k,]) %*% (matrix(W[k,],N,p,byrow=TRUE) - t(X))
  }

DW <- as.numeric(optimProto)*DW/N

T <- matrix(beta2,n,M)
Dbeta <- 2*BETA / T^2 * (dU *(T - BETA2) - matrix(rowSums(dU*BETA2),n,M) + dU * BETA2)

Dgamma <- - 2 * rowMeans(Ds * d * s) * Gamma
Dalpha <- (rowMeans(Ds * expo)+ mu) * 0.99 * (1-alphap) * alphap

G = c(as.vector(DW),Dgamma,as.vector(Dbeta),Dalpha)

return(list(E=E,G=G))

}


