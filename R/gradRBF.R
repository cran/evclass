# Gradient of the log-likelihood objective function. Used by RBFfit.
gradRBF<-function(theta,X,Y,lambda,Popt,P0){
  K<-ncol(Y)
  n<-nrow(X)
  J<-ncol(X)
  N<-length(theta)
  R<-N/(J+1+K)
  P<-matrix(theta[1:(R*J)],R,J)
  Gamma<-theta[(R*J+1):(R*J+R)]
  W<-matrix(theta[((J+1)*R+1):N],R,K)

  # Propagation
  d2<-matrix(0,n,R)
  for(r in 1:R) d2[,r] <- rowSums((X - matrix(P[r,],n,J,byrow=TRUE))^2)
  s <- exp(-0.5*matrix(Gamma^2,n,R,byrow=TRUE)*d2)

  a<-s %*% W
  Pi<-exp(a)
  Pi<-Pi/matrix(rowSums(Pi),n,K)
  SYPi<-rowSums(Y*Pi)

  # Gradient calculation
  Delta1<-matrix(0,n,K)
  delta<-diag(K)
  for(k in 1:K) Delta1[,k]<-rowSums(Y*Pi*(matrix(delta[k,],n,K,byrow=TRUE)-matrix(Pi[,k],n,K)))/SYPi
  DlDW<-t(s)%*%Delta1/n - 2*lambda*W

  Delta2<-Delta1%*%t(W)
  DlDGamma<-colMeans(-Delta2*matrix(Gamma,n,R,byrow=TRUE)*d2*s)

  DlDP<-matrix(0,R,J)
  if(Popt){
    for(r in 1:R){
      DlDP[r,]<-colMeans(matrix(Delta2[,r]*rep(Gamma[r]^2)*s[,r],n,J)*(X-matrix(P[r,],n,J,byrow=TRUE)))#-2*mu*(P[r,]-P0[r,])
    }
  }

  grad<-c(as.vector(DlDP),DlDGamma,as.vector(DlDW))

  return(grad)
}
