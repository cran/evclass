# Log-likelihood objective function. Used by RBFfit.
likRBF<-function(theta,X,Y,lambda,Popt,P0){
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
  return(mean(log(rowSums(Y*Pi))) - lambda * sum(W^2)) #- mu*sum((P-P0)^2))
}
