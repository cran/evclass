#' Training of a radial basis function classifier
#'
#'\code{RBFfit} performs parameter optimization for a radial basis function (RBF) classifier.
#'
#'The RBF neural network is trained by maximizing the conditional log-likelihood (or, equivalently,
#'by minimizing the cross-entropy loss function). The optimization procedure is the BFGS
#'algorithm implemented in function \code{optim}.
#'
#' @param x Input matrix of size n x d, where n is the number of objects and d the number of
#' attributes.
#' @param y Vector of class labels (of length n). May be a factor, or a vector of
#' integers from 1 to M (number of classes).
#' @param param Initial parameters (see \code{\link{RBFinit}}).
#' @param lambda Regularization hyperparameter (default=0).
#' @param control Parameters passed to function \code{\link{optim}}.
#' @param optimProto Boolean. If TRUE, the prototypes are optimized (default). Otherwise,
#' they are fixed.
#'
#'
#' @return A list with three elements:
#'   \describe{
#'   \item{param}{Optimized network parameters.}
#'   \item{loglik}{Final value of the log-likelihood objective function.}
#'   \item{err}{Training error rate.}
#'  }
#'
#'@author Thierry Denoeux.
#'
#' @export
#' @importFrom stats optim
#'
#' @seealso \code{\link{proDSinit}}, \code{\link{proDSval}}
#'
#' @examples ## Glass dataset
#' data(glass)
#' xapp<-glass$x[1:89,]
#' yapp<-glass$y[1:89]
#' ## Initialization
#' param0<-RBFinit(xapp,yapp,nproto=7)
#' ## Training
#' fit<-RBFfit(xapp,yapp,param0,control=list(fnscale=-1,trace=2))
RBFfit<-function(x,y,param,lambda=0,control=list(fnscale=-1,trace=2,maxit=1000),optimProto=TRUE){
  Popt<-optimProto
  X<-x
  param0<-param
  theta<-c(as.vector(param0$P),param0$Gamma,as.vector(param0$W))
  N<-length(theta)
  if(is.vector(y)){
    y<-as.integer(as.factor(y))
    K <- max(y)
    Id <- diag(K)
    Y<-Id[y,]
  }else{
    Y<-y
    K<-ncol(Y)
  }
  n<-nrow(X)
  J<-ncol(X)
  R<-length(param0$Gamma)

  opt<-optim(par=theta,fn=likRBF,gr=gradRBF,method="BFGS",
             control=control,X,Y,lambda,Popt,param0$P)
  theta<-opt$par
  P<-matrix(theta[1:(R*J)],R,J)
  Gamma<-theta[(R*J+1):(R*J+R)]
  W<-matrix(theta[((J+1)*R+1):N],R,K)

  d2<-matrix(0,n,R)
  for(r in 1:R) d2[,r] <- rowSums((X - matrix(P[r,],n,J,byrow=TRUE))^2)
  s <- exp(-0.5*matrix(Gamma^2,n,R,byrow=TRUE)*d2)

  a<-s %*% W
  Pi<-exp(a)
  Pi<-Pi/matrix(rowSums(Pi),n,K)

  yhat<-max.col(Pi)
  err<-mean(yhat != y)

#  if(K==2){
#    B<-W
#    A<-rep(0,R)
#    else{
#      B<-W-matrix(rowMeans(W),R,K)
#      A<-matrix(0,R,K)
 #   }


  param<-list(P=P,Gamma=Gamma,W=W) #,A=A,B=B)
  return(list(param=param,loglik=opt$value,err=err))
}
