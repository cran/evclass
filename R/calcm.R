#' Determination of optimal coefficients for computing weights of evidence in logistic regression
#'
#'\code{calcAB} transforms coefficients alpha and beta computed by \code{calcm} into weights of
#'evidence, and then into mass and contour (plausibility) functions. These mass functions
#' can be used to express uncertainty about the prediction of logistic regression or multilayer
#' neural network classifiers (See Denoeux, 2019).
#'
#' An error may occur if the absolute values of some coefficients are too high. It is then advised
#' to recompute these coefficients by training the logistic regression or neural network classifier
#' with L2 regularization. With M classes, the output mass functions have 2^M focal sets.
#' Using this function with large M may cause memory issues.
#'
#' @param x Matrix (n,d) of feature values, where d is the number of features, and n is the number
#' of observations. Can be a vector if $d=1$.
#' @param A Vector of length d (for M=2) or matrix of size (d,M) (for M>2) of coefficients alpha.
#' @param B Vector of length d (for M=2) or matrix of size (d,M) (for M>2) of coefficients beta
#'
#' @return A list with six elements:
#'   \describe{
#'   \item{F}{Matrix (2^M,M) of focal sets.}
#'   \item{mass}{Matrix (n,2^M) of mass functions (one in each row).}
#'   \item{pl}{Matrix (n,M) containing the plausibilities of singletons.}
#'   \item{bel}{Matrix (n,M) containing the degrees of belief of singletons.}
#'   \item{prob}{Matrix (n,M) containing the normalized plausibilities of singletons.}
#'   \item{conf}{Vector of length n containing the degrees of conflict.}
#'  }
#'
#'@references
#'T. Denoeux. Logistic Regression, Neural Networks and Dempster-Shafer Theory: a New Perspective.
#'Knowledge-Based Systems, Vol. 176, Pages 54–67, 2019.
#'
#'@author Thierry Denoeux.
#'
#' @export
#' @import ibelief R.utils
#'
#' @seealso \code{\link{calcAB}}
#'
#' @examples ## Example with 2 classes and logistic regression
#' data(ionosphere)
#' x<-ionosphere$x[,-2]
#' y<-ionosphere$y-1
#' fit<-glm(y ~ x,family='binomial')
#' AB<-calcAB(fit$coefficients,colMeans(x))
#' Bel<-calcm(x,AB$A,AB$B)
#' Bel$focal
#' Bel$mass[1:5,]
#' Bel$pl[1:5,]
#' Bel$conf[1:5]
#' ## Example with K>2 classes and multilayer neural network
#' library(nnet)
#' data(glass)
#' K<-max(glass$y)
#' d<-ncol(glass$x)
#' n<-nrow(x)
#' x<-scale(glass$x)
#' y<-as.factor(glass$y)
#' p<-3 # number of hidden units
#' fit<-nnet(y~x,size=p)  # training a neural network with 3 hidden units
#' W1<-matrix(fit$wts[1:(p*(d+1))],d+1,p) # Input-to-hidden weights
#' W2<-matrix(fit$wts[(p*(d+1)+1):(p*(d+1) + K*(p+1))],p+1,K) # hidden-to-output weights
#' a1<-cbind(rep(1,n),x)%*%W1  # hidden unit activations
#' o1<-1/(1+exp(-a1)) # hidden unit outputs
#' AB<-calcAB(W2,colMeans(o1))
#' Bel<-calcm(o1,AB$A,AB$B)
#' Bel$focal
#' Bel$mass[1:5,]
#' Bel$pl[1:5,]
#' Bel$conf[1:5]
calcm<-function(x,A,B){

  calcw1<-function(x,A,B){
    w<-B*matrix(x,nrow(B),ncol(B))+A
    Wp<-matrix(pmax(0,w),nrow(w),ncol(w))
    Wm<-matrix(-pmin(0,w),nrow(w),ncol(w))
    return(c(colSums(Wp),colSums(Wm)))
  }

  if(is.vector(x)) x<-t(as.matrix(x))
  n<-nrow(x)
  J<-ncol(x)
  if(is.vector(B)){
    K<-2
    wj<-x*matrix(B,n,J,byrow=TRUE)+ matrix(A,n,J,byrow=TRUE)
    wp<-rowSums(matrix(pmax(0,wj),n,J))
    wm<-rowSums(matrix(-pmin(0,wj),n,J))
  }else{
    wpm <- t(apply(x,1,calcw1,A,B))
    K<-ncol(B)
  }

  ii<-1:2^K
  N<-length(ii)
  F<-matrix(0,N,K)
  CC<-intToBin(0:(N-1))
  for(i in 1:N) F[i,]<-as.numeric(substring(CC[i],1:K,1:K))
  F<-F[,K:1]

  card<-rowSums(F)
  ii1<-which(card==1)
  ii2<-which(card==K-1)

  if(K==2){
    Wpm<-cbind(rep(1,n),pmax(exp(-wp),0),pmax(exp(-wm),0),rep(1,n))
  }else{
    Wpm<-matrix(1,nrow(wpm),2^K)
    Wpm[,ii1]<-exp(-wpm[,1:K])
    Wpm[,ii2]<-exp(-wpm[,(2*K):(K+1)])
  }
  Wpm<-pmax(Wpm,.Machine$double.xmin)
  mass<-t(apply(Wpm,1,wtom))
  if(any(is.nan(mass))) stop("Infinite masses: numerical problem.")
  conf<-mass[,1]
  mass[,1]<-0
  mass<-mass/matrix(rowSums(mass),nrow(mass),ncol(mass))
  Pl<-t(apply(mass,1,mtopl))
  pl<-Pl[,ii1]
  bel<-mass[,ii1]
  p<-pl/matrix(rowSums(pl),nrow(pl),ncol(pl))

  return(list(focal=F,mass=mass,pl=pl,bel=bel,prob=p,conf=conf))
}
