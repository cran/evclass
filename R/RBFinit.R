#' Initialization of parameters for a Radial Basis Function classifier
#'
#'\code{RBFinit} returns initial parameter values for a Radial Basis Function classifier.
#'
#' The prototypes are initialized by the k-means algorithms. The hidden-to-output weights are initialized
#' by linear regression. The scale parameter for each prototype is computed as the inverse of the square
#' root of the mean squared distances to this prototype. The final number of prototypes may be different
#' from the desired number \code{nproto} depending on the result of the k-means clustering (clusters
#' composed of only one input vector are discarded).
#'
#' @param x Input matrix of size n x d, where n is the number of objects and d the number of
#' attributes.
#' @param y Vector of class labels (of length n). May be a factor, or a vector of
#' integers from 1 to M (number of classes).
#' @param nproto Number of prototypes
#'
#' @return A list with three elements containing the initialized network parameters
#'   \describe{
#'   \item{P}{Matrix of size (R,d), containing the R prototype coordinates.}
#'   \item{Gamma}{Vector of length R, containing the scale parameters.}
#'   \item{W}{Matrix of size (R,M), containing the hidden-to-output weights.}
#'  }
#'
#'@author Thierry Denoeux.
#'
#' @export
#' @importFrom stats lm kmeans
#'
#' @seealso \code{\link{RBFfit}}, \code{\link{RBFval}}
#'
#' @examples ## Glass dataset
#' data(glass)
#' xapp<-glass$x[1:89,]
#' yapp<-glass$y[1:89]
#' param0<-RBFinit(xapp,yapp,nproto=7)
#' param0
RBFinit<-function(x,y,nproto){
  R<-nproto
  X<-x
  y<-as.integer(as.factor(y))
  K <- max(y)
  Id <- diag(K)
  Y<-Id[y,]
  X<-as.matrix(X)
  n <- nrow(X)
  J<-ncol(X)

  clus <- kmeans(X,R)
  ii<-which(clus$size>1)
  R<-length(ii)
  P0<-clus$centers[ii,]

  d2<-matrix(0,n,R)
  for(r in 1:R) d2[,r] <- rowSums((X - matrix(P0[r,],n,J,byrow=TRUE))^2)
  gamma0<-1/sqrt(colMeans(d2))
  s <- exp(-0.5*matrix(gamma0^2,n,R,byrow=TRUE)*d2)
  W0<-matrix(0,R,K)
  for(k in 1:K){
    reg<-lm(Y[,k] ~s+0)
    W0[,k]<-reg$coef
  }

  return(list(P=P0,Gamma=gamma0,W=W0))
}
