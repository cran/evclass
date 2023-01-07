#' Classification of a test set by the evidential neural network classifier
#'
#'\code{proDSval} classifies instances in a test set using the evidential neural network classifier.
#'
#' If class labels for the test set are provided, the test error rate is also returned.
#'
#' @param x Matrix of size n x d, containing the values of the d attributes for the test data.
#' @param param Neural network parameters, as provided by \code{\link{proDSfit}}.
#' @param y Optional vector of class labels for the test data. May be a factor, or a vector of
#' integers from 1 to M (number of classes).
#'
#' @return A list with three elements:
#'   \describe{
#'   \item{m}{Predicted mass functions for the test data. The first M columns correspond
#'   to the mass assigned to each class. The last column corresponds to the mass
#'   assigned to the whole set of classes.}
#'   \item{ypred}{Predicted class labels for the test data.}
#'   \item{err}{Test error rate (if the class label of test data has been provided).}
#'  }
#'
#'@references
#'T. Denoeux. A neural network classifier based on Dempster-Shafer theory.
#'IEEE Trans. on Systems, Man and Cybernetics A, 30(2):131--150, 2000.
#'
#'@author Thierry Denoeux.
#'
#' @export
#'
#' @seealso \code{\link{proDSinit}}, \code{\link{proDSfit}}
#'
#' @examples ## Glass dataset
#' data(glass)
#' xapp<-glass$x[1:89,]
#' yapp<-glass$y[1:89]
#' xtst<-glass$x[90:185,]
#' ytst<-glass$y[90:185]
#' ## Initialization
#' param0<-proDSinit(xapp,yapp,nproto=7)
#' ## Training
#' fit<-proDSfit(xapp,yapp,param0)
#' ## Test
#' val<-proDSval(xtst,fit$param,ytst)
#' ## Confusion matrix
#' table(ytst,val$ypred)
proDSval<-function(x,param,y=NULL){
  n<-nrow(param$W)
  p<-ncol(param$W)
  x<-as.matrix(x)
  y<-as.integer(as.factor(y))
  N <- nrow(x)
  M <- ncol(param$beta)
  x<-t(x)

  BETA2<-param$beta^2
  beta2 <- rowSums(BETA2)
  U <- BETA2 / matrix(beta2,n,M)
  alphap<-0.99 / (1 + exp(-param$alpha))

  mk <- rbind(matrix(0,M,N),rep(1,N))
  for(k in 1:n){
    d <- 0.5*colSums((x - matrix(param$W[k,],p,N))^2)
    s <- alphap[k]*exp(- param$gamma[k]^2 * d)
    m = rbind(U[k,] %o% s,1-s)
    mk = rbind( mk[1:M,] * (m[1:M,] + matrix(m[M+1,],M,N,byrow=TRUE)) +
                  m[1:M,]*matrix(mk[M+1,],M,N,byrow=TRUE), mk[M+1,]*m[M+1,])
  }
  # normalisation
  K <- colSums(mk)
  mk<-t(mk/matrix(K,M+1,N,byrow=TRUE))
  ypred<-max.col(mk[,1:M])

  if(!is.null(y)){
    err<-length(which(y!=ypred))/N
  } else err<-NULL

  return(list(m=mk,ypred=ypred,err=err))
}




