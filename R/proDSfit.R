#' Training of the evidential neural network classifier
#'
#'\code{proDSfit} performs parameter optimization for the evidential neural network classifier.
#'
#'If \code{optimProto=TRUE} (default), the prototypes are optimized. Otherwise, they are fixed to
#'their initial value.
#'
#' @param x Input matrix of size n x d, where n is the number of objects and d the number of
#' attributes.
#' @param y Vector of class lables (of length n). May be a factor, or a vector of
#' integers from 1 to M (number of classes).
#' @param param Initial parameters (see \code{link{proDSinit}}).
#' @param lambda Parameter of the cost function. If \code{lambda=1}, the
#' cost function measures the error between the plausibilities and the 0-1 target values.
#' If \code{lambda=1/M}, where M is the number of classes (default), the piginistic probabilities
#' are considered in the cost function. If \code{lambda=0}, the beliefs are used.
#' @param mu Regularization hyperparameter (default=0).
#' @param optimProto Boolean. If TRUE, the prototypes are optimized (default). Otherwise, they are fixed.
#' @param options A list of parameters for the optimization algorithm: maxiter
#' (maximum number of iterations), eta (initial step of gradient variation),
#' gain_min (minimum gain in the optimisation loop), disp (integer; if >0, intermediate
#' results are displayed every disp iterations).
#'
#'
#' @return A list with three elements:
#'   \describe{
#'   \item{param}{Optimized network parameters.}
#'   \item{cost}{Final value of the cost function.}
#'   \item{err}{Training error rate.}
#'  }
#'
#'@references T. Denoeux. A neural network classifier based on Dempster-Shafer theory.
#'IEEE Trans. on Systems, Man and Cybernetics A, 30(2):131--150, 2000.
#'
#'@author Thierry Denoeux.
#'
#' @export
#'
#' @seealso \code{\link{proDSinit}}, \code{\link{proDSval}}
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
proDSfit <- function(x,y,param,lambda=1/max(as.numeric(y)),mu=0,optimProto=TRUE,
                     options=list(maxiter=500,eta=0.1,gain_min=1e-4,disp=10)){
  x<-as.matrix(x)
  y<-as.integer(as.factor(y))
  M<-max(y)
  n<-nrow(param$W)
  p<-ncol(param$W)
  Id <- diag(M)
  t<-Id[y,]
  N<-nrow(x)
  w <- c(as.vector(param$W),param$gamma,as.vector(param$beta),param$alpha)
  l <- length(w)

  opt<-harris(w,t(x),t(t),lambda,mu,optimProto,options)

  W<-matrix(opt$w[1:(n*p)],n,p)
  Gamma <- opt$w[(n*p+1):(n*p+n)]
  BETA <- matrix(opt$w[(n*(p+1)+1):(n*(p+1+M))],n,M)
  alpha <- opt$w[(n*(p+1+M)+1):l]
  param<-list(alpha=alpha,gamma=Gamma,beta=BETA,W=W)

  propag<-proDSval(x,param,y)

  return(list(param=param,cost=opt$E,err=propag$err))
}


