#' Training of the EkNN classifier
#'
#'\code{EkNNfit} optimizes the parameters of the EkNN classifier.
#'
#'If the argument \code{param} is not supplied, the function \code{\link{EkNNinit}} is called.
#'
#' @param x Input matrix of size n x d, where n is the number of objects and d the number of
#' attributes.
#' @param y Vector of class labels (of length n). May be a factor, or a vector of
#' integers.
#' @param K Number of neighbors.
#' @param param Initial parameters (default: NULL).
#' @param alpha Parameter \eqn{\alpha} (default: 0.95)
#' @param lambda Parameter of the cost function. If \code{lambda=1}, the
#' cost function measures the error between the plausibilities and the 0-1 target values.
#' If \code{lambda=1/M}, where M is the number of classes (default), the piginistic probabilities
#' are considered in the cost function. If \code{lambda=0}, the beliefs are used.
#' @param optimize Boolean. If TRUE (default), the parameters are optimized.
#' @param options A list of parameters for the optimization algorithm: maxiter
#' (maximum number of iterations), eta (initial step of gradient variation),
#' gain_min (minimum gain in the optimisation loop), disp (Boolean; if TRUE, intermediate
#' results are displayed during the optimization).
#'
#' @return A list with five elements:
#'   \describe{
#'   \item{param}{The optimized parameters.}
#'   \item{cost}{Final value of the cost function.}
#'   \item{err}{Leave-one-out error rate.}
#'   \item{ypred}{Leave-one-out predicted class labels.}
#'   \item{m}{Leave-one-out predicted mass functions. The first M columns correspond
#'   to the mass assigned to each class. The last column corresponds to the mass
#'   assigned to the whole set of classes.}
#'  }
#'
#'@references T. Denoeux. A k-nearest neighbor classification rule based on Dempster-Shafer
#'theory. IEEE Transactions on Systems, Man and Cybernetics, 25(05):804--813, 1995.
#'
#' L. M. Zouhal and T. Denoeux. An evidence-theoretic k-NN rule with parameter
#' optimization. IEEE Transactions on Systems, Man and Cybernetics Part C,
#' 28(2):263--271,1998.
#'
#'Available from \url{https://www.hds.utc.fr/~tdenoeux}.
#'
#'@author Thierry Denoeux.
#'
#' @export
#' @import FNN
#'
#' @seealso \code{\link{EkNNinit}}, \code{\link{EkNNval}}
#'
#' @examples ## Iris dataset
#' data(iris)
#' x<-iris[,1:4]
#' y<-iris[,5]
#' fit<-EkNNfit(x,y,K=5)
EkNNfit<-function(x,y,K,param=NULL,alpha=0.95,lambda=1/max(as.numeric(y)),optimize=TRUE,
                  options=list(maxiter=300,eta=0.1,gain_min=1e-6,disp=TRUE)){
  y<-as.integer(y)
  x<-as.matrix(x)
  if(is.null(param)) param<-EkNNinit(x,y,alpha)
  knn<-get.knn(x,k=K)
  knn$nn.dist<-knn$nn.dist^2
  if(optimize) opt<-optimds(x,y,param,knn,K,lambda,options)
  class <- classds(opt$param,knn,y,K)
  return(list(param=opt$param,cost=opt$cost,err=class$err,ypred=class$ypred,m=class$m))
}
