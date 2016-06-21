#' Initialization of parameters for the EkNN classifier
#'
#'\code{EkNNinit} returns initial parameter values for the EkNN classifier.
#'
#'Each parameter \eqn{\gamma_k} is set ot the inverse of the square root of the mean
#'Euclidean distances wihin class k. Note that \eqn{\gamma_k} here is the square root
#'of the \eqn{\gamma_k} as defined in (Zouhal and Denoeux, 1998). By default, parameter alpha is set
#'to 0.95. This value normally does not have to be changed.
#'
#'
#' @param x Input matrix of size n x d, where n is the number of objects and d the number of
#' attributes.
#' @param y Vector of class lables (of length n). May be a factor, or a vector of
#' integers.
#' @param alpha Parameter \eqn{\alpha}.
#'
#' @return A list with two elements:
#'   \describe{
#'   \item{gamma}{Vector of parameters \eqn{\gamma_k}, of length c, the number of classes.}
#'   \item{alpha}{Parameter \eqn{\alpha}, set to 0.95.}
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
#' @importFrom stats dist
#'
#' @seealso \code{\link{EkNNfit}}, \code{\link{EkNNval}}
#'
#' @examples ## Iris dataset
#' data(iris)
#' x<-iris[,1:4]
#' y<-iris[,5]
#' param<-EkNNinit(x,y)
#' param
EkNNinit<-function(x,y,alpha=0.95){
  y<-as.numeric(y)
  M<-max(y)
  gamm<-rep(0,M)
  for(k in 1:M){
    D<-dist(x[y==k,])
    gamm[k]<-1/sqrt(mean(D))
  }
  return(list(gamma=gamm,alpha=alpha))
}
