#' Classification of a test set by the EkNN classifier
#'
#'\code{EkNNval} classifies instances in a test set using the EkNN classifier.
#'
#' If class labels for the test set are provided, the test error rate is also returned.
#' If parameters are not supplied, they are given default values by \code{\link{EkNNinit}}.
#'
#'
#' @param xtrain Matrix of size ntrain x d, containing the values of the d attributes for the
#' training data.
#' @param ytrain Vector of class labels for the training data (of length ntrain). May
#' be a factor, or a vector of integers.
#' @param xtst Matrix of size ntst x d, containing the values of the d attributes for the
#' test data.
#' @param K Number of neighbors.
#' @param ytst Vector of class labels for the test data (optional). May
#' be a factor, or a vector of integers.
#' @param param Parameters, as returned by \code{\link{EkNNfit}}.
#'
#' @return A list with three elements:
#'   \describe{
#'   \item{m}{Predicted mass functions for the test data. The first M columns correspond
#'   to the mass assigned to each class. The last column corresponds to the mass
#'   assigned to the whole set of classes.}
#'   \item{ypred}{Predicted class labels for the test data.}
#'   \item{err}{Test error rate.}
#'  }
#'
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
#' @seealso \code{\link{EkNNinit}}, \code{\link{EkNNfit}}
#'
#' @examples ## Iris dataset
#' data(iris)
#' train<-sample(150,100)
#' xtrain<-iris[train,1:4]
#' ytrain<-iris[train,5]
#' xtst<-iris[-train,1:4]
#' ytst<-iris[-train,5]
#' K<-5
#' fit<-EkNNfit(xtrain,ytrain,K)
#' test<-EkNNval(xtrain,ytrain,xtst,K,ytst,fit$param)
EkNNval <- function(xtrain,ytrain,xtst,K,ytst=NULL,param=NULL){

  ytrain<-as.numeric(ytrain)
  if(!is.null(ytst)) ytst<-as.numeric(ytst)

  if(is.null(param)) param<-EkNNinit(xtrain,ytrain)

  Napp<-nrow(xtrain)
  M<-max(ytrain)
  N<-nrow(xtst)

  knn<-get.knnx(xtrain, xtst, k=K)
  knn$nn.dist<-knn$nn.dist^2
  is<-t(knn$nn.index)
  ds<-t(knn$nn.dist)

  m = rbind(matrix(0,M,N),rep(1,N))


  for(i in 1:N){
    for(j in 1:K){
      m1 <- rep(0,M+1)
      m1[ytrain[is[j,i]]] <- param$alpha*exp(-param$gamma[ytrain[is[j,i]]]^2*ds[j,i])
      m1[M+1] <- 1 - m1[ytrain[is[j,i]]]
      m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
      m[M+1,i] <- m1[M+1] * m[M+1,i]
      m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
    }
  }
  m<-t(m)
  ypred<-max.col(m[,1:M])
  if(!is.null(ytst)) err<-length(which(ypred != ytst))/N else err<-NULL

  return(list(m=m,ypred=ypred,err=err))

}



