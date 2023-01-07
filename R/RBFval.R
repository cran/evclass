#' Classification of a test set by a radial basis function classifier
#'
#'\code{RBFval} classifies instances in a test set using a radial basis function classifier. Function
#'\code{\link{calcm}} is called for computing output belief functions. It is recommended to set
#'\code{calc.belief=FALSE} when the number of classes is very large, to avoid memory problems.
#'
#' If class labels for the test set are provided, the test error rate is also returned.
#'
#' @param x Matrix of size n x d, containing the values of the d attributes for the test data.
#' @param param Neural network parameters, as provided by \code{\link{RBFfit}}.
#' @param y Optional vector of class labels for the test data. May be a factor, or a vector of
#' integers from 1 to M (number of classes).
#' @param calc.belief If TRUE (default), output belief functions are calculated.
#'
#' @return A list with four elements:
#'   \describe{
#'   \item{ypred}{Predicted class labels for the test data.}
#'   \item{err}{Test error rate (if the class label of test data has been provided).}
#'   \item{Prob}{Output probabilities.}
#'   \item{Belief}{If \code{calc.belief=TRUE}, output belief function, provided as a list
#'   output by function \code{\link{calcm}}.}
#'  }
#'
#'@references
#'T. Denoeux. Logistic Regression, Neural Networks and Dempster-Shafer Theory: a New Perspective.
#'Knowledge-Based Systems, Vol. 176, Pages 54–67, 2019.
#'
#'Ling Huang, Su Ruan, Pierre Decazes and Thierry Denoeux. Lymphoma segmentation from 3D PET-CT
#' images using a deep evidential network. International Journal of Approximate Reasoning,
#' Vol. 149, Pages 39-60, 2022.
#'
#'@author Thierry Denoeux.
#'
#' @export
#' @import ibelief
#'
#' @seealso \code{\link{RBFinit}}, \code{\link{RBFfit}}, \code{\link{calcm}}
#'
#' @examples ## Glass dataset
#' data(glass)
#' xapp<-glass$x[1:89,]
#' yapp<-glass$y[1:89]
#' xtst<-glass$x[90:185,]
#' ytst<-glass$y[90:185]
#' ## Initialization
#' param0<-RBFinit(xapp,yapp,nproto=7)
#' ## Training
#' fit<-RBFfit(xapp,yapp,param0)
#' ## Test
#' val<-RBFval(xtst,fit$param,ytst)
#' ## Confusion matrix
#' table(ytst,val$ypred)
RBFval<-function(x,param,y=NULL,calc.belief=TRUE){
  X<-x
  K<-ncol(param$W)
  n<-nrow(X)
  J<-ncol(X)
  R<-length(param$Gamma)

  d2<-matrix(0,n,R)
  for(r in 1:R) d2[,r] <- rowSums((X - matrix(param$P[r,],n,J,byrow=TRUE))^2)
  s <- exp(-0.5*matrix(param$Gamma^2,n,R,byrow=TRUE)*d2)

  a<-s %*% param$W
  Pi<-exp(a)
  Pi<-Pi/matrix(rowSums(Pi),n,K)

  B<-param$W-matrix(rowMeans(param$W),R,K)
  A<-matrix(0,R,K)

  # Mass calculation
  if(calc.belief){
    if(K==2){
      Belief<-calcm(s,as.vector(A[,1]),as.vector(B[,1]))
    }  else Belief<-calcm(s,A,B)
  } else Belief<-NULL

  ypred<-max.col(Pi)

  if(!is.null(y)){
    err<-mean((y!=ypred))
  } else err<-NULL

  return(list(ypred=ypred,err=err,Prob=Pi,Belief=Belief))
}
