#' Determination of optimal coefficients for computing weights of evidence in logistic regression
#'
#'\code{calcAB} computes optimal coefficients alpha and beta needed to transform coefficients
#' from logistic regression (or connections weights between the last hidden layer and the output
#' layer of multilayer neural networks) into weights of evidence. These weights of evidence
#' can then be used to express the outputs of logistic regression or multilayer neural networks
#' as "latent" mass functions.
#'
#' @param W Vector of coefficients of length (d+1), where d is the number of features, in the
#' case of M=2 classes, or (d+1,M) matrix of coefficients (or connection weights) in the case
#' of M>2 classes.
#' @param mu Optional vector containing the means of the d features.
#'
#' @return A list with two elements:
#'   \describe{
#'   \item{A}{Vector of length d (M=2) or matrix of size (d,M) (for M>2) of coefficients alpha.}
#'   \item{B}{Vector of length d (M=2) or matrix of size (d,M) (for M>2) of coefficients beta.}
#'  }
#'
#'@references
#'T. Denoeux. Logistic Regression, Neural Networks and Dempster-Shafer Theory: a New Perspective.
#'Knowledge-Based Systems, Vol. 176, Pages 54–67, 2019.
#'
#'@author Thierry Denoeux.
#'
#' @export
#'
#' @seealso \code{\link{calcm}}
#'
#' @examples ## Example with 2 classes and logistic regression
#' data(ionosphere)
#' x<-ionosphere$x[,-2]
#' y<-ionosphere$y-1
#' fit<-glm(y ~ x,family='binomial')
#' AB<-calcAB(fit$coefficients,colMeans(x))
#' AB
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
#' AB
calcAB<-function(W,mu=NULL){
  if(is.vector(W)){ # W is a vector of length J+1
    K<-2
    J<-length(W)-1
    if(is.null(mu)) mu<-rep(0,J)
    B<-W[2:(J+1)]
    A<-rep(W[1]/J,J)+rep(mean(B*mu),J)-B*mu
  }else{
    # W is a (J+1)xK weight matrix
    K<- ncol(W)
    J<-nrow(W)-1
    if(is.null(mu)) mu<-rep(0,J)
    W<-W-matrix(rowMeans(W),J+1,K)
    B<-as.matrix(W[2:(J+1),])
    if(J==1) B<-t(B)
    MU<-matrix(mu,J,K)
    BMU<-B*MU
    A<-matrix(W[1,]/J,J,K,byrow=TRUE)+matrix(colMeans(BMU),J,K,byrow=TRUE)-BMU
  }
  return(list(A=A,B=B))
}
