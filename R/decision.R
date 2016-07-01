#' Decision rules for evidential classifiers
#'
#'\code{decision} returns decisions from a loss matrix and mass functions computed
#'by an evidential classifier.
#'
#'This function implements the decision rules described in Denoeux (1997), with an
#'arbitrary loss function. The decision rules are the minimization of the lower,
#'upper or pignistic expectation, and Jaffray's decision rule based on minimizing a
#'convex combination of the lower and upper expectations. The function also handles
#'the case where there is an "unknown" class, in addition to the classes represented
#'in the training set.
#'
#' @param m Matrix of masses for n test cases. Each row is a mass function. The first M columns
#' correspond to the mass assigned to each of the M classes. The last column
#' corresponds to the mass assigned to the whole set of classes.
#' @param L The loss matrix of dimension (na,M) or (na,M+1), where na is the set
#' of actions. L[k,j] is the loss incurred of action j is chosen and the true class
#' if \eqn{\omega_k}. If L has M+1 rows, the last row corresponds to the unknown
#' class.
#' @param rule Decision rule to be used. Must be one of these: 'upper' (upper
#' expectation), 'lower' (lower expectations), 'pignistic' (pignistic expectation),
#' 'hurwicz' (weighted sum of the lower and upper expectations).
#' @param rho Parameter between 0 and 1. Used only is rule='rho'.
#'
#' @return A n-vector with the decisions (integers between 1 and na).
#'
#'@references
#'T. Denoeux. Analysis of evidence-theoretic decision rules for pattern
#'classification. Pattern Recognition, 30(7):1095--1107, 1997.
#'
#'Available from \url{https://www.hds.utc.fr/~tdenoeux}.
#'
#'@author Thierry Denoeux.
#'
#' @export
#'
#' @seealso \code{\link{EkNNval}}, \code{\link{proDSval}}
#'
#' @examples ## Example with M=2 classes
#' m<-matrix(c(0.9,0.1,0,0.4,0.6,0,0.1,0.1,0.8),3,3,byrow=TRUE)
#' ## Loss matrix with na=4 acts: assignment to class 1, assignment to class2,
#' # rejection, and assignment to the unknown class.
#' L<-matrix(c(0,1,1,1,0,1,0.2,0.2,0.2,0.25,0.25,0),3,4)
#' d<-decision(m,L,'upper') ## instances 2 and 3 are rejected
#' d<-decision(m,L,'lower') ## instance 2 is rejected, instance 3 is
#' # assigned to the unknown class
#'
decision<-function(m,L=1-diag(ncol(m)-1),rule=c('upper','lower','pignistic','hurwicz'),
                   rho=0.5){
  M<-ncol(m)-1
  n<-nrow(m)
  if(nrow(L)==(M+1)) m<-cbind(m[,1:M],rep(0,n),m[,M+1]) # unknown class
  if(rule=='upper'){
    L1<-apply(L,2,max)
  }else if(rule=='lower'){
    L1<-apply(L,2,min)
  }else if(rule=='pignistic'){
    L1<-colMeans(L)
  }else if(rule=='hurwicz'){
    L1<-rho*apply(L,2,min)+ (1-rho)*apply(L,2,max)
  }
  C<-m %*% rbind(L,L1)
  return(max.col(-C))
}
