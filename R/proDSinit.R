#' Initialization of parameters for the evidential neural network classifier
#'
#'\code{proDSinit} returns initial parameter values for the evidential neural network classifier.
#'
#'The prototypes are initialized by the k-means algorithms. The initial membership values \eqn{u_{ik}} of
#'each prototype \eqn{p_i} to class \eqn{\omega_k} are normally defined as the proportion of training samples
#'from class \eqn{\omega_k} in the neighborhood of prototype \eqn{p_i}. If arguments \code{crisp} and
#'\code{nprotoPerClass} are set to TRUE, the prototypes are assigned to one and only one class.
#'
#' @param x Input matrix of size n x d, where n is the number of objects and d the number of
#' attributes.
#' @param y Vector of class lables (of length n). May be a factor, or a vector of
#' integers from 1 to M (number of classes).
#' @param nproto Number of prototypes.
#' @param nprotoPerClass Boolean. If TRUE, there are \code{nproto} prototypes per class. If
#' FALSE (default), the total number of prototypes is equal to \code{nproto}.
#' @param crisp Boolean. If TRUE, the prototypes have full membership to only one class. (Available only is
#' nprotoPerClass=TRUE).
#'
#' @return A list with four elements containg the initialized network parameters
#'   \describe{
#'   \item{alpha}{Vector of length r, where r is the number of prototypes.}
#'   \item{gamma}{Vector of length r}
#'   \item{beta}{Matrix of size (r,M), where M is the number fo classes.}
#'   \item{W}{Matrix of size (r,d), containing the prototype coordinates.}
#'  }
#'
#'@references T. Denoeux. A neural network classifier based on Dempster-Shafer theory.
#'IEEE Trans. on Systems, Man and Cybernetics A, 30(2):131--150, 2000.
#'
#'Available from \url{https://www.hds.utc.fr/~tdenoeux}.
#'
#'@author Thierry Denoeux.
#'
#' @export
#' @importFrom stats runif kmeans
#'
#' @seealso \code{\link{proDSfit}}, \code{\link{proDSval}}
#'
#' @examples ## Glass dataset
#' data(glass)
#' xapp<-glass$x[1:89,]
#' yapp<-glass$y[1:89]
#' param0<-proDSinit(xapp,yapp,nproto=7)
#' param0
proDSinit<- function(x,y,nproto,nprotoPerClass=FALSE,crisp=FALSE){
  y<-as.integer(as.factor(y))
  x<-as.matrix(x)
  M <- max(y)
  N <- nrow(x)
  Id <- diag(M)
  t<-Id[y,]

  if(nprotoPerClass){
    W0<-NULL
    BETA0<-NULL
    for(i in 1:M){
      ii<-which(y==i)
      clus <- kmeans(x[ii,],nproto)
      W0 <- rbind(W0,clus$centers)
      BETA0 <- rbind(BETA0,t[ii[1:nproto],])
    }
    n<-nproto*M
    if(!crisp) BETA0 <- BETA0 + 0.1
  } else{
    n<-nproto
    BETA0<-matrix(0,n,M)
    clus <- kmeans(x,n)
    class<-clus$cluster
    W0<-clus$centers
    for(i in 1:n){
      ii <- which(class == i)
      if(length(ii)==1){
        BETA0[i,] = t[ii,]
      } else if(length(ii)>0){
        BETA0[i,] = sqrt(colMeans(t[ii,]))
      } else BETA0[i,] <- runif(M)
    } # end for
  } # end if
  alpha0 <- rep(0,n)
  gamma0 <- rep(0.1,n)
  return(list(alpha=alpha0,gamma=gamma0,beta=BETA0,W=W0))

}



