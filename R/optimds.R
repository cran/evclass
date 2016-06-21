# Optimization algorithm for EkNNfit
optimds<- function(x,y,param,knn,K,lambda,options){
  M <- max(y)
  alpha <-param$alpha
  P<-param$gamma

  a <-1.2
  b<-0.8
  c <- 0.5
  mi <- 1e-4
  mx <- 1e6
  pas <- options$eta * rep(1,M)
  it <- 0
  gain <- 1
  costgrad<-gradientds(y,knn,P,alpha,K,lambda)
  Errcou1<-costgrad$cost
  Dp <- costgrad$grad
  Pp <- P
  Errcop <- Errcou1 + 1

  while((gain >= options$gain_min) & (it <= options$maxiter)){
    it = it + 1;
    costgrad<-gradientds(y,knn,P,alpha,K,lambda)
    Errco<-costgrad$cost
    D <- costgrad$grad
    if(is.nan(gain) | is.infinite(gain)) gain <- 1
    if(options$disp) print(c(it,Errco,gain))
    if(Errco > Errcop){
      P <- Pp
      D <- Dp
      pas <- pas * c
      P <- P - pas * D
    }else{
      gain <- .9*gain + .1*abs(Errcop - Errco)
      Pp <- P
      test <-  ((D * Dp) >= 0)
      pas <- ((test * a) + ((!test) * b)) * pas
      pas <- (pas <= mx) * pas + (pas > mx) * mx
      pas <- (pas >= mi) * pas + (pas < mi) * mi
      Dp <- D
      P <- P - pas * D
      Errcop <- Errco
    }
  }
  return(list(param=list(gamma=P,alpha=alpha),cost=Errco))
}




