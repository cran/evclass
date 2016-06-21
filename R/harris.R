# optimization algorithm used by proDSfit
harris <-function(x0,z,t,lambda,mu,optimProto,options){
  pas <- rep(options$eta,length(x0))
  a <- 1.2
  b <- 0.8
  c <- 0.5
  ovf <- 1e4
  unf <- 1e-6
  x<-x0

  it <- 0
  gain<-1
  yg<- foncgradr2n(x,z,t,lambda,mu,optimProto)
  yp<-yg$E
  gp<-yg$G


  while((gain >= options$gain_min) & (it <= options$maxiter)){
    it<-it+1;
    yg<-foncgradr2n(x,z,t,lambda,mu,optimProto);
    y<-yg$E
    g<-yg$G
      if(options$disp >0){
        if(it-floor(it/options$disp)*options$disp ==1) print(c(it, y, 10*gain))
      }
    if(y > yp){
      x <- xp
      g <- gp
      pas <- pas * c
      x <- x - pas * g
    } else{
      gain<-0.9*gain+0.1*abs(yp-y)
      xp <- x
      test <- (g * gp) >= 0
      pas <- ((test * a) + ((!test) * b)) * pas
      pas <- (pas<=ovf) * pas + (pas>ovf) * ovf
      pas <- (pas>=unf) * pas + (pas<unf) * unf
      gp <- g
      x <- x - pas * g
      yp <- y
    } # end if then else
  } # end while

  return(list(w=x,E=y))
}







