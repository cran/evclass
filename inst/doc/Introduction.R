## ---- message=FALSE------------------------------------------------------
library(evclass)

## ------------------------------------------------------------------------
data(ionosphere)
xtr<-ionosphere$x[1:176,]
ytr<-ionosphere$y[1:176]
xtst<-ionosphere$x[177:351,]
ytst<-ionosphere$y[177:351]


## ------------------------------------------------------------------------
param0<- EkNNinit(xtr,ytr)
options=list(maxiter=300,eta=0.1,gain_min=1e-5,disp=FALSE)
fit<-EkNNfit(xtr,ytr,param=param0,K=5,options=options)

## ------------------------------------------------------------------------
print(fit$err)
table(fit$ypred,ytr)

## ------------------------------------------------------------------------
val<-EkNNval(xtrain=xtr,ytrain=ytr,xtst=xtst,K=5,ytst=ytst,param=fit$param)
print(val$err)
table(val$ypred,ytst)

## ------------------------------------------------------------------------
err<-rep(0,15)
i<-0
for(K in 1:15){
  fit<-EkNNfit(xtr,ytr,K,options=list(maxiter=100,eta=0.1,gain_min=1e-5,disp=FALSE))
  err[K]<-fit$err
}
plot(1:15,err,type="b",xlab='K',ylab='LOO error rate')

## ------------------------------------------------------------------------
fit<-EkNNfit(xtr,ytr,K=8,options=list(maxiter=100,eta=0.1,gain_min=1e-5,disp=FALSE))
val<-EkNNval(xtrain=xtr,ytrain=ytr,xtst=xtst,K=8,ytst=ytst,param=fit$param)
print(val$err)
table(val$ypred,ytst)

## ------------------------------------------------------------------------
data(glass)
xtr<-glass$x[1:89,]
ytr<-glass$y[1:89]
xtst<-glass$x[90:185,]
ytst<-glass$y[90:185]

## ------------------------------------------------------------------------
param0<-proDSinit(xtr,ytr,nproto=7,nprotoPerClass=FALSE,crisp=FALSE)

## ------------------------------------------------------------------------
options<-list(maxiter=500,eta=0.1,gain_min=1e-5,disp=20)
fit<-proDSfit(x=xtr,y=ytr,param=param0,options=options)

## ------------------------------------------------------------------------
val<-proDSval(xtst,fit$param,ytst)
print(val$err)
table(ytst,val$ypred)

