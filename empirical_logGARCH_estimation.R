
library(Rsolnp)


1997

n<-length(y)

 
z[, 2:11] <- scale(z[, 2:11])
z_raw = z

z<- cbind(z_raw[,c(1, 2, 7,8)],M%*%z_raw[,c(5, 8)])

dim(z)


ml_splogGARCHX_sato_exo<-function(par,y,z){
  rho<-par[1]
  lambda <- par[2]

  
  n<-length(y)
  # logM<--sum(log(rep(1,n)-rho*Meig))
  logR <- sum(log(rep(1,n)-lambda*Meig))
  logR_12 <- sum(log(rep(1,n)-(lambda + rho)*Meig))
  
  logy=log(y^2)
  R <- diag(n) - lambda * M
  R_ <- solve(R)
  
  R_12 <- diag(n) - (lambda + rho) * M
  
   
  gammaa <- solve(t(z)%*%t(R_)%*%R_%*%z)%*%t(z)%*%t(R_)%*%R_%*%R_12%*%logy
  
  vmat <- R_%*%(R_12%*%logy-z%*%gammaa)
  sigma2 <- t(vmat)%*%vmat
  
  
  value<- -n/2*log(sigma2) - logR + logR_12
  #    sum(log(abs(det(diag(rep(1,n))-rho*M))))
  #   sum(log(abs(det(diag(rep(1,n))-rho*diag(c(y/h))%*%M%*%diag(c(y))))))
  return(-value)
}

ml_splogGARCHX_sato_exo(c(0.5, 0),y,z)


fun1<-function(x){
  return(max(x,-10))
}


fun2<-function(x){
  return(min(max(x,0),0.999))
  
}


 

rho<-0.03
lambda <- 0.6
gamma0<-0.4
ml_splogGARCHX_sato_exo(c(rho,lambda ),y,z)
estlim<-c(rho,lambda )
lowr<-estlim-0.02 -c(0, 0.1)
uppr<-estlim+0.1 + c(0, 0.1)



estml<- optim(estlim,fn=ml_splogGARCHX_sato_exo,y=y,z=z,method = "L-BFGS-B",hessian = T,
              lower =lowr,
              upper =uppr)

estml$par
rhoest<-estml$par[1]
lambdaest<-estml$par[2]

logy=log(y^2) 

R <- diag(n) - lambdaest * M
R_ <- solve(R)
R_12 <- diag(n) - (lambdaest + rhoest) * M 

 
gammaa_est <- solve(t(z)%*%t(R_)%*%R_%*%z)%*%t(z)%*%t(R_)%*%R_%*%R_12%*%logy


z_ <- z[,-1]

 

cest<-R_%*%(R_12%*%logy- z_%*%gammaa_est[-1,])


gamma0est<-(1-lambdaest) * log(sum(exp(cest))/n)

gamma0est

ml_splogGARCHX_exo<-function(par,y,z){
  rho<-par[1]
  lambda<-par[2]
  gamma0<-par[3]
  gamma1<-par[4]
  gamma2<-par[5]
  gamma3<-par[6]
  gamma4<-par[7]
  gamma5<-par[8]
  # gamma6<-par[9]
  # gamma7<-par[10]
  # gamma8<-par[11]
  # gamma9<-par[12]
  # gamma10<-par[13]
  # 
  
  n<-length(y)
  logM<-sum(log(rep(1,n)-rho*Meig))
  
  R<- diag(n) - lambda * M
  R_ <- solve(R)
  
  z<-z
  gammaa<- c(gamma0,gamma1,gamma2,gamma3,gamma4,gamma5)
  # ,gamma7,gamma8,gamma9,gamma10
  
  h<-exp(R_%*%(rho * M %*% log(y^2)+as.matrix(z)%*%cbind(gammaa)))
  
  v<-y/sqrt(h)
  
  logR <- sum(log(rep(1,n)-lambda*Meig))
  logR_12 <- sum(log(rep(1,n)-(lambda + rho)*Meig))
  
  value<-sum(log(dnorm(v)))-sum(log(sqrt(h))) - logR + logR_12
  #    sum(log(abs(det(diag(rep(1,n))-rho*M))))
  #   sum(log(abs(det(diag(rep(1,n))-rho*diag(c(y/h))%*%M%*%diag(c(y))))))
  return(-value)
}

gammaa<-c(gamma0est, gammaa_est[-1,])
 

ml_value <- ml_splogGARCHX_exo(c(estml$par, c(gamma0est), gammaa_est[-1,]), y, z)

2*ml_value+log(dim(ydata)[1]) * (dim(z)[2] + 2)

dim(z)[2] + 2
