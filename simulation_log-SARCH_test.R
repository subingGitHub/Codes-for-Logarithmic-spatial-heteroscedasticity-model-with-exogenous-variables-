library(MASS)

rho<-0

gamma0<-1

gamma1<-3

gamma2<-3
 
n<-dim(M)[1]

a<-c(rho,gamma0,gamma1,gamma2)


powvalueWald1<-powvalueWald2<-c()
powvalueLM1<-powvalueLM2<-c()
powvalueLR1<-powvalueLR2<-c()

#Meig<-as.matrix(eigen(M)$values)
A<-solve(diag(n)-rho*M)
W<-M


W2<-M2<-M%*%M  

M3<-M2%*%M 
M4<-M3%*%M 
#########

p1<-M 
p2<-M2-sum(diag(M2))/n*diag(n) 
p3<-M3-sum(diag(M3))/n*diag(n) 
p4<-M4-sum(diag(M4))/n*diag(n) 

# p1<-M
# p2<-M+diag(rep(c(-1,1),n/2))
# p3<-M+diag(rep(c(1,-1),n/2)) 
# p4<-M+diag(c(rep(c(1),n/2),rep(c(-1),n/2)))
#p2<-M2<-M%*%M-diag(n)*(sum(diag(M%*%M))/n)
#p3<-M3<-M%*%M%*%M-diag(n)*(sum(diag(M%*%M%*%M))/n)
#p4<-M4<-M%*%M%*%M%*%M-diag(n)*(sum(diag(M%*%M%*%M%*%M))/n)
#p2<-W+diag(rep(c(-1,1),n/2)) 
#p3<-W+diag(rep(c(1,-1),n/2)) 
#p4<-M+diag(c(rep(c(1),n/2),rep(c(-1),n/2)))

rlo<-1

while(length(powvalueWald1)<1000){
  
  eps<-rnorm(n)
  #rbin<-rbinom(n,1,0.5)
  #eps<-(rbin*rnorm(n,-2,1)+(1-rbin)*rnorm(n,2,1))/sqrt(5)
  #eps<-runif(n,-sqrt(3),sqrt(3))
  
  x0<-rep(1,n)
  x1<-rnorm(n,sd=1)
  z<-cbind(x0,x1,M%*%x1)
  x<-z
  
  gamma<-rbind(gamma0,gamma1,gamma2)
  
  y<-diag(rbinom(n,1,0.5)*2-1)%*%sqrt(exp(A%*%(as.matrix(x)%*%gamma+matrix(log(eps^2)))))
  
  
  z<-x
  #########################################
  #########################################
  #Q<-cbind(z,W%*%z,W%*%W%*%z)
  Q<-cbind(x0,x1,W%*%x1,W2%*%x1)
  
  gmm_splogARCHX1<-function(par,y,z){
    rho<-par[1]
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    n<-length(y)
    
    h<-exp(rho*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0,gamma1,gamma2))
    
    v<-y/sqrt(h)
    
    v<-v^2-1
    
    g<-t(cbind(p1%*%v,p2%*%v,p3%*%v,p4%*%v,Q))%*%v
    
    value<-t(g)%*%g
    
    #    sum(log(abs(det(diag(rep(1,n))-rho*M))))
    #   sum(log(abs(det(diag(rep(1,n))-rho*diag(c(y/h))%*%M%*%diag(c(y))))))
    return(value)
  }
  
  estml1<-try(suppressWarnings(optim(c(rho,gamma0,gamma1,gamma2),fn=gmm_splogARCHX1,y=y,z=z,method = "L-BFGS-B",hessian = T,
                                     lower =c(-0.999,-Inf,-Inf,-Inf),
                                     upper =c(0.999,Inf,Inf,Inf)))
              ,silent = TRUE)
  if('try-error' %in% class(estml1)){next}
  
  
  
  ########################################################
  covR<-function(par,y,z){
    
    rho<-par[1]
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    h<-exp(rho*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0,gamma1,gamma2))
    
    v<-y/sqrt(h)
    
    v<-v^2-1
    
    
    mu4<-c(mean(v^4))
    mu3<-c(mean(v^3))
    sigma2<-c(var(v))
    
    
    ooomga<-cbind(diag(p1),diag(p2),diag(p3),diag(p4))
    dimQ<-dim(Q)[2]
    
    omatrx1<-cbind((mu4-3*sigma2^2)*t(ooomga)%*%ooomga,mu3*t(ooomga)%*%Q)
    
    omatrx2<-cbind(mu3*t(Q)%*%ooomga,matrix(rep(0,dimQ*dimQ),nrow =dimQ))
    
    omatrx11<-rbind(omatrx1,omatrx2)
    
    omatrx3<-rbind(c(sum(diag(p1%*%(p1+t(p1)))),sum(diag(p1%*%(p2+t(p2)))),sum(diag(p1%*%(p3+t(p3)))),sum(diag(p1%*%(p4+t(p4)))),rep(0,dimQ)),
                   c(sum(diag(p2%*%(p1+t(p1)))),sum(diag(p2%*%(p2+t(p2)))),sum(diag(p2%*%(p3+t(p3)))),sum(diag(p2%*%(p4+t(p4)))),rep(0,dimQ)),
                   c(sum(diag(p3%*%(p1+t(p1)))),sum(diag(p3%*%(p2+t(p2)))),sum(diag(p3%*%(p3+t(p3)))),sum(diag(p3%*%(p4+t(p4)))),rep(0,dimQ)),
                   c(sum(diag(p4%*%(p1+t(p1)))),sum(diag(p4%*%(p2+t(p2)))),sum(diag(p4%*%(p3+t(p3)))),sum(diag(p4%*%(p4+t(p4)))),rep(0,dimQ)))
    
    omatrx<-omatrx11+sigma2^2*rbind(omatrx3,cbind(rep(0,dimQ),rep(0,dimQ),rep(0,dimQ),rep(0,dimQ),t(Q)%*%Q/sigma2))
    
    return(omatrx)
  }
  
  
  ########################################
  omatrx<-covR(estml1$par,y,z)
  
  
  gmm_splogARCHX2<-function(par,y,z){
    rho<-par[1]
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    n<-length(y)
    #logM<-sum(log(rep(1,n)-rho*Meig))
    
    #v<-(diag(n)-rho*M)%*%log(y^2)-as.matrix(z)%*%rbind(gamma0,gamma1,gamma2)+1.26
    
    
    h<-exp(rho*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0,gamma1,gamma2))
    
    v<-y/sqrt(h)
    
    v<-v^2-1
    
    g<-t(cbind(p1%*%v,p2%*%v,p3%*%v,p4%*%v,Q))%*%v     
    
    
    value<-t(g)%*%ginv(omatrx)%*%g
    
    #    sum(log(abs(det(diag(rep(1,n))-rho*M))))
    #   sum(log(abs(det(diag(rep(1,n))-rho*diag(c(y/h))%*%M%*%diag(c(y))))))
    return(value)
  }
  
  estml2<-try(suppressWarnings(optim(c(rho,gamma0,gamma1,gamma2),fn=gmm_splogARCHX2,y=y,z=z,method = "L-BFGS-B",hessian = T,
                                     lower =c(-0.999,-Inf,-Inf,-Inf),
                                     upper =c(0.999,Inf,Inf,Inf)))
              ,silent = TRUE)
  if('try-error' %in% class(estml2)){next}
  
  ############################################################
  ############################################################
  
  gmm_splogARCHX10<-function(par,y,z){
    rho<-0
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    n<-length(y)
    
    h<-exp(rho*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0,gamma1,gamma2))
    
    v<-y/sqrt(h)
    
    v<-v^2-1
    
    g<-t(cbind(p1%*%v,p2%*%v,p3%*%v,p4%*%v,Q))%*%v
    
    value<-t(g)%*%g
    
    #    sum(log(abs(det(diag(rep(1,n))-rho*M))))
    #   sum(log(abs(det(diag(rep(1,n))-rho*diag(c(y/h))%*%M%*%diag(c(y))))))
    return(value)
  }
  
  estml10<-try(suppressWarnings(optim(c(0,gamma0,gamma1,gamma2),fn=gmm_splogARCHX10,y=y,z=z,method = "L-BFGS-B",hessian = T,
                                      lower =c(-0.999,-Inf,-Inf,-Inf),
                                      upper =c(0.999,Inf,Inf,Inf)))
               ,silent = TRUE)
  if('try-error' %in% class(estml10)){next}
  
  
  
  #########################################################
  
  covR0<-function(par,y,z){
    
    rho<-0
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    h<-exp(rho*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0,gamma1,gamma2))
    
    v<-y/sqrt(h)
    
    v<-v^2-1
    
    
    mu4<-c(mean(v^4))
    mu3<-c(mean(v^3))
    sigma2<-c(var(v))
    
    ooomga<-cbind(diag(p1),diag(p2),diag(p3),diag(p4))
    dimQ<-dim(Q)[2]
    
    omatrx1<-cbind((mu4-3*sigma2^2)*t(ooomga)%*%ooomga,mu3*t(ooomga)%*%Q)
    
    omatrx2<-cbind(mu3*t(Q)%*%ooomga,matrix(rep(0,dimQ*dimQ),nrow =dimQ))
    
    omatrx11<-rbind(omatrx1,omatrx2)
    
    omatrx3<-rbind(c(sum(diag(p1%*%(p1+t(p1)))),sum(diag(p1%*%(p2+t(p2)))),sum(diag(p1%*%(p3+t(p3)))),sum(diag(p1%*%(p4+t(p4)))),rep(0,dimQ)),
                   c(sum(diag(p2%*%(p1+t(p1)))),sum(diag(p2%*%(p2+t(p2)))),sum(diag(p2%*%(p3+t(p3)))),sum(diag(p2%*%(p4+t(p4)))),rep(0,dimQ)),
                   c(sum(diag(p3%*%(p1+t(p1)))),sum(diag(p3%*%(p2+t(p2)))),sum(diag(p3%*%(p3+t(p3)))),sum(diag(p3%*%(p4+t(p4)))),rep(0,dimQ)),
                   c(sum(diag(p4%*%(p1+t(p1)))),sum(diag(p4%*%(p2+t(p2)))),sum(diag(p4%*%(p3+t(p3)))),sum(diag(p4%*%(p4+t(p4)))),rep(0,dimQ)))
    
    omatrx<-omatrx11+sigma2^2*rbind(omatrx3,cbind(rep(0,dimQ),rep(0,dimQ),rep(0,dimQ),rep(0,dimQ),t(Q)%*%Q/sigma2))
    
    return(omatrx)
  }
  
  
  ########################################
  omatrx0<-covR0(estml10$par,y,z)
  
  
  gmm_splogARCHX20<-function(par,y,z){
    rho<-0
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    n<-length(y)
    #logM<-sum(log(rep(1,n)-rho*Meig))
    
    #v<-(diag(n)-rho*M)%*%log(y^2)-as.matrix(z)%*%rbind(gamma0,gamma1,gamma2)+1.26
    
    
    h<-exp(rho*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0,gamma1,gamma2))
    
    v<-y/sqrt(h)
    
    v<-v^2-1
    
    g<-t(cbind(p1%*%v,p2%*%v,p3%*%v,p4%*%v,Q))%*%v     
    
    
    value<-t(g)%*%ginv(omatrx0)%*%g
    
    #    sum(log(abs(det(diag(rep(1,n))-rho*M))))
    #   sum(log(abs(det(diag(rep(1,n))-rho*diag(c(y/h))%*%M%*%diag(c(y))))))
    return(value)
  }
  
  estml20<-try(suppressWarnings(optim(c(0,gamma0,gamma1,gamma2),fn=gmm_splogARCHX20,y=y,z=z,method = "L-BFGS-B",hessian = T,
                                      lower =c(-0.999,-Inf,-Inf,-Inf),
                                      upper =c(0.999,Inf,Inf,Inf)))
               ,silent = TRUE)
  if('try-error' %in% class(estml20)){next}
  
  
  
  powvalueLR1<-c(powvalueLR1,(estml20$value-estml2$value))
  
  ################################################
  ################################################
  
  Dmatrix0<-function(par,y,z){
    rho<-0
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    h<-exp(rho*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0,gamma1,gamma2))
    
    v<-y/sqrt(h)
    
    Evlogv<-c(mean((v^2-1)*log(v^2)))
    Ev2<-c(mean((v^2)))
    Elogy<-c(mean(log(v^2)))
    Eylogy<-c(mean(v^2*log(v^2)))
    
    B<-M%*%solve(diag(n)-rho*M)
    vecB<-matrix(c(diag(B)),nrow=n)
    
    matrix1<-matrix(apply(B,1,sum))*Elogy*Ev2+vecB*(Eylogy-Elogy*Ev2)+B%*%z%*%gamma
    
    
    D<-t(rbind(cbind(Evlogv*Ev2*sum(diag((p1+t(p1))%*%B)),Evlogv*Ev2*sum(diag((p2+t(p2))%*%B)),Evlogv*Ev2*sum(diag((p3+t(p3))%*%B)),Evlogv*Ev2*sum(diag((p4+t(p4))%*%B)), t(matrix1)%*%Q),
               cbind(matrix(rep(0,3*4),nrow=3),t(z)%*%Q)))
    
    return(D)
    
  }
  
  Rfun0<-function(par,y,z){
    rho<-0
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    n<-length(y)
    #logM<-sum(log(rep(1,n)-rho*Meig))
    
    #v<-(diag(n)-rho*M)%*%log(y^2)-as.matrix(z)%*%rbind(gamma0,gamma1,gamma2)+1.26
    
    
    h<-exp(rho*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0,gamma1,gamma2))
    
    v<-y/sqrt(h)
    
    v<-v^2-1
    
    g<-t(cbind(p1%*%v,p2%*%v,p3%*%v,p4%*%v,Q))%*%v     
    
    
    
    return(g)
  }
  
  omatrx0<-covR0(estml20$par,y,z)
  D0<-Dmatrix0(estml20$par,y,z)
  R0<-Rfun0(estml20$par,y,z)
  
  value10<-t(D0)%*%ginv(omatrx0)%*%R0
  value20<-t(D0)%*%ginv(omatrx0)%*%D0
  
  powvalueLM1<-c(powvalueLM1,  t(value10)%*%ginv(value20)%*%value10)
  
  
  ##########################################
  ########################################
  #####################################################
  #####################################################
  
  Dmatrix<-function(par,y,z){
    rho<-par[1]
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    h<-exp(rho*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0,gamma1,gamma2))
    
    v<-y/sqrt(h)
    
    Evlogv<-c(mean((v^2-1)*log(v^2)))
    Ev2<-c(mean((v^2)))
    Elogy<-c(mean(log(v^2)))
    Eylogy<-c(mean(v^2*log(v^2)))
    
    B<-M%*%solve(diag(n)-rho*M)
    vecB<-matrix(c(diag(B)),nrow=n)
    
    matrix1<-matrix(apply(B,1,sum))*Elogy*Ev2+vecB*(Eylogy-Elogy*Ev2)+B%*%z%*%gamma
    
    D<-t(rbind(cbind(Evlogv*Ev2*sum(diag((p1+t(p1))%*%B)),Evlogv*Ev2*sum(diag((p2+t(p2))%*%B)),Evlogv*Ev2*sum(diag((p3+t(p3))%*%B)),Evlogv*Ev2*sum(diag((p4+t(p4))%*%B)),t(matrix1)%*%Q),
               cbind(matrix(rep(0,3*4),nrow=3),t(z)%*%Q)))
    
    return(D)
  }
  
  
  omatrx<-covR(estml2$par,y,z)
  D<-Dmatrix(estml2$par,y,z)
  
  est<-matrix(estml2$par)
  est0<-matrix(c(0,gamma0,gamma1,gamma2))
  
  rmatrix<-matrix(c(1,0,0,0),ncol=1)
  powvalueWald1<-c(powvalueWald1, t(est[1])%*%ginv(t(rmatrix)%*%ginv(t(D)%*%ginv(omatrx)%*%D)%*%rmatrix)%*%(est[1]))
  
  #################################
  ##################
  
  
  if(length(powvalueWald1)%%50==0){
    print(rlo)
    rlo<-rlo+1
  }
  
  
}



c(round(sum(abs(powvalueWald1)>qchisq(0.99,1))/length(powvalueWald1),3),
  round(sum(abs(powvalueLM1)>qchisq(0.99,1))/length(powvalueLM1),3),
  round(sum(abs(powvalueLR1)>qchisq(0.99,1))/length(powvalueLR1),3))

c(round(sum(abs(powvalueWald1)>qchisq(0.95,1))/length(powvalueWald1),3), 
  round(sum(abs(powvalueLM1)>qchisq(0.95,1))/length(powvalueLM1),3),
  round(sum(abs(powvalueLR1)>qchisq(0.95,1))/length(powvalueLR1),3))

c(round(sum(abs(powvalueWald1)>qchisq(0.90,1))/length(powvalueWald1),3),
  round(sum(abs(powvalueLM1)>qchisq(0.90,1))/length(powvalueLM1),3),
  round(sum(abs(powvalueLR1)>qchisq(0.90,1))/length(powvalueLR1),3))


 
