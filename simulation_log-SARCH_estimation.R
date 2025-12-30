library(MASS)



rho<-0.3

gamma0<-1

gamma1<-3

gamma2<-3

# W <- W[1:50, 1:50]
# W <- W / apply(W,1, sum)
M<-as.matrix(W)
Meig <- as.matrix(eigen(M)$values)

n<-dim(M)[1] 

W2<-W%*%W  

M2<-M%*%M  
M3<-M2%*%M 
M4<-M3%*%M 
#########
# p1<-M+diag(c(rep(c(1),n/2),rep(c(-1),n/2)))
# p2<-M2-diag(diag(M2))+diag(c(rep(c(1),n/2),rep(c(-1),n/2)))
# p3<-M3-diag(diag(M3))+diag(c(rep(c(1),n/2),rep(c(-1),n/2)))
# p4<-M4-diag(diag(M4))+diag(c(rep(c(1),n/2),rep(c(-1),n/2)))

# q1<-M
# p2<-M2-diag(diag(M2))
# p3<-M3-diag(diag(M3))
# p4<-M4-diag(diag(M4))

p1<-M
p2<-M2-sum(diag(M2))/n*diag(n)
p3<-M3-sum(diag(M3))/n*diag(n)
p4<-M4-sum(diag(M4))/n*diag(n)
#p2<-p2/apply(p2,1,sum)
#p3<-p3/apply(p3,1,sum)
#p4<-p4/apply(p4,1,sum)
a<-b<-c<-d<-cc<-matrix(c(rho,gamma0,gamma1,gamma2),nrow=1)


a_ad <- b_ad <- c_ad <- matrix(c(0,0,0,0),nrow=1) 

rlo<-1

A<-(diag(n)-rho*M)

A_bar <-solve(A)

A_dot <- - M

A_dotdot<- 0

A_star = A_dot %*% A_bar


trace_A_star <- sum(diag(A_star))

  

l_n <- matrix(1, nrow = n, ncol = 1)

while(dim(a)[1]<1000){
  
  #eps<-rnorm(n)                                  
  rbin<-rbinom(n,1,0.5)
  eps<-(rbin*rnorm(n,-2,1)+(1-rbin)*rnorm(n,2,1))/sqrt(5)
  #eps<-runif(n,-sqrt(3),sqrt(3))
  #rbin<-rbinom(n,1,0.5)
  #eps<-(rbin*rnorm(n,-1,1)+(1-rbin)*rnorm(n,1,1))/sqrt(2)
  #eps<-(runif(n,-2,2))*sqrt(3/4)
  
  
  x0<-rep(1,n)
  x1<-rnorm(n,sd=1)
  z<-cbind(x0,x1,M%*%x1)
  x<-z
  
  A<-(diag(n)-rho*M)
  
  A_bar <-solve(A)
  
  gamma<-rbind(gamma0,gamma1,gamma2)
  
  y<-diag(sign(eps))%*%sqrt(exp(A_bar%*%as.matrix(z)%*%gamma+A_bar%*%matrix(log(eps^2))))
  
  ml_splogARCHX<-function(par,y,z){
    rho<-par[1]
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    n<-length(y)
    logM<-sum(log(rep(1,n)-rho*Meig))
    
    h<-exp(rho*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0,gamma1,gamma2))
    
    v<-y/sqrt(h)
    
    
    value<-sum(log(dnorm(v)))-sum(log(sqrt(h)))+logM
    #    sum(log(abs(det(diag(rep(1,n))-rho*M))))
    #   sum(log(abs(det(diag(rep(1,n))-rho*diag(c(y/h))%*%M%*%diag(c(y))))))
    return(-value)
  }
  
  estml<-try(suppressWarnings(optim(c(rho,gamma0,gamma1,gamma2),fn=ml_splogARCHX,y=y,z=z,method = "L-BFGS-B",hessian = T,
                                    lower =c(-0.99,-Inf,-Inf,-Inf),
                                    upper =c(0.99,Inf,Inf,Inf)))
             ,silent = TRUE)
  if('try-error' %in% class(estml)){next}
  
  
  
  ######################################################
  
  rhoest<-estml$par[1]
  gamma0est<-estml$par[2]
  gamma1est<-estml$par[3]
  gamma2est<-estml$par[4]
  
  hest<-exp(rhoest*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0est,gamma1est,gamma2est))
  
  vest<-y/sqrt(hest)
  
  
  mu4<-c(mean(vest^4))
  mu3<-c(mean(vest^3))
  sigma2<-c(var(vest^2))
  
  a_e <- c(mean((vest^2 - 1)^2 * log(vest^2)))
  
  b_e <- c(mean(log(vest^2)))
  
  c_e <- c(mean((vest^2 - 1)^2 * (log(vest^2))^2 )) 
  
  d_e <- c(mean((vest^2 - 1) * log(vest^2) ))
  
  e_e <- c(mean((log(vest^2))^2 ))
  
  f_e <- a_e - sigma2 * b_e
  
  
  sigma2_star <- c(mean(vest^2))
  
  a_e_star <- c(mean(vest^2 * log(vest^2)))
  
  c_e_star <- c(mean(vest^2 * (log(vest^2))^2 )) 
  
  d_e_star <- c(mean(vest^2 * log(vest^2) ))
  
  f_e_star <- a_e_star - sigma2_star * b_e 
  
  
  
  
  
  
  # E_rhorho_star <- 1/ 2 * sigma2_star * b_e^2 * c(t(A_star%*%l_n)%*%(A_star%*%l_n)) - sigma2_star * b_e * c( t(A_star%*%l_n)%*%(A_star%*%z%*%gamma) ) +
  #   1/2 * sigma2_star * c(t(A_star%*%z%*%gamma)%*%(A_star%*%z%*%gamma)) + 1/2 * sigma2_star * (e_e - b_e^2 ) * sum(diag(A_star%*%t(A_star)))  - sum(diag(A_star%*%A_star))
  # 
  # 
  # 
  # E_rho_gamma_star <- sigma2_star * t(b_e * A_star %*%l_n  + A_star %*%z%*%gamma)%*%z
  #   
  # 
  # E_gammagamma_star <- 1/4* sigma2_star * t(z)%*%z
  # 
  #   E_rhorho<- 1/2 * sigma2 * b_e * c(t(l_n) %*%t(A_star) %*% A_star %*% z %*%gamma) +
  #     1/4 * t(A_star%*%z%*%gamma) %*% (A_star%*%z%*%gamma) + 
  #     1/4 * ( (d_e^2 - 4) * trace_A_star^2+ sigma2 * (e_e - b_e^2) * sum(diag(A_star%*%t(A_star))) + d_e^2 *  sum(diag(A_star%*%A_star)) )
  #   
  #   
  #   E_rho_gamma <- sigma2 * t(b_e * A_star %*%l_n  + A_star %*%z%*%gamma)%*%z
  #   
  #   
  #   E_gammagamma<- 1/4* sigma2 * t(z)%*%z
  
  
  
  
  
  Omgea_ML<-function(par){
    rho<-par[1]
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    gamma<-rbind(gamma0,gamma1,gamma2)
    
    A<-(diag(n)-rho*M)
    
    A_bar <-solve(A)
    
    A_dot <- - M
    
    A_dotdot<- 0 * M
    
    A_star <- A_dot %*% A_bar 
    
    vec_A_star <- as.matrix(diag(A_star))
    
    trace_A_star <- sum(diag(A_star))
    
    E_rhorho<-  1/4 *(c_e - 2*a_e * b_e - 2*d_e^2 - sigma2 * e_e + 3* sigma2 * b_e^2) * t(vec_A_star)%*%vec_A_star + 
      1/2 *(a_e - sigma2* b_e) * b_e * t(vec_A_star)%*%A_star%*%l_n + 
      1/2 * (a_e - sigma2 *b_e) * t(vec_A_star)%*%A_star%*%z%*%gamma +
      1/4 * sigma2 * b_e^2 * t(A_star%*%l_n)%*%(A_star%*%l_n) + 
      1/2 * sigma2 * b_e * c(t(l_n) %*%t(A_star) %*% A_star %*% z %*%gamma) +
      1/4 * sigma2 * t(A_star%*%z%*%gamma) %*% (A_star%*%z%*%gamma) + 
      1/4 * ( (d_e^2 + 4 * (1- d_e)) * trace_A_star^2 + sigma2 * (e_e - b_e^2) * sum(diag(A_star%*%t(A_star))) ) + 
      1/4 * d_e^2 *  sum(diag(A_star%*%A_star)) 
    
    
    
    E_rho_gamma <- 1/4 * f_e * t(vec_A_star)%*%z + 1/4 * sigma2 * t(b_e * A_star %*%l_n  + A_star %*%z%*%gamma)%*%z
    
    
    E_gammagamma<- 1/4* sigma2 * t(z)%*%z
    
    
    
    return(rbind(cbind(E_rhorho, E_rho_gamma ),
                 cbind(t(E_rho_gamma), E_gammagamma)))
    
  }
  
  Omgea_mat <- Omgea_ML(estml$par)
  
  
  
  
  
  
  Sigma_ML <-function(par){
    
    rho<-par[1]
    gamma0<-par[2]
    gamma1<-par[3]
    gamma2<-par[4]
    
    gamma<-rbind(gamma0,gamma1,gamma2)
    
    A<-(diag(n)-rho*M)
    
    A_bar <-solve(A)
    
    A_dot <- - M
    
    A_dotdot<- 0 * M
    
    A_star <- A_dot %*% A_bar 
    
    vec_A_star <- as.matrix(diag(A_star))
    
    E_rhorho_star <-   -1/2 * (c_e_star  - 2* a_e_star * b_e - sigma2_star * e_e + 2*sigma2_star *b_e^2) * t(vec_A_star)%*%vec_A_star -
      (a_e_star * b_e - sigma2_star * b_e^2)* t(vec_A_star)%*%A_star%*%l_n - 
      (a_e_star - sigma2 * b_e) * t(vec_A_star)%*%A_star%*%z%*%gamma  - 
      1/ 2 * sigma2_star * b_e^2 * c(t(A_star%*%l_n)%*%(A_star%*%l_n)) -
      sigma2_star * b_e * c( t(A_star%*%l_n)%*%(A_star%*%z%*%gamma) ) -
      1/2 * sigma2_star * c(t(A_star%*%z%*%gamma)%*%(A_star%*%z%*%gamma)) -
      1/2 * sigma2_star * (e_e - b_e^2 ) * sum(diag(A_star%*%t(A_star))) -
      sum(diag(A_star%*%A_star))  + (1 - 1/2 *a_e) * sum(diag(A_dotdot%*%A_bar))
    
    
    
    E_rho_gamma_star <- -1/2 * f_e_star * t(vec_A_star)%*%z - 1/2 * sigma2_star * t(b_e * A_star %*%l_n  + A_star %*%z%*%gamma)%*%z  
    
    
    E_gammagamma_star <- - 1/2* sigma2_star * t(z)%*%z
    
    
    
    return(-rbind(cbind(E_rhorho_star, E_rho_gamma_star ),
                  cbind(t(E_rho_gamma_star), E_gammagamma_star)))
    
  }
  
  
  
  
  Sigma_mat <- Sigma_ML(estml$par)
  ml_ad <- sqrt(diag(ginv(Sigma_mat)%*%Omgea_mat%*%ginv(Sigma_mat)))
   
  #######################################################################  
  #sec<-secparfun(c(estml1$par),y,x)
  # bias<-solve(-sec)%*%matrix(c(biasfun(estml1$par,y,x),0,0,0),ncol=1)
  
  
  # b<-rbind(b,estml1$par-t(bias))
  ###########################################################################
  
  
  ##########################################################
  
  #Q<-cbind(z,W%*%z,W2%*%z)
  Q<-cbind(x0,x1,W%*%x1,W2%*%x1)
  
  gmm_splogARCHX1<-function(par,y,z){
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
  
  
  
  #########################################################
 
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
  
  D<-Dmatrix(estml1$par,y,z)
  omatrx<-covR(estml1$par,y,z)
  
  Sigma_gmm <- ginv(t(D)%*%D)%*%(t(D)%*%omatrx%*%D)%*%ginv(t(D)%*%D)
  
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
  
  
  
  D<-Dmatrix(estml2$par,y,z)
  omatrx<-covR(estml2$par,y,z)
  
  Sigma_ogmm <- ginv(t(D)%*%ginv(omatrx)%*%D)
  
  
  ##########################################################
  
  
  
  
  #  B<-M%*%solve(diag(n)-rhoest*M)
  # vecB<-matrix(c(diag(B)),nrow=n)
  
  # mulogy<-c(mean(log(vest^2)))
  # muylogy<-c(mean(vest^2*log(vest^2)))
  
  # 
  # p5<-B-diag(B)
  
  # Q<-cbind(matrix(apply(B,1,sum))*mulogy+vecB*(muylogy-mulogy)+B%*%z%*%gammaest,z)
  
  # dimQ<-dim(Q)[2]
  # ooomga<-cbind(diag(p5))
  
  # omatrx1<-cbind((mu4-3*sigma2^2)*t(ooomga)%*%ooomga,mu3*t(ooomga)%*%Q)
  
  # omatrx2<-cbind(mu3*t(Q)%*%ooomga,matrix(rep(0,dimQ*dimQ),nrow =dimQ))
  
  # omatrx11<-rbind(omatrx1,omatrx2)
  
  
  # omatrx3<-rbind(c(sum(diag(p5%*%(p5+t(p5)))),rep(0,dimQ)))
  
  # omatrx<-omatrx11+sigma2^2*rbind(omatrx3,cbind(rep(0,dimQ),t(Q)%*%Q/sigma2))
  
  
  
  
  #  gmm_splogARCHX4<-function(par,y,z){
  #   rho<-par[1]
  #    gamma0<-par[2]
  #   gamma1<-par[3]
  #   gamma2<-par[4]
  
  #   n<-length(y)
  #logM<-sum(log(rep(1,n)-rho*Meig))
  
  #v<-(diag(n)-rho*M)%*%log(y^2)-as.matrix(z)%*%rbind(gamma0,gamma1,gamma2)+1.26
  
  
  #   h<-exp(rho*M%*%log(y^2)+as.matrix(z)%*%rbind(gamma0,gamma1,gamma2))
  
  #   v<-y/sqrt(h)
  
  #   v<-v^2-1
  
  #   g<-t(cbind(p5%*%v,Q))%*%v     
  
  
  #   value<-t(g)%*%ginv(omatrx)%*%g
  
  #    sum(log(abs(det(diag(rep(1,n))-rho*M))))
  #   sum(log(abs(det(diag(rep(1,n))-rho*diag(c(y/h))%*%M%*%diag(c(y))))))
  #  return(value)
  # }
  
  # estml4<-try(suppressWarnings(optim(c(rho,gamma0,gamma1,gamma2),fn=gmm_splogARCHX4,y=y,z=z,method = "L-BFGS-B",hessian = T,
  #                                    lower =c(-0.999,-Inf,-Inf,-Inf),
  #                                   upper =c(0.999,Inf,Inf,Inf)))
  #             ,silent = TRUE)
  #if('try-error' %in% class(estml4)){next}
  
  
  #########################
  
  ml_splogARCHX_sato<-function(par,y,z){
    rho<-par[1]
    n<-length(y)
    logM<-sum(log(rep(1,n)-rho*Meig))
    
    logy=log(y^2)
    
    
    gammaa<-solve(t(z)%*%z)%*%t(z)%*%(diag(n)-rho*M)%*%logy
    vmat<-(diag(n)-rho*M)%*%logy-z%*%gammaa
    sigma2<-t(vmat)%*%vmat
    
    
    value<- -n/2*log(sigma2)+logM
    #    sum(log(abs(det(diag(rep(1,n))-rho*M))))
    #   sum(log(abs(det(diag(rep(1,n))-rho*diag(c(y/h))%*%M%*%diag(c(y))))))
    return(-value)
  }
  
  estml_sato<-optim(rho,fn=ml_splogARCHX_sato,y=y,z=z,method = "L-BFGS-B",hessian = T,
                    lower = -0.999,
                    upper = 0.999)
  
  rhoest<-estml_sato$par
  logy=log(y^2)
  gammaest<-solve(t(z)%*%z)%*%t(z)%*%(diag(n)-rhoest*M)%*%logy
  
  cest<-(diag(n)-rhoest*M)%*%logy-z%*%gammaest+gammaest[1]
  
  gamma0est<-log(sum(exp(cest))/n)
  
  
  
  a<-rbind(a,estml$par)
  a_ad <-rbind(a_ad, ml_ad)
  
  
  b<-rbind(b,estml1$par)
  b_ad <-rbind(b_ad, sqrt(diag(Sigma_gmm)))
  c<-rbind(c,estml2$par)
  c_ad <-rbind(c_ad, sqrt(diag(Sigma_ogmm)))
  
  
  d<-rbind(d,c(rhoest,gamma0est,c(gammaest)[2:3]))
  
  if(dim(a)[1]%%50==0){
    print(rlo)
    rlo<-rlo+1
  }
  
}





a[1,]

apply(a,2,mean)-a[1,]
apply(b,2,mean)-a[1,]
apply(c,2,mean)-a[1,]
apply(d,2,mean)-a[1,]




sqrt(apply(a,2,var))
sqrt(apply(b,2,var))
sqrt(apply(c,2,var))
sqrt(apply(d,2,var))



apply(a_ad,2,mean)
apply(b_ad,2,mean)
apply(c_ad,2,mean)                 

