load("F:/_论文/_paper/reviewed paper/logARCH/logSARCH (paper and coding)/coding/empirical/empirical data/house data/house1993_SAR.RData")
library(Rsolnp)



z_raw = z 

dim(z)

ml_splogARCHX_durbin<-function(par,y,z){
  rho<-par[1]
 
  gamma0<-par[2]
  gamma1<-par[3]
  gamma2<-par[4]
  gamma3<-par[5]
  gamma4<-par[6]
  gamma5<-par[7]
  
  gamma6<-par[8]
  gamma7<-par[9]
  gamma8<-par[10]
  gamma9<-par[11]
  gamma10<-par[12]
  
  
  n<-length(y)
  logM<-sum(log(rep(1,n)-rho*Meig))
  
  h<-exp(rho*M%*%log(y^2)+z%*%rbind(gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8,gamma9,gamma10))
  
  v<-y/sqrt(h)
  
  
  value<-sum(log(dnorm(v)))-sum(log(sqrt(h)))+logM
  #    sum(log(abs(det(diag(rep(1,n))-rho*M))))
  #   sum(log(abs(det(diag(rep(1,n))-rho*diag(c(y/h))%*%M%*%diag(c(y))))))
  return(-value)
}


fun1<-function(x){
  return(max(x,-10))
}

fun2<-function(x){
  return(min(max(x,-0.999),0.999))
  
}


# 假设你的数据是 z
z_scaled <- z  # 创建副本以保留原始数据

# 对第2到第10列进行标准化（列均值为0，标准差为1）
z_scaled[, 2:11] <- scale(z[, 2:11])


z<-z_scaled

rho<-0.3
gamma0<-0
gamma1<-0
gamma2<-0
gamma3<-0
gamma4<-0
gamma5<-0
gamma6<-0
gamma7<-0
gamma8<-0
gamma9<-0
gamma10<-0


ml_splogARCHX_durbin(c(rho,gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8,gamma9,gamma10),y,z)
estlim<-c(rho,gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8,gamma9,gamma10)
alim<-0.2
lowr<-c(rho,gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8,gamma9,gamma10)-alim
uppr<-c(rho,gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8,gamma9,gamma10)+alim



for(i in 1:1000){
  estml<-#solnp(estlim,fun=ml_splogARCHX_durbin,y=y,z=z,
    #      LB=lowr,
    #      UB=uppr)
    optim(estlim,fn=ml_splogARCHX_durbin,y=y,z=z,method = "L-BFGS-B",hessian = T,
          lower =lowr,
          upper =uppr)
  print(estml$par)
  print(i)
  print(2*estml$value+2*12)
  estlim<-estml$par
  lowr<-c(apply(t(as.matrix(estlim[1]-alim)),2,fun2),c(estlim[2:12]-alim))
  uppr<-c(apply(t(as.matrix(estlim[1]+alim)),2,fun2),c(estlim[2:12]+alim))
}


estml$par
 



##############

n<-length(y)

round(2*estml$value+log(n)*length(lowr),3)

#############



rhoest<-estml$par[1]
gammaest<-estml$par[2:12]
 dim(z)
 dim(cbind(c(estml$par[2:12])))

hest<-exp(rhoest*M%*%log(y^2)+ z%*%cbind(c(estml$par[2:12])))

vest<-y/sqrt(hest)

l_n <- matrix(1, nrow = n, ncol = 1)

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




Omgea_ML<-function(par){
  rho<-par[1]
  gamma0<-par[2]
  gamma1<-par[3]
  gamma2<-par[4]
  gamma3<-par[5]
  gamma4<-par[6]
  gamma5<-par[7]
  
  gamma6<-par[8]
  gamma7<-par[9]
  gamma8<-par[10]
  gamma9<-par[11]
  gamma10<-par[12]
  
  
  n<-length(y)
  logM<-sum(log(rep(1,n)-rho*Meig))
  
  gamma<-rbind(gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8,gamma9,gamma10)
  
  h<-exp(rho*M%*%log(y^2)+z%*%rbind(gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8,gamma9,gamma10)) 
  
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
  gamma3<-par[5]
  gamma4<-par[6]
  gamma5<-par[7]
  
  gamma6<-par[8]
  gamma7<-par[9]
  gamma8<-par[10]
  gamma9<-par[11]
  gamma10<-par[12]
  
  
  n<-length(y)
  logM<-sum(log(rep(1,n)-rho*Meig))
  
  gamma<-rbind(gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8,gamma9,gamma10)
  
  h<-exp(rho*M%*%log(y^2)+z%*%rbind(gamma0,gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8,gamma9,gamma10)) 
  
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


library(MASS)

Sigma_mat <- Sigma_ML(estml$par)
ml_ad <- sqrt(diag(ginv(Sigma_mat)%*%Omgea_mat%*%ginv(Sigma_mat)))


estml$par / ml_ad

tvalue<- estml$par / ml_ad

2 * (1 - pnorm(abs(tvalue))) * 10


0, 1, 3,4,6, 10 




save.image("C:/Users/subing/Desktop/house1998_SARCH_nonexo.RData")

