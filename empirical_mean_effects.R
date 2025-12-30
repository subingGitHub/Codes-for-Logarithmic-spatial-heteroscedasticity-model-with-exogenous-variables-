
library(spatialreg)

load("C:/Users/subing/Desktop/log-SHE two new empirical data/house data/house1998_SAR.RData")

summary(model.SLM.ML)

 
gamma1<-0.539
gamma2<-0.027
gamma3<-0.103
gamma4<--0.504


gamma1M<--0.257
gamma2M<-0.031
gamma3M<-0
gamma4M<-0





n<-dim(M)[1]

effect<-solve(diag(n)-c(model.SLM.ML$rho)*M)


effect1<-effect%*%(gamma1*diag(n)+gamma1M*M)
effect2<-effect%*%(gamma2*diag(n)+gamma2M*M)
effect3<-effect%*%(gamma3*diag(n)+gamma3M*M)
effect4<-effect%*%(gamma4*diag(n)+gamma4M*M)


dir_1<-round(mean(diag(effect1)),3)
dir_2<-round(mean(diag(effect2)),3)
dir_3<-round(mean(diag(effect3)),3)
dir_4<-round(mean(diag(effect4)),3)
 



undir_1<-round((sum(effect1)-sum(diag(effect1)))/n,3)
undir_2<-round((sum(effect2)-sum(diag(effect2)))/n,3)
undir_3<-round((sum(effect3)-sum(diag(effect3)))/n,3)
undir_4<-round((sum(effect4)-sum(diag(effect4)))/n,3)

total_1<-round(sum(effect1)/n,3)
total_2<-round(sum(effect2)/n,3)
total_3<-round(sum(effect3)/n,3)
total_4<-round(sum(effect4)/n,3)

cat(total_1,'&',total_2,'&',total_3,'&',total_4)
cat(dir_1,'&',dir_2,'&',dir_3,'&',dir_4)
cat(undir_1,'&',undir_2,'&',undir_3,'&',undir_4)

# cat(dir_1,'&',undir_1,'&',total_1)
# cat(dir_2,'&',undir_2,'&',total_2)
# cat(dir_3,'&',undir_3,'&',total_3)
# cat(dir_4,'&',undir_4,'&',total_4)





