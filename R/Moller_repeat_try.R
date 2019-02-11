require('gtools')
set.seed(42)
source('Moller_repeat.R')
#source("misc.R")
nlat = 20
siteposi = 1.00 * permutations(n=nlat,r=2,v=(1:nlat),repeats.allowed = T)

distM = as.matrix((dist(siteposi)))
distM = 1*(distM==1)

ones = rep(1,times = nlat*nlat)
X = cbind(ones)
theta = matrix(c(-0,0.3))

set.seed(12345)
Zsample = rIsing(X,distM,theta,method = "CFTP",nIter=100,n=1)

var_prop = c(rep(2.5e-5,2))

kk=Moller.sampler_repeat(X=X,distM=distM,
                                      #detmat = detmat, 
                                      #detX=detX, 
                                      Z=Zsample,
                                      mcmc.save = 1000, burn.in = 100 , 
                                      vars_prior = 100000,
                                      vars_prop = var_prop
                                      ,seed = 42
                                      ,init = theta
                                      , thin.by = 1)
plot(kk$theta.mcmc[,1])
plot(kk$theta.mcmc[,2])
