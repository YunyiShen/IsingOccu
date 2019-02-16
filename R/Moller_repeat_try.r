require('gtools')
set.seed(42)
source('Moller_repeat.R')
#source("misc.R")
nlat = 20
siteposi = 1.00 * permutations(n=nlat,r=2,v=(1:nlat),repeats.allowed = T)

distM = as.matrix((dist(siteposi)))
distM=distM-1
diag(distM)=0
dist_thr = 5
#distM = 1*(distM==1)

ones = rep(1,times = nlat*nlat)
X = cbind(ones)
theta = matrix(c(-0,0.15,2))


set.seed(42)
Zsample = rIsing(X,distM,theta,method = "CFTP",nIter=100,n=1,dist_thr)
raster::plot(raster::raster(matrix(Zsample,nlat,nlat)))

var_prop = c(2.5e-5,2.5e-5,1e-4)

kk=Moller.sampler_repeat(X=X,distM=distM,
                                      #detmat = detmat, 
                                      #detX=detX, 
                                      Z=Zsample,
                                      mcmc.save = 1500, burn.in = 1000 , 
                                      vars_prior = 100000,
                                      vars_prop = var_prop
                                      ,seed = 42
                                      ,init = theta
                                      , thin.by = 1,dist_thr = 5)
plot(kk$theta.mcmc[,1])
plot(kk$theta.mcmc[,2])
plot(kk$theta.mcmc[,3])
sqrt(apply(kk$theta.mcmc,2,var))
