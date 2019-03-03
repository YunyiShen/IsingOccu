require('gtools')
set.seed(42)
source('Moller_repeat_spin_glass.r')
#source("misc.R")
nlat = 12

nrep = 4

siteposi = 1.00 * permutations(n=nlat,r=2,v=(1:nlat),repeats.allowed = T)

distM = as.matrix((dist(siteposi)))
distM=distM==1
diag(distM)=0
distM = distM * 1.0

ones = rep(1,times = nlat*nlat)
X = cbind(ones)
theta = matrix(c(-.0,0,.15,.15,-.2))


set.seed(45)
Zsample = rIsing(X,distM,theta,method = "CFTP",nIter=100,n=nrep)
raster::plot(raster::raster(matrix(Zsample,nlat,nlat)))

Hamiltonian(theta, X, distM, Zsample)

#var_prop = c(2.5e-5,2.5e-5,2.5e-5,2.5e-5,2.5e-5)
#var_prop = c(5e-5,5e-5,7e-5,7e-5,4e-4)

var_prop = c(rep(5e-5,2),rep(1e-4,2),5e-4)

kk=Moller.sampler_repeat(X=X,distM=distM,
                                      #detmat = detmat, 
                                      #detX=detX, 
                                      Z=Zsample,
                                      mcmc.save = 5000, burn.in = 1000 , 
                                      vars_prior = 100000,
                                      vars_prop = var_prop
                                      ,seed = 42
                                      ,init = theta
                                      , thin.by = 1)
plot(kk$theta.mcmc[,1])
plot(kk$theta.mcmc[,2])
plot(kk$theta.mcmc[,3])
plot(kk$theta.mcmc[,4])
plot(kk$theta.mcmc[,5])
sqrt(apply(kk$theta.mcmc,2,var))
