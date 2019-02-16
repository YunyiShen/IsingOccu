require('gtools')
set.seed(42)
source('Sampler_withZ.R')
source("misc.R")
nlat = 20
siteposi = 1.00 * permutations(n=nlat,r=2,v=(1:nlat),repeats.allowed = T)

distanceM = as.matrix((dist(siteposi)))
#distanceM = distanceM-1
distanceM= 1*( distanceM==1)
diag(distanceM) = 0

ones = rep(1,times = nlat*nlat)
X = cbind(ones)

detX = list()
nperiod = 10
for (i in 1:nperiod){
  temp = matrix(runif(nlat^2),nrow = nlat^2,ncol = 1)
	detX[[i]] = temp
}

theta = matrix(c(-0, # env reaction of 1
                 -0,  # env reaction of 2
                 0,1,  # detection beta of 1
                 0,1,   # detection beta of 2
                 .15,3,        # eta01 d1
                 0.15,3,		  # eta02 d2
                 -0.5))

p = length(theta)
sites = nrow(distanceM)
ncov = ncol(X)
zeros = matrix(0,nrow=nlat^2,ncol=ncov)
beta1 = as.numeric( matrix(c(theta[1:(2*ncol(X))])))
Xfull = cbind(rbind(X,zeros),rbind(zeros,X))
thr = Xfull%*%beta1

set.seed(42)
Zsample = rIsingOccu(X,distanceM,theta,method = "CFTP",nIter=100,n=1,int_range = "nn")

raster::plot(raster::raster(
  matrix(
    Zsample
    [1:nlat^2],
    nrow=nlat,ncol=nlat)))

raster::plot(raster::raster(
  matrix(
    Zsample
    [1:nlat^2 + nlat^2],
    nrow=nlat,ncol=nlat)))

raster::plot(raster::raster(matrix(
  Zsample[1:nlat^2 + nlat^2] + 
  Zsample[1:nlat^2] ,
  nrow=nlat,ncol=nlat
  
)))
#raster::plot(raster::raster(matrix(detmat[401:800,1],nrow = 20,ncol=20)))
detmat = matrix(nrow = length(Zsample),ncol = nperiod)
# Test the detection function
detSample = IsingOccu_sample.detection(theta, X,  Z=Zsample, detmat, detX)
detmat = detSample

distM = distanceM
MPLE = optim(((theta)),IsingOccu.logPL,envX=X, distM=distanceM, Z=Zsample ,detmat=detmat, detX=detX, int_range = "nn",control  = list(maxit=5000))
IsingOccu.logPL(theta, X, distM, Zsample ,detmat, detX, int_range = "nn")



var_prop = c(rep(2.5e-5,2),rep(2.5e-3,4),rep(2.5e-5,5))

kk=IsingOccu.fit.Moller.sampler_withZ(X=X,distM=distanceM,
                                detmat = detmat, 
                                detX=detX, 
                                Z=Zsample,
                                mcmc.save = 50000, burn.in = 1000 , 
                                vars_prior = rep(1000000,4*ncol(X)+2*ncol(detX[[1]])+5),
                                vars_prop = var_prop,
                                int_range = "nn",seed = 42
                                ,init = MPLE$par
                                , thin.by = 1)


plot(kk$theta.mcmc[,2])
