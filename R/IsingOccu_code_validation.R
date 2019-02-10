require('gtools')
set.seed(42)
source('IsingOccu.R')
source("misc.R")
nlat = 20
siteposi = 1.00 * permutations(n=nlat,r=2,v=(1:nlat),repeats.allowed = T)

distanceM = as.matrix((dist(siteposi)))
distanceM = distanceM-1
diag(distanceM) = 0

set.seed(42)
#x = rep(0:(nlat - 1) / (nlat - 1), times = nlat) - 0.5
#y = rep(0:(nlat - 1) / (nlat - 1), each = nlat) - 0.5
x = runif(nlat^2,-1,1)
y = runif(nlat^2,-1,1)
ones = rep(1,times = nlat*nlat)
X = cbind(ones,x,y)

raster::plot(raster::raster(
  matrix(
    x
    [1:nlat^2],
    nrow=nlat,ncol=nlat)))

# 2 env var and 2 det var, 5 repeats
detX = list()
nperiod = 5
for (i in 1:nperiod){
  temp = matrix(runif(nlat^2 * 2),nrow = nlat^2,ncol = 2)
	detX[[i]] = temp
}


theta = matrix(c(-0.35,-1,-1, # env reaction of 1
                 -.15,1,1,  # env reaction of 2
                 0,-1,1,-1,1,  # detection beta of 1
                 0,1,-1,1,-1,   # detection beta of 2
                 0.15,3,        # eta01 d1
                 0.15,3,		  # eta02 d2
                 -.1))
# first 3, all environmental factor for spc.1, 4-6, environment for spc.2, 7-11, detection for spc.1
#   12-16 detection for spc.2, 17, spatial for spc.1, 18 spatial for spc.2, 19 interspecies
detmat = matrix(0,nrow = 2*nlat^2,ncol = nperiod) # a sample detection matrix



p = length(theta)
sites = nrow(distanceM)
ncov = ncol(X)
zeros = matrix(0,nrow=nlat^2,ncol=ncov)
beta1 = as.numeric( matrix(c(theta[1:(2*ncol(X))])))
Xfull = cbind(rbind(X,zeros),rbind(zeros,X))
thr = Xfull%*%beta1

#plot(1:200,exp(thr)/((exp(-thr))+(exp(thr))))
raster::plot(raster::raster(
  matrix(
    (exp(thr)/(exp(-thr)+exp(thr)))
    [1:(nlat^2)],
    nrow=nlat,ncol=nlat)))

raster::plot(raster::raster(
  matrix(
    (exp(thr)/(exp(-thr)+exp(thr)))
     [1:nlat^2 + nlat^2],
    nrow=nlat,ncol=nlat)))

# Test the sampler Z function
set.seed(12345)
Zsample = rIsingOccu(X,distanceM,theta,method = "CFTP",nIter=100,n=1,int_range = "exp")

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

# Test the detection function
detSample = IsingOccu_sample.detection(theta, X,  Z=Zsample, detmat, detX)
detmat = detSample

raster::plot(raster::raster(
  matrix(
    detSample
    [1:nlat^2,1],
    nrow=nlat,ncol=nlat)))

# log Pseudo-Likelihood
#IsingOccu.logPL(theta, X, distanceM, Z=Zsample, detmat, detX,int_range = "exp")

optPLwithZ = optim(par=((theta)),fn=IsingOccu.logPL,NULL,envX=X,distM=distanceM,Z=Zsample,detmat=detmat,detX=detX,int_range = "exp",control=list(maxit=15000))
#optPLZ = optim(par=((theta)),fn=logPL,NULL,envX=X,distM=distanceM,Z=Zsample,int_range = "exp",control=list(maxit=5000))
logPL(theta,envX=X,distM=distanceM,Z=Zsample,int_range = "exp")
#optPLwithZ$par
#abs((theta-optPLwithZ$par)/theta)

# Test MCMC helper functions
## test negPotential
envX = X
distM=distanceM
Z=Zsample
int_range = "exp"

# theta_simple = c(0,0,100,100,0,2,0,2,-1)
# X_simple = matrix(1)
# detX_simple = list() 
# detX_simple[[1]]=NULL
# Z_simple = matrix(c(1,-1),2,1)
# detmat_simple = matrix(c(1,-1),2,1)
# distM_simple = matrix(0)
# IsingOccu.logL.innorm(theta_simple, envX=X_simple, distM=distM_simple, Z=Z_simple ,detmat_simple, detX_simple, int_range = "exp")
## Initial value
Ini = Initial_MPLE(detmat,envX,detX,distM,"exp")


IsingOccu.logL.innorm(theta, envX=X, distM=distanceM, Z=Zsample ,detmat, detX, int_range = "exp")
IsingOccu.logL.innorm(theta+runif(length(theta)), envX=X, distM=distanceM, Z=Zsample ,detmat, detX, int_range = "exp")
# since normalizing function is different, is not campairable 

## test Moller ratio
#(2*(runif(length(Z))>0.5)-1)

Moller.ratio(theta_curr=theta 
                        ,theta_prop=theta+.1*runif(length(theta))
                        #,theta_prop = theta
                        ,Z_curr=Z
                        #,Z_prop=(2*(runif(length(Z))>0.5)-1)
                        ,Z_prop = Z
                        ,Z_temp_curr = Z
                        ,Z_temp_prop = Z
                        #,x_curr=detmat
                        #,x_prop=IsingOccu_sample.detection(theta, X,  Z=Zsample, detmat, detX)
                        ,detmat=detmat
                        ,vars_prior=1
                        ,theta_tuta=theta
                        ,envX, detX, distM,int_range="exp" )

## test sampler

var_prop = c(rep(2.5e-5,6),rep(2.5e-3,10),1e-6,4e-8,1e-6,4e-8,1e-6)

kk=IsingOccu.fit.Moller.sampler(X=X,distM=distanceM,
                                detmat = detmat, 
                                detX=detX, 
                                mcmc.save = 20000, burn.in = 300 , 
                                vars_prior = rep(1000000,4*ncol(X)+2*ncol(detX[[1]])+5),
                                vars_prop = var_prop,
                                int_range = "exp",seed = 42
                                ,init = optPLwithZ$par
                                , thin.by = 1
                                )
 
plot(kk$theta.mcmc[,21])
