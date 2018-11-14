require('gtools')
set.seed(125)
source('autooccu 20181001ver.R')
nlat = 20
siteposi = 1.00 * permutations(n=nlat,r=2,v=(1:nlat),repeats.allowed = T)

distanceM = as.matrix((dist(siteposi)))
distanceM = distanceM-1
diag(distanceM) = 0

set.seed(12345)
#x = rep(0:(nlat - 1) / (nlat - 1), times = nlat) - 0.5
#y = rep(0:(nlat - 1) / (nlat - 1), each = nlat) - 0.5
x = runif(nlat^2,-1,1)
y = runif(nlat^2,-1,1)
ones = rep(1,times = nlat*nlat)
X = cbind(ones,x,y)

raster::plot(raster::raster(
  matrix( 
    x
    [1:400],
    nrow=20,ncol=20)))

# 2 env var and 2 det var, 5 repeats
detX = list()
nperiod = 7
for (i in 1:nperiod){
  temp = matrix(runif(nlat^2 * 2),nrow = nlat^2,ncol = 2)
	detX[[i]] = temp
}


theta = matrix(c(0.2,1,1, # env reaction of 1
                 .5,1,1,  # env reaction of 2
                 1,-1,1,-1,2,  # detection beta of 1
                 1,1,1,-1,2,   # detection beta of 2
                 0.1,2,        # eta01 d1
                 0.1,2,		  # eta02 d2
                 -1))
# first 3, all environmental factor for spc.1, 4-6, environment for spc.2, 7-11, detection for spc.1
#   12-16 detection for spc.2, 17, spatial for spc.1, 18 spatial for spc.2, 19 interspecies
detmat = matrix(0,nrow = 2*nlat^2,ncol = nperiod) # a sample detection matrix



p = length(theta)
sites = nrow(distanceM)
ncov = ncol(X)
zeros = matrix(0,nrow=400,ncol=ncov)
beta1 = as.numeric( matrix(c(theta[1:(2*ncol(X))])))
Xfull = cbind(rbind(X,zeros),rbind(zeros,X))
thr = Xfull%*%beta1

#plot(1:200,exp(thr)/((exp(-thr))+(exp(thr))))

raster::plot(raster::raster(
  matrix( 
    (exp(thr)/(exp(-thr)+exp(thr)))
     [401:800],
    nrow=20,ncol=20)))

set.seed(12345)
Zsample = rautoccu(X,distanceM,theta,method = "MH",nIter=300,n=1,int_range = "exp")

raster::plot(raster::raster(
  matrix( 
    Zsample
    [401:800],
    nrow=20,ncol=20)))
#raster::plot(raster::raster(matrix(detmat[401:800,1],nrow = 20,ncol=20)))

detSample = autooccu_sample.detection(theta, X, distM, Z=Zsample, detmat, detX)
detmat = detSample

raster::plot(raster::raster(
  matrix( 
    detSample
    [1:400,1],
    nrow=20,ncol=20)))

autooccu.logPL(theta, X, distanceM, Z=Zsample, detmat, detX,int_range = "exp")

optPLwithZ = optim(par=(theta),fn=autooccu.logPL,NULL,envX=X,distM=distanceM,Z=Zsample,detmat=detmat,detX=detX,int_range = "exp")

optPLwithZ$par
abs((theta-optPLwithZ$par)/theta)


### NOT RUN

Zsample_MCEM = autooccuMCEM_sampleZ(theta, X, A, A1, A2, 1*(runif(length(Zsample))>=0.5), detmat, detX, 10000)

cl <- parallel::makePSOCKcluster(detectCores())
ptm = proc.time()
autooccu.ElogPL(theta, X, A, A1 , A2, Zsample_MCEM, detmat, detX)
proc.time()-ptm
ptm = proc.time()
autooccu.ElogPLnotpr(theta, X, A, A1 , A2, Zsample_MCEM, detmat, detX)
proc.time()-ptm
#ptm = proc.time()
#cl <- parallel::makePSOCKcluster(4)
clusterExport(cl,list("Zsample_MCEM","theta","X","A","A1","A2","detmat","detX","autooccu.logPL"),envir = environment())
www = parallel::parCapply(cl=cl,x=Zsample_MCEM,FUN='autooccu.logPL',theta = theta,
                envX = X, A = A, A1 = A1, A2 = A2,  detmat = detmat, detX = detX)
#proc.time()-ptm

#cl <- parallel::makePSOCKcluster(detectCores())
autooccu.fit(X, A, A1, A2, detSample, detX, MCEMsample = 5000 ,hessian = F, method = 'BFGS')
#closeNode(cl)
