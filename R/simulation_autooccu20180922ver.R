require('gtools')
set.seed(125)
source('autooccu 20180922ver.R')
nlat = 5
siteposi = 1.00 * permutations(n=nlat,r=2,v=(1:nlat),repeats.allowed = T)

distanceM = exp(-as.matrix(dist(siteposi)))
distanceM = distanceM/max(distanceM)
diag(distanceM) = 0
zeros = matrix(0,nrow = nlat^2,ncol = nlat^2)
eyes = diag(1,nrow = nlat^2,ncol = nlat^2)
A = rbind(cbind(distanceM,zeros),cbind(zeros,zeros))
A1 = rbind(cbind(zeros,zeros),cbind(zeros,distanceM))
A2 = rbind(cbind(zeros,eyes),cbind(eyes,zeros))
rm(zeros)
rm(eyes)

x = rep(0:(nlat - 1) / (nlat - 1), times = nlat) - 0.5
y = rep(0:(nlat - 1) / (nlat - 1), each = nlat) - 0.5
ones = rep(1,times = nlat*nlat)
X1 = cbind(ones,x,y)
X2 = cbind(X1,matrix(0,nrow(X1),ncol(X1)))
X3 = cbind(matrix(0,nrow(X1),ncol(X1)),X1)
X = rbind(X2,X3)
rm(X1)
rm(X2)
rm(X3)
# 2 env var and 2 det var, 5 repeats
detX = list()
nperiod = 5
for (i in 1:nperiod){
  temp = matrix(runif(nlat^2 * 2),nrow = nlat^2,ncol = 2)
	detX[[i]] = rbind(cbind(temp,0*temp),cbind(0*temp,temp))
}

theta = matrix(c(-2,2,2,-1,1,4,1,-1,1,-1,2,1,1,1,-1,2,0.6,0.4,-0.5)) # real theta
# first 3, all environmental factor for spc.1, 4-6, environment for spc.2, 7-11, detection for spc.1
#   12-16 detection for spc.2, 17, spatial for spc.1, 18 spatial for spc.2, 19 interspecies
detmat = matrix(0,nrow = 2*nlat^2,ncol = nperiod) # a sample detection matrix

Zsample = rautoccu(X, A, A1, A2,theta)

detSample = autooccu_sample.detection(theta, X, A, A1, A2, Zsample, detmat, detX)
detmat = detSample

autooccu.logPL(theta, X, A, A1, A2, Zsample, detSample, detX)

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
