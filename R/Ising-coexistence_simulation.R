require('gtools')
require("IsingSampler")
require("raster")
set.seed(125)
source('autooccu 20181001ver.R')
nlat = 20
siteposi = 1.00 * permutations(n=nlat,r=2,v=(1:nlat),repeats.allowed = T)

distanceM = as.matrix((dist(siteposi)))
distanceM = distanceM-1
diag(distanceM) = 0

#neutral
outfield_neu = matrix(0,nrow=2*nlat^2,ncol=1)
theta = c(0.1,1,        # eta01 d1
          0.1,1,		  # eta02 d2
          -1
) #eta1

Graph = getGraph(distanceM,theta)
Zsample = IsingSampler(n=3000,graph = Graph,thresholds=outfield_neu, responses = c(-1L, 1L),nIter=100,method="CFTP")
cosine = rowSums(Zsample*Zsample)/(nlat^2)

# different half case
outfield1_half = matrix(-1,nrow = nlat,ncol = nlat)
outfield1_half[,nlat+1:nlat] = -outfield1_half[,nlat+1:nlat]
plot(raster(outfield1))
outfield1_half = reshape(outfield1_half,nrow=nlat^2,ncol=1,byrow=T)
outfield2_half = -outfield1_half
outfield_half=rbind(outfield1_half,outfield2_half)

# random case
outfield_rand = matrix(runif(nlat^2,min=-1,max=1),nrow=2*nlat^2,ncol=1)

