source("misc_island.R")
#source("Moller_island.R")
source("Murray_island.R")

###### graph data ######
link = "C:/Users/yshen99/Documents/GitHub/RFIBM_MCMC/Island/"
link = "E:/UW Lab jobs/2. ISING Occupancy model/2. RFIBMs MCMC/RFIBM/island/"
island = read.csv(paste0(link,"CT_posi_only_island.csv"))

link_inner = as.matrix( read.csv(paste0(link, "link_inner.csv"),row.names = 1))
#link_outer = as.matrix( read.csv("link_outer.csv",row.names = 1))
link_outer = as.matrix( read.csv(paste0(link,"link_outer_full.csv"),row.names = 1))
link_mainland = as.matrix( read.csv(paste0(link,"link_mainland.csv")))
#link_outer = 0 * link_outer # this makes it a mainland-island system
#link_mainland = matrix(0,155,1)

distM_full = as.matrix( read.csv(paste0(link,"distM_full.csv"),row.names = 1))
distM_mainland = as.matrix( read.csv(paste0(link,"dist_to_mainland.csv"),row.names = 1))

intcd = min(min((distM_mainland*link_mainland)[(distM_mainland*link_mainland)>0]),
            min((link_outer*distM_full)[(link_outer*distM_full)>0]))
normd = max(max(distM_mainland*link_mainland),max(link_outer*distM_full))-intcd


distM_full = (distM_full-intcd)/normd # normalizing the distance
distM_mainland = (distM_mainland-intcd)/normd

###### simulation ######

spp_mat = matrix(c(0,1,1,0),2,2)
envX = matrix(1,155,1)

theta = list(beta_occu = c(0,0),
             beta_det = c(0,1,-1,0,1,-1),
             eta_intra = c(.15,.15),
             eta_inter = c(.2,.2),
             #d_inter = c(.2,.2),
             spp_mat = -0.15 * spp_mat)

link_map = 
  list(inter = link_outer * exp(-distM_full),
       intra = link_inner)

nrep = 4
set.seed(42)
Z_sample = rIsingOccu_multi(theta,
                            envX,
                            distM_full,link_map ,
                            distM_mainland,link_mainland * exp(-distM_mainland) ,
                            int_range_intra="nn",int_range_inter="nn",
                            n=nrep,method = "CFTP",nIter = 100)

require(ggplot2)

tempdata = data.frame(island[,6:7],
                      Z_1 = Z_sample[1:155,1],
                      Z_2 = Z_sample[156:310,1])

ggplot(data = tempdata,aes(x=X,y=Y,color = Z_2))+
  geom_point()

Hamiltonian(theta,envX,distM_full,link_map,distM_mainland,link_mainland*exp(-distM_mainland),int_range_intra="nn",int_range_inter="nn",(Z_sample))


detX = list()
nperiod = 5
for(j in 1:nrep){
  detX[[j]] = list()
  for (i in 1:nperiod){
    temp = matrix(runif(155 * 2),nrow = 155,ncol = 2)
    detX[[j]][[i]] = temp
  }
} # detX should be a list with length of # of repeats, with each element has a list has length of nperiod 
  #    then in the second level list, it is a matrix with nrow = site ncol = ncov, 

Pdet = Pdet_multi(nperiod, envX,detX[[1]], theta$beta_det, nspp=nrow(spp_mat))

detmat = Sample_detection(nrep,nperiod,envX,detX,theta$beta_det,nspp = nrow(spp_mat),Z=Z_sample)
# detmat should be a list with length nrep, each element is a matrix with nrow=nspp*nsite, ncol = nperiod
#  if there is no detection history for a site, simple give all 0, and mark the site in no_obs as 1
#  makesure the order of species should be same as spp_mat, say if 1:nsite is 1st species, in spp_mat the first row
#   (certainly also first col) should be the same species

tempdata = data.frame(island[,6:7],
                      Z_1 = detmat[[1]][1:155,2],
                      Z_2 = detmat[[1]][156:310,2])

ggplot(data = tempdata,aes(x=X,y=Y,color = Z_1))+
  geom_point()



H = IsingOccu_multi.logL.innorm(theta, envX, distM_full,link_map,distM_mainland,link_mainland,int_range_intra="nn",int_range_inter="nn", Z_sample ,detmat, detX)


nspp = 2

vars_prop = list( beta_occu = rep(1e-4,nspp * ncol(envX))
                  ,beta_det = rep(2.5e-3,2 * (ncol(detX[[1]][[1]]) + ncol(envX)) )
                  ,eta_intra = rep(1e-4,nspp)
                  ,eta_inter = rep(1e-4,nspp)
                  #,d_intra=rep(2.5e-5,nspp)
                  #,d_inter = rep(1e-4,nspp)
                  ,spp_mat = 1e-4)


Z_absolute = (sapply(detmat,rowSums)>0) * 2 - 1

no_obs = 0*Z_sample
no_obs[11]=1


kk = IsingOccu.fit.Murray.sampler(X = envX, detmat =  detmat,no_obs = no_obs
                                  , detX =  detX
                                  , mcmc.iter = 5e3, burn.in = 2e2
                                  , vars_prop = vars_prop
                                  , vars_prior = 200000
                                  , Zprop_rate = 0.1
                                  , Zprop_rate_missing_obs = 0.3
                                  , distM=distM_full,link_map=link_map
                                  , dist_mainland =  distM_mainland , link_mainland =  link_mainland * exp(-distM_mainland)
                                  , int_range_intra="nn",int_range_inter="nn"
                                  #, Z = Z_sample # just used in formating, if assuming perfect detection, simple giving Z and set Zprop_rate=0
                                  , Z = Z_absolute
                                  , seed = 42
                                  , ini = theta,thin.by = 1,report.by = 100)




