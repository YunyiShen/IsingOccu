source("misc_island.R")
#source("Moller_island.R")
source("Murray_Ising_det.R")

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

theta = list(beta_occu = c(-.25,-.25),
             beta_det = c(-.3,-.3),
             eta_intra = c(.3,.3),
             eta_inter = c(.3,.3),
             #d_inter = c(.2,.2),
             spp_mat = .1* spp_mat,
             spp_mat_det = -.3 * spp_mat)

link_map = 
  list(inter = link_outer * exp(-distM_full),
       intra = link_inner)

nrep = 1
set.seed(12345)
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

  #    then in the second level list, it is a matrix with nrow = site ncol = ncov, 

#Pdet = Pdet_multi(nperiod, envX,detX[[1]], theta$beta_det, nspp=nrow(spp_mat))

#detmat = Sample_detection(nrep,nperiod,envX,detX,theta$beta_det,nspp = nrow(spp_mat),Z=Z_sample)
#detmat = lapply(detmat,function(w){w*2-1}) 

#sppmat_det = -0.1 * spp_mat
#Pdet_Ising(nperiod,envX,detX[[1]],beta_det = theta$beta_det,sppmat_det,Z = Z_sample,detmat[[1]])

nperiod = 10
nsite = 155
nspp = 2
detmat = list(matrix(-1,nsite*nspp,nperiod))

no_obs=11
set.seed(42)
detmat = Sample_Ising_detection_rep(1,nperiod,envX,detX=NULL,beta_det = theta$beta_det,theta$spp_mat_det,Z = ( Z_sample),detmat,nIter=100,n=1, method = "CFTP")
detmat[[1]][1]=NA # get some NAs
detmat[[1]][c(2,155+2),]=NA
Pdet_Ising_rep(1,nperiod,envX,detX=NULL,beta_det = theta$beta_det,theta$spp_mat_det,Z = Z_sample,detmat)

#detmat[[1]][c(no_obs,no_obs+155),]=-1

tempdata = data.frame(island[,6:7],
                      Z_1 = detmat[[1]][1:155,2],
                      Z_2 = detmat[[1]][156:310,2])

ggplot(data = tempdata,aes(x=X,y=Y,color = Z_1))+
  geom_point()


nspp = 2

vars_prop = list( beta_occu = rep(1e-4,nspp * ncol(envX))
                  ,beta_det = rep(2.5e-3,nspp * ( ncol(envX)) ) # no extra det thing
                  ,eta_intra = rep(1e-4,nspp)
                  ,eta_inter = rep(1e-4,nspp)
                  #,d_intra=rep(2.5e-5,nspp)
                  #,d_inter = rep(1e-4,nspp)
                  ,spp_mat = 1e-4
                  ,spp_mat_det = 1e-3)


Z_absolute = (sapply(detmat,function(detmat_i){rowSums((detmat_i+1)/2)>0})) * 2 - 1


no_obs = c(no_obs,no_obs+155)


kk = IsingOccu.fit.Murray.sampler_Ising_det(X = envX, detmat =  detmat,no_obs = NULL
                                  , detX =  NULL
                                  , mcmc.iter = 5000, burn.in = 500
                                  , vars_prop = vars_prop
                                  , vars_prior = 200000
                                  , Zprop_rate = .8
                                  , Zprop_rate_missing_obs = 0
                                  , distM=distM_full,link_map=link_map
                                  , dist_mainland =  distM_mainland , link_mainland =  link_mainland * exp(-distM_mainland)
                                  , int_range_intra="nn",int_range_inter="nn"
                                  #, Z = Z_sample # just used in formating, if assuming perfect detection, simple giving Z and set Zprop_rate=0
                                  , Z = Z_absolute
                                  , seed = 42
                                  , ini = theta,thin.by = 1,report.by = 500)







