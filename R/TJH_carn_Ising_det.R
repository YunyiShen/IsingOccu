source("misc_island.R")
#source("Moller_island.R")
source("Murray_Ising_det.R")

###### graph data ######
link = "C:/Users/yshen99/Documents/GitHub/RFIBM_MCMC/TJH/"
link = "E:/UW Lab jobs/2. ISING Occupancy model/2. RFIBMs MCMC/RFIBM/TJH/"
island = read.csv(paste0(link,"TJH_unique_grids.csv"))

link_outer = as.matrix( read.csv(paste0(link, "link_TJH.csv")))
#link_inner = as.matrix( read.csv("link_outer.csv",row.names = 1))
link_inner = 0 * link_outer # this makes it a mainland-island system
link_mainland = matrix(0,97,1)

distM_full = as.matrix( read.csv(paste0(link,"distM_TJH.csv")))
distM_mainland = matrix(100,97,1)

intcd = min(min((link_outer*distM_full)[(link_outer*distM_full)>0]))
normd = max(max(link_outer*distM_full))-intcd


distM_full = (distM_full-intcd)/normd # normalizing the distance
distM_mainland = (distM_mainland-intcd)/normd

detmat = list(as.matrix(read.csv(paste0(link,"BM_GZL_HHD_ZH_20dayfull.csv")))[98:(97*4),])
#full = read.csv(paste0(link,"PA_all_full.csv"),row.names=1)
#Z_sample = matrix(c(full$Coyote,full$Fox_red,full$Bobcat))

###### simulation ######

spp_mat = matrix(1,3,3)
diag(spp_mat)=0
envX = cbind(matrix(1,97,1),island$ELE)
envX = apply(envX,2,function(k){(k-min(k))/(max(k)-min(k))})
envX[,1]=1

theta = list(beta_occu = rep(0,6),
             beta_det = rep(0,6),
             eta_intra = c(0,0,0),
             eta_inter = c(.1,.1,.1),
             #d_inter = c(.2,.2),
             spp_mat = 0.1 * spp_mat,
             spp_mat_det = -0.1 * spp_mat)

link_map = 
  list(inter = link_outer, # * exp(-distM_full),
       intra = link_inner)

nrep = 1
  #    then in the second level list, it is a matrix with nrow = site ncol = ncov, 

#Pdet = Pdet_multi(nperiod, envX,detX[[1]], theta$beta_det, nspp=nrow(spp_mat))

#detmat = Sample_detection(nrep,nperiod,envX,detX,theta$beta_det,nspp = nrow(spp_mat),Z=Z_sample)
#detmat = lapply(detmat,function(w){w*2-1}) 

#sppmat_det = -0.1 * spp_mat
#Pdet_Ising(nperiod,envX,detX[[1]],beta_det = theta$beta_det,theta$sppmat_det,Z = Z_sample,detmat[[1]])
#Pdet_Ising_rep(1,52,envX,NULL,theta$beta_det,theta$spp_mat_det,Z = Z_sample,detmat)

#no_obs=150:155
#no_obs = c(no_obs, no_obs + 155, no_obs + 310)

nspp = 3

vars_prop = list( beta_occu = rep(3e-3,nspp * ncol(envX))
                  ,beta_det = rep(2.5e-3,nspp * ( ncol(envX)) ) # no extra det thing
                  ,eta_intra = rep(2e-4,nspp)
                  ,eta_inter = rep(5e-4,nspp)
                  #,d_intra=rep(2.5e-5,nspp)
                  #,d_inter = rep(1e-4,nspp)
                  ,spp_mat = 2.5e-3
                  ,spp_mat_det = 3e-3)

detmat_nona = lapply(detmat,function(mat){
  mat[is.na(mat)]=-1
  return(mat)
})
Z_absolute = (sapply(detmat_nona,function(detmat_i){rowSums((detmat_i+1)/2)>0})) * 2 - 1

rm(detmat_nona)
datatemp  = data.frame(island,
                       Z1 = Z_absolute[1:97,],
                       Z2 = Z_absolute[1:97+97],
                       Z3 = Z_absolute[1:97+2*97],
                       Z4 = Z_absolute[1:97+3*97])

#Z_absolute = (sapply(detmat,function(detmat_i){rowSums((detmat_i+1)/2)>0})) * 2 - 1
require(ggplot2)
ggplot(data = datatemp,aes(x=LONG,y=LAT,color = Z1))+
  geom_point()



kk = IsingOccu.fit.Murray.sampler_Ising_det(X = envX, detmat =  detmat
                                  , detX =  NULL
                                  , mcmc.iter = 20000, burn.in = 1500
                                  , vars_prop = vars_prop
                                  , vars_prior = 200000
                                  , Zprop_rate = 0.05
                                  #, Zprop_rate_missing_obs = 0
                                  , distM=distM_full,link_map=link_map
                                  , dist_mainland =  distM_mainland , link_mainland =  link_mainland * exp(-distM_mainland)
                                  , int_range_intra="nn",int_range_inter="nn"
                                  #, Z = Z_sample # just used in formating, if assuming perfect detection, simple giving Z and set Zprop_rate=0
                                  #, Z = Z_absolute
                                  , seed = 42
                                  , ini = theta,thin.by = 5,report.by = 100,nIter = 30)







