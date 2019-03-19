source("misc_island.R")
source("Moller_island.R")

###### graph data ######
link = "C:/Users/yshen99/Documents/GitHub/RFIBM_MCMC/Island/"
#link = "E:/UW Lab jobs/2. ISING Occupancy model/2. RFIBMs MCMC/RFIBM/island/"
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
             eta_inter = c(.3,.3),
             d_inter = c(.2,.2),
             spp_mat = -0.1 * spp_mat)

link_map = 
  list(inter = link_outer,
       intra = link_inner)

nrep = 4
set.seed(42)
Z_sample = rIsingOccu_multi(theta,
                            envX,
                            distM_full,link_map ,
                            distM_mainland,link_mainland,
                            int_range_intra="nn",int_range_inter="exp",
                            n=nrep,method = "CFTP",nIter = 100)

require(ggplot2)

tempdata = data.frame(island[,6:7],
                      Z_1 = Z_sample[1:155,1],
                      Z_2 = Z_sample[156:310,1])

ggplot(data = tempdata,aes(x=X,y=Y,color = Z_1))+
  geom_point()

Hamiltonian(theta,envX,distM_full,link_map,distM_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp",Z_sample)


detX = list()
nperiod = 10
for(j in 1:nrep){
  detX[[j]] = list()
  for (i in 1:nperiod){
    temp = matrix(runif(155 * 2),nrow = 155,ncol = 2)
    detX[[j]][[i]] = temp
  }
}

Pdet = Pdet_multi(nperiod, envX,detX[[1]], theta$beta_det, nspp=nrow(spp_mat))

detmat = Sample_detection(nrep,nperiod,envX,detX,theta$beta_det,nspp = nrow(spp_mat),Z=Z_sample)


tempdata = data.frame(island[,6:7],
                      Z_1 = detmat[[1]][1:155,1],
                      Z_2 = detmat[[1]][156:310,1])

ggplot(data = tempdata,aes(x=X,y=Y,color = Z_1))+
  geom_point()



H = IsingOccu_multi.logL.innorm(theta, envX, distM_full,link_map,distM_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp", Z_sample ,detmat, detX)

theta_prop = list(beta_occu = c(0.01,0.01),
                  beta_det = c(0,1,-1,0,1,-1),
                  eta_intra = c(.15,.151),
                  eta_inter = c(.31,.3),
                  d_inter = c(.21,.2),
                  spp_mat = -0.12 * spp_mat)
Z_temp_prop = rIsingOccu_multi(theta_prop,
                               envX,
                               distM_full,link_map ,
                               distM_mainland,link_mainland,
                               int_range_intra="nn",int_range_inter="exp",
                               n=nrep,method = "CFTP",nIter = 100)
M_ratio = Moller.ratio(theta_curr = theta ,theta_prop
                        ,Z_curr = Z_sample ,Z_prop = Z_sample
                        ,Z_temp_curr = Z_sample, Z_temp_prop
                        ,detmat
                        ,vars_prior = 10000
                        ,theta_tuta = theta
                        ,envX, detX
                        ,distM_full,link_map
                        ,distM_mainland,link_mainland
                        ,int_range_intra="nn",int_range_inter="exp")

nspp = 2
vars_prop = list( beta_occu = rep(1e-5,2 * ncol(envX))
    ,beta_det = rep(1e-5,2 * (ncol(detX[[1]][[1]]) + ncol(envX)) )
    ,eta_intra = rep(1e-5,nspp)
    ,eta_inter = rep(1e-5,nspp*(nspp-1)/2)
    ,d_intra=rep(1e-5,nspp)
    ,d_inter = rep(1e-5,nspp)
    ,spp_mat = 1e-5)


kk = IsingOccu.fit.Moller.sampler(envX,detmat,detX
                              , mcmc.save = 2000, burn.in = 100 
                              , vars_prop = vars_prop
                              , vars_prior = 2000
                              , Zprop_rate = 0
                              , distM=distM_full,link_map=link_map
                              , distM_mainland , link_mainland
                              , int_range_intra="nn",int_range_inter="exp"
                              , Z = Z_sample
                              , seed = 42
                              , ini = theta)
