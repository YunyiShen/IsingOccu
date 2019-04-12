source("misc_island.R")
#source("Moller_island.R")
source("Murray_island.R")

###### graph data ######
link = "C:/Users/yshen99/Documents/GitHub/RFIBM_MCMC/Island/"
link = "E:/UW Lab jobs/2. ISING Occupancy model/2. RFIBMs MCMC/RFIBM/island/"
island = read.csv(paste0(link,"CT_posi_only_island.csv"))
squ = read.csv(paste0(link,"Squirrel2016.csv"),row.names = 1)+read.csv(paste0(link,"Squirrel2015.csv"),row.names = 1)
squ = apply(as.matrix(squ),1,function(x){min(x,70)})

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

spp_mat = matrix(1,2,2)
diag(spp_mat)=0
envX = matrix(1,155,1)
envX = cbind(envX,as.matrix(full$Squirrel))

theta = list(beta_occu = c(0,0.1,0,0.1),
             beta_det = c(0,0,1,-1,0,0,1,-1),
             eta_intra = c(.15,.15),
             eta_inter = c(.15,.15),
             #d_inter = c(.2,.2),
             spp_mat = -0.15 * spp_mat)

link_map = 
  list(inter = link_outer * exp(-distM_full),
       intra = link_inner)

nrep = 1
set.seed(42)

#rep1 = read.csv(paste0(link,"PA_Coyote_Foxes_1415.csv"),row.names=1)
#rep1_vec = matrix(c(rep1$Coyote,rep1$Fox_red))
#rep2 = read.csv(paste0(link,"PA_Coyote_Foxes_1617.csv"),row.names=1)
#rep2_vec = matrix(c(rep2$Coyote,rep2$Fox_red))

#Z_sample = cbind(rep1_vec,rep2_vec)

full = read.csv(paste0(link,"PA_all_full.csv"),row.names=1)
Z_sample = matrix(c(full$Fisher,full$Marten))

require(ggplot2)

tempdata = data.frame(island[,6:7],
                      Z_1 = full$Coyote,
                      Z_2 = full$Fox_red,
                      Z_3 = squ
                      )


ggplot(data = tempdata,aes(x=X,y=Y,color = log(Z_3+1)))+
  geom_point()

detX = list()
nperiod = 3
for(j in 1:nrep){
  detX[[j]] = list()
  for (i in 1:nperiod){
    temp = matrix(runif(155 * 2),nrow = 155,ncol = 2)
    detX[[j]][[i]] = temp
  }
} # detX should be a list with length of # of repeats, with each element has a list has length of nperiod 
  #    then in the second level list, it is a matrix with nrow = site ncol = ncov, 
detmat = Sample_detection(nrep,nperiod,envX,detX,theta$beta_det,nspp = nrow(spp_mat),Z=Z_sample)

nspp = 2
vars_prop = list( beta_occu = rep(2.5e-3,nspp * ncol(envX))
                  ,beta_det = rep(2.5e-3,2 * (ncol(detX[[1]][[1]]) + ncol(envX)) )
                  ,eta_intra = rep(5e-4,nspp)
                  ,eta_inter = rep(5e-4,nspp)
                  #,d_intra=rep(2.5e-5,nspp)
                  #,d_inter = rep(1e-4,nspp)
                  ,spp_mat = 5e-4)

no_obs = 0*Z_sample
no_obs[c(150:155,150:155+155),]=1				  
				  
kk = IsingOccu.fit.Murray.sampler(X = envX, detmat =  detmat,no_obs = no_obs
                                  , detX =  detX
                                  , mcmc.iter = 5e4, burn.in = 4e3
                                  , vars_prop = vars_prop
                                  , vars_prior = 200000
                                  , Zprop_rate = 0
                                  , Zprop_rate_missing_obs = 0.3
                                  , distM=distM_full,link_map=link_map
                                  , dist_mainland =  distM_mainland , link_mainland =  link_mainland * exp(-distM_mainland)
                                  , int_range_intra="nn",int_range_inter="nn"
                                  , Z = Z_sample # just used in formating, if assuming perfect detection, simple giving Z and set Zprop_rate=0
                                  #, Z = Z_absolute
                                  , seed = 42
                                  , ini = theta,thin.by = 10,report.by = 100,nIter = 100)


