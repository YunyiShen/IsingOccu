source("./R/misc.R")
#source("Moller_island.R")
source("./R/Main_Sampler.R")
require(Matrix)
require(Rcpp)
Rcpp::sourceCpp("src/IsingCpp_CFTP_sparse.cpp")

###### graph data ######
link = "./data/APIS/"
island = read.csv(paste0(link,"CT_posi_only_island.csv"))

link_inner = as.matrix( read.csv(paste0(link, "link_inner.csv"),row.names = 1))
link_inner = as(link_inner,'dsCMatrix')
link_outer = as.matrix( read.csv(paste0(link,"link_outer_full.csv"),row.names = 1))
link_outer = 0 * link_outer # this makes it a mainland-island system
link_outer = as(link_outer,'dsCMatrix')
link_mainland = matrix(1,155,1)

distM_full = as.matrix( read.csv(paste0(link,"distM_full.csv"),row.names = 1))
distM_mainland = as.matrix( read.csv(paste0(link,"dist_to_mainland.csv"),row.names = 1))

intcd = min(min((distM_mainland*link_mainland)[(distM_mainland*link_mainland)>0]),
            min((link_outer*distM_full)[(link_outer*distM_full)>0]))
normd = max(max(distM_mainland*link_mainland),max(link_outer*distM_full))-intcd


distM_full = (distM_full-intcd)/normd # normalizing the distance
distM_mainland = (distM_mainland-intcd)/normd

detmat = list(as.matrix(read.csv(paste0(link,"Coyote_Fox_Bobcat_90dfull_by_islands.csv"),header = F)[1:310,]))
full = read.csv(paste0(link,"PA_all_full.csv"),row.names=1)
Z_sample = matrix(c(full$Coyote,full$Fox_red))

###### analysis ######

spp_mat = matrix(1,2,2)
diag(spp_mat)=0
spp_mat = as(spp_mat,'dsCMatrix')
envX = matrix(1,155,1)

theta = list(beta_occu = c(-.3,-.3),
             beta_det = c(-.3,-.3),
             eta_intra = c(0,0),
             eta_inter = c(.0,.2),
             #d_inter = c(.2,.2),
             spp_mat = 0.3 * spp_mat,
             spp_mat_det = -0.3 * spp_mat)

link_map = 
  list(inter = link_outer * exp(-distM_full),
       intra = link_inner)

nrep = 1
nspp = 2

vars_prop = list( beta_occu = c(2.5e-3,2.5e-3)
                  ,beta_det = rep(1e-2,nspp * ( ncol(envX)) ) # no extra det thing
                  ,eta_intra = c(2.5e-3,2.5e-3)
                  ,eta_inter = c(5e-3,5e-3)
                  ,d_intra=rep(2.5e-5,nspp)
                  ,d_inter = rep(2.5e-3,nspp)
                  ,spp_mat = 5e-3
                  ,spp_mat_det = 5e-3)
detX = NULL

para_prior = list( beta_occu = rep(1000,2 * ncol(envX))
                   ,beta_det = rep(.1,2 * (ncol(envX)) )
                   ,eta_intra = rep(1000,nspp)
                   ,eta_inter = rep(1000,nspp)
                   ,d_intra=rep(1000,nspp)
                   ,d_inter = rep(1000,nspp)
                   ,spp_mat = 1000
                   ,spp_mat_det = 1000)

detmat_0 = lapply(detmat,function(ww){ww[is.na(ww)]=-1;return(ww)})
Z_absolute = (sapply(detmat_0,function(detmat_i){rowSums((detmat_i+1)/2)>0})) * 2 - 1





kk = IsingOccu.fit.Murray.sampler_Ising_det(X = envX, detmat =  detmat
                                  , detX =  NULL
                                  , mcmc.iter = 3000000, burn.in = 50000
                                  , vars_prop = vars_prop
                                  , para_prior = para_prior
                                  , Zprop_rate = 1
                                  , uni_prior = F
                                  , distM=distM_full,link_map=link_map
                                  , dist_mainland =  distM_mainland , link_mainland =  link_mainland * exp(-2*distM_mainland)
                                  , int_range_intra="nn",int_range_inter="nn"
                                  
                                  , seed = 42
                                  , ini = theta,thin.by = 300,report.by = 5000
                                  , nIter = 150,method = "MH",Gibbs = T)


save.image("CF_Mainland_island_3000k_norm_prior_high_det_0.1_90d_MH_Gibbs.RData")
# latest tuned parameter in 20191112



