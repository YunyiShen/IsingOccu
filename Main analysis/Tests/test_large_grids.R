source("./R/misc.R")
source("./R/Main_Sampler.R")
source("./R/Simu_data_Sampling.R")
require(Matrix)
require(Rcpp)
Rcpp::sourceCpp("src/IsingCpp_CFTP_sparse.cpp")


## generate graph 
n_grids = 15# 15 by 15 grid system
link_inner = adjacency.matrix(n_grids) # nearest neighborhood 
link_outer = Matrix(0,n_grids^2,n_grids^2,sparse = T)
link_mainland = matrix(0,n_grids^2,1)

###### True Parameter Setting ######

spp_mat = matrix(1,2,2)
diag(spp_mat)=0
spp_mat = as(spp_mat,'dsCMatrix')
envX = matrix(1,n_grids^2,1)
envX = cbind(envX ,rnorm(n_grids^2))

theta = list(beta_occu = c(-.5,-.5,-.5,.5),
             beta_det = c(-.3,.3,-.3,.3),
             eta_intra = c(0.15,0.15),
             eta_inter = c(1,1),
             spp_mat = 0 * spp_mat,
             spp_mat_det = -.2 * spp_mat)

link_map = 
  list(inter = link_outer,
       intra = link_inner)

nrep = 1
nspp = 2
nperiod = 6
nsite = n_grids^2


distM_mainland = matrix(0,nsite,1)

###### Simulate Data ######
set.seed(42)
MRF = getMRF(theta,envX,distM = 0*link_map[[1]],link_map,link_mainland, link_mainland = link_mainland ,
			 int_range_intra="nn",int_range_inter="nn")

Z_simu = IsingSamplerCpp(1, MRF$A, MRF$thr, 1, 30, c(-1,1), F,NA+MRF$thr) ## take true occupancy

detmat_format = lapply(1:nrep,function(dummy,nperiod,nsite,nspp){
	matrix(-1,nrow = nsite*nspp,ncol = nperiod)
},nperiod,nsite,nspp) # just for formating, hopefully it works

detmat_simu = Sample_Ising_detection_rep(nrep,nperiod,envX,NULL,
										 theta$beta_det,theta$spp_mat_det,Z_simu,
										 detmat_format,nIter=100,n=1, method = "CFTP")

#year1_dethis = detmat_simu[[1]]
Z_absolute = (sapply(detmat_simu,function(detmat_i){rowSums((detmat_i+1)/2)>0})) * 2 - 1


###### Run the Model! ######

vars_prop = list( beta_occu = c(1e-3,1e-3)
                  ,beta_det = rep(5e-3,nspp * ( ncol(envX)) ) # no extra det thing
                  ,eta_intra = rep(1e-3,nspp)
                  ,eta_inter = c(2e-3,2e-3)
                  ,spp_mat = 1e-3
                  ,spp_mat_det = 5e-3)

para_prior = list( beta_occu = rep(1000,nspp * ncol(envX))
                   ,beta_det = rep(10,nspp * (ncol(envX)) )
                   ,eta_intra = rep(1000,nspp)
                   ,eta_inter = rep(1000,nspp)
                   ,d_intra=rep(1000,nspp)
                   ,d_inter = rep(1000,nspp)
                   ,spp_mat = 1000
                   ,spp_mat_det = 1000)


kk = IsingOccu.fit.Murray.sampler_Ising_det(X = envX, detmat =  detmat_simu
                                            , detX =  NULL
                                            , mcmc.iter = 150000, burn.in = 5000
                                            , vars_prop = vars_prop
                                            , para_prior = para_prior
                                            , Zprop_rate = 1
                                            , uni_prior = F
                                            , distM=link_map[[1]],link_map=link_map
                                            , dist_mainland =  distM_mainland , link_mainland =  link_mainland 
                                            , int_range_intra="nn",int_range_inter="nn"                                          
                                            , seed = 42
                                            , ini = theta,thin.by = 1,report.by = 500,nIter = 50,method = "CFTP",Gibbs = T)


save.image("Test_large_grid20by20_niche_diff_100K_with_intra.RData")
