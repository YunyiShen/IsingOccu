source("./R/misc.R")
source("./R/Main_Sampler.R")
source("./R/Simu_data_Sampling.R")
require(Matrix)
require(Rcpp)
require(RcppArmadillo)
library(lineprof)
Rcpp::sourceCpp("src/IsingCpp_CFTP_sparse.cpp")

# test Pdet_Ising_single_siteCpp:
thr <- matrix(rnorm(10),5,2)
Z=matrix(1,2,1)
dethis <- 2*(matrix(runif(10),5,2)<0.5)+1
dethis[1,]=NA

sppmat_det <- matrix(c(0,0.3,0.3,0),2,2)

Pdet_Ising_single_siteCpp(thr,Z,dethis,sppmat_det,c(-1L,1L))
Pdet_Ising_single_site(thr,Z,dethis,sppmat_det)

lineprof(lapply(1:1000,function(dummy){ Pdet_Ising_single_siteCpp(thr,Z,dethis,sppmat_det,c(-1L,1L))}))
lineprof(lapply(1:1000,function(dummy){ Pdet_Ising_single_site(thr,Z,dethis,sppmat_det)}))

# test extract_thrCpp
set.seed(42)
thr_list <- list()
thr_list[[1]] <- matrix(runif(6),3,2)
thr_list[[2]] <- matrix(runif(6),3,2)

extract_thr(1,thr_list)
extract_thrCpp(0,thr_list,2,2,3)

lineprof(lapply(1:10000,function(dummy) extract_thr(1,thr_list)))
lineprof(lapply(1:10000,function(dummy) extract_thrCpp(0,thr_list,2,2,3)))

## test Gibbs helper
set.seed(42)
spp_mat = matrix(1,2,2)
diag(spp_mat)=0
spp_mat = as(spp_mat,'dsCMatrix')
envX = matrix(1,n_grids^2,1)
envX = cbind(envX ,rnorm(n_grids^2))

theta = list(beta_occu = c(-.5,-.5,-.5,.5),
             beta_det = c(-.3,.5,-.3,.5),
             eta_intra = c(0.1,0.1),
             eta_inter = c(1,1),
             spp_mat = 0 * spp_mat,
             spp_mat_det = -.2 * spp_mat)

link_map = 
  list(inter = link_outer,
       intra = link_inner)

nrep = 1
nspp = 2
nperiod = 5
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


det_thr = list(matrix(rnorm(nsite*nperiod),nsite,nperiod),
               matrix(rnorm(nsite*nperiod),nsite,nperiod))





Z_new = Gibbs_Z_helperCpp(Z_absolute,which(Z_absolute==-1),detmat_simu[[1]],MRF$A,MRF$thr,as.matrix(spp_mat),det_thr,nsite,c(-1L,1L))
