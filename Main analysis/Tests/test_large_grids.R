source("./R/misc.R")
source("./R/Main_Sampler.R")
source("./R/Simu_data_Sampling.R")
require(Matrix)
require(Rcpp)
Rcpp::sourceCpp("src/IsingCpp_CFTP_sparse.cpp")


## generate graph 
n_grids = 15 # 15 by 15 grid system
link_inner = adjacency.matrix(n_grids) # nearest neighborhood 
link_outer = matrix(0,n_grids^2,1)

###### True Parameter Setting ######

spp_mat = matrix(1,2,2)
diag(spp_mat)=0
spp_mat = as(spp_mat,'dsCMatrix')
envX = matrix(1,n_grids^2,1)
envX = cbind(envX ,rnorm(n_grids^2))

theta = list(beta_occu = c(-.5,.5,-.5,-.5),
             beta_det = c(-.3,-.3),
             eta_intra = c(0.1,0.1),
             eta_inter = c(1,1),
             spp_mat = -.4 * spp_mat,
             spp_mat_det = -.2 * spp_mat)

link_map = 
  list(inter = link_outer,
       intra = link_inner)

nrep = 1
nspp = 2
nperiod = 8
nsite = n_grid^2

