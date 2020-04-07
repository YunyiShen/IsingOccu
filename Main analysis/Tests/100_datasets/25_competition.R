source("./R/misc.R")
source("./R/Main_Sampler.R")
source("./R/Simu_data_Sampling.R")
require(Matrix)
require(Rcpp)
Rcpp::sourceCpp("src/IsingCpp_CFTP_sparse.cpp")


## generate graph 
n_grids = 25# 25 by 25 grid system
link_inner = adjacency.matrix(n_grids) # nearest neighborhood 
link_outer = Matrix(0,n_grids^2,n_grids^2,sparse = T)
link_mainland = matrix(0,n_grids^2,1)

###### True Parameter Setting ######

spp_mat = matrix(1,2,2)
diag(spp_mat)=0
spp_mat = as(spp_mat,'dsCMatrix')
envX = matrix(1,n_grids^2,1)
envX = cbind(envX ,rnorm(n_grids^2))

theta = list(beta_occu = c(-.3,.5,-.3,.5),
             beta_det = c(-.3,.3,-.3,.3),
             eta_intra = c(0.25,0.25),
             eta_inter = c(1,1),
             spp_mat = -.3 * spp_mat,
             spp_mat_det = -.2 * spp_mat)


link_map = 
  list(inter = link_outer,
       intra = link_inner)

nrep = 1
nspp = 2
nperiod = 6
nsite = n_grids^2
distM_mainland = matrix(0,nsite,1)

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

detmat_format = lapply(1:nrep,function(dummy,nperiod,nsite,nspp){
  matrix(-1,nrow = nsite*nspp,ncol = nperiod)
},nperiod,nsite,nspp) # just for formating, hopefully it works

n_dataset = 100

eta_intra_1 = matrix(NA,nrow = n_dataset,ncol = 5000)
eta_intra_2 = matrix(NA,nrow = n_dataset,ncol = 5000)
beta_1 = matrix(NA,nrow = n_dataset,ncol = 5000)
beta_2 = matrix(NA,nrow = n_dataset,ncol = 5000)
gamma_oc = matrix(NA,nrow = n_dataset,ncol = 5000)
gamma_de = matrix(NA,nrow = n_dataset,ncol = 5000)

###### Simulate Data ######
set.seed(42)

for(i in 1:n_dataset){
	cat(i,"\n\n")
  MRF = getMRF(theta,envX,distM = 0*link_map[[1]],link_map,link_mainland, link_mainland = link_mainland ,
			 int_range_intra="nn",int_range_inter="nn")

  Z_simu = IsingSamplerCpp(1, MRF$A, MRF$thr, 1, 30, c(-1,1), F,NA+MRF$thr) ## take true occupancy


  detmat_simu = Sample_Ising_detection_rep(nrep,nperiod,envX,NULL,
										 theta$beta_det,theta$spp_mat_det,Z_simu,
										 detmat_format,nIter=100,n=1, method = "CFTP")


  kk = IsingOccu.fit.Murray.sampler_Ising_det(X = envX, detmat =  detmat_simu
                                            , detX =  NULL
                                            , mcmc.iter = 50000, burn.in = 5000
                                            , vars_prop = vars_prop
                                            , para_prior = para_prior
                                            , Zprop_rate = 1
                                            , uni_prior = F
                                            , distM=link_map[[1]],link_map=link_map
                                            , dist_mainland =  distM_mainland , link_mainland =  link_mainland 
                                            , int_range_intra="nn",int_range_inter="nn"                                          
                                            , seed = 42
                                            , ini = theta,thin.by = 10,report.by = 2500,nIter = 50,method = "CFTP",Gibbs = T)
  eta_intra_1[i,] = kk$theta.mcmc$eta_intra[,1]
  eta_intra_2[i,] = kk$theta.mcmc$eta_intra[,2]
  beta_1[i,] = kk$theta.mcmc$beta_occu[,2]
  beta_2[i,] = kk$theta.mcmc$beta_occu[,4]
  gamma_oc[i,] = kk$theta.mcmc$spp_mat[,2]
  gamma_de[i,] = kk$theta.mcmc$spp_mat_det[,2]
}

write.csv(eta_intra_1,"./Main analysis/Results/Big_simulation/25/C/eta_intra_1.csv")
write.csv(eta_intra_2,"./Main analysis/Results/Big_simulation/25/C/eta_intra_2.csv")

write.csv(beta_1,"./Main analysis/Results/Big_simulation/25/C/beta_1.csv")
write.csv(beta_2,"./Main analysis/Results/Big_simulation/25/C/beta_2.csv")

write.csv(gamma_oc,"./Main analysis/Results/Big_simulation/25/C/gamma_oc.csv")
write.csv(gamma_de,"./Main analysis/Results/Big_simulation/25/C/gamma_de.csv")









