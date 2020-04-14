source("./R/misc.R")
source("./R/Main_Sampler.R")
source("./R/Simu_data_Sampling.R")
require(Matrix)
require(Rcpp)
Rcpp::sourceCpp("src/IsingCpp_CFTP_sparse.cpp")



###### Graph ######



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


###### True Parameter Setting ######

spp_mat = matrix(1,2,2)
diag(spp_mat)=0
spp_mat = as(spp_mat,'dsCMatrix')
envX = matrix(1,155,1)

theta = list(beta_occu = c(-.3,-.3),
             beta_det = c(-.2,-.2),
             eta_intra = c(0.2,0.2),
             eta_inter = c(.3,.3),
             spp_mat = 0.25 * spp_mat,
             spp_mat_det = 0.2 * spp_mat)


link_map = 
  list(inter = link_outer * exp(-2*distM_full),
       intra = link_inner)

nrep = 1
nspp = 2
nperiod = 6
nsite = 155





vars_prop = list( beta_occu = c(5e-3,5e-3)
                  ,beta_det = rep(1e-3,nspp * ( ncol(envX)) ) # no extra det thing
                  ,eta_intra = rep(1e-3,nspp)
                  ,eta_inter = c(5e-3,5e-3)
                  ,spp_mat = 1e-2
                  ,spp_mat_det = 1e-2)

para_prior = list( beta_occu = rep(1000,2 * ncol(envX))
                   ,beta_det = rep(10,2 * (ncol(envX)) )
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

n_dataset = 100

eta_intra_1 = matrix(NA,nrow = n_dataset,ncol = 5000)
eta_intra_2 = matrix(NA,nrow = n_dataset,ncol = 5000)
beta_1 = matrix(NA,nrow = n_dataset,ncol = 5000)
beta_2 = matrix(NA,nrow = n_dataset,ncol = 5000)
gamma_oc = matrix(NA,nrow = n_dataset,ncol = 5000)
gamma_de = matrix(NA,nrow = n_dataset,ncol = 5000)


eta_intra_1 = read.csv("./Main analysis/Results/Big_simulation/APIS/S/eta_intra_1.csv",row.names = 1)
eta_intra_2 = read.csv("./Main analysis/Results/Big_simulation/APIS/S/eta_intra_2.csv",row.names = 1)

beta_1 = read.csv("./Main analysis/Results/Big_simulation/APIS/S/beta_1.csv",row.names = 1)
beta_2 = read.csv("./Main analysis/Results/Big_simulation/APIS/S/beta_2.csv",row.names = 1)
  
gamma_oc = read.csv("./Main analysis/Results/Big_simulation/APIS/S/gamma_oc.csv",row.names = 1)
gamma_de = read.csv("./Main analysis/Results/Big_simulation/APIS/S/gamma_de.csv",row.names = 1)


###### Simulate Data ######
set.seed(6797)


for(i in 34:n_dataset){
  cat(i,"\n\n")
  MRF = getMRF(theta,envX,distM_full,link_map,distM_mainland,link_mainland = link_mainland * exp(-2*distM_mainland),
	  	 int_range_intra="nn",int_range_inter="nn")

  Z_simu = IsingSamplerCpp(1, MRF$A, MRF$thr, 1, 30, c(-1,1), T,NA+MRF$thr) ## take true occupancy


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
                                            , distM=distM_full,link_map=link_map
                                                , dist_mainland =  distM_mainland , link_mainland =  link_mainland * exp(-2*distM_mainland)
                                            , int_range_intra="nn",int_range_inter="nn"                                          
                                            , seed = NULL
                                            , ini = theta,thin.by = 10,report.by = 5000,nIter = 150,method = "MH",Gibbs = T)


  eta_intra_1[i,] = kk$theta.mcmc$eta_intra[,1]
  eta_intra_2[i,] = kk$theta.mcmc$eta_intra[,2]
  beta_1[i,] = kk$theta.mcmc$eta_inter[,1]
  beta_2[i,] = kk$theta.mcmc$eta_inter[,2]
  gamma_oc[i,] = kk$theta.mcmc$spp_mat[,2]
  gamma_de[i,] = kk$theta.mcmc$spp_mat_det[,2]

  write.csv(eta_intra_1,"./Main analysis/Results/Big_simulation/APIS/S/eta_intra_1.csv")
  write.csv(eta_intra_2,"./Main analysis/Results/Big_simulation/APIS/S/eta_intra_2.csv")
  
  write.csv(beta_1,"./Main analysis/Results/Big_simulation/APIS/S/beta_1.csv")
  write.csv(beta_2,"./Main analysis/Results/Big_simulation/APIS/S/beta_2.csv")
  
  write.csv(gamma_oc,"./Main analysis/Results/Big_simulation/APIS/S/gamma_oc.csv")
  write.csv(gamma_de,"./Main analysis/Results/Big_simulation/APIS/S/gamma_de.csv")
  boxplot(t(as.matrix(gamma_oc)))
  abline(0.25,0)
}
            




















