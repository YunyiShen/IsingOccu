source("./R/misc.R")
source("./R/Main_Sampler.R")
source("./R/Simu_data_Sampling.R")
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


###### True Parameter Setting ######

spp_mat = matrix(1,2,2)
diag(spp_mat)=0
spp_mat = as(spp_mat,'dsCMatrix')
envX = matrix(1,155,1)

theta = list(beta_occu = c(-.5,.5),
             beta_det = c(-.5,-.5),
             eta_intra = c(0.15,0.15),
             eta_inter = c(1,-1),
             spp_mat = 0 * spp_mat,
             spp_mat_det = -.2 * spp_mat)

link_map = 
  list(inter = link_outer * exp(-2*distM_full),
       intra = link_inner)

nrep = 1
nspp = 2
nperiod = 8
nsite = 155



###### Simulate Data ######
set.seed(42)
MRF = getMRF(theta,envX,distM_full,link_map,distM_mainland,link_mainland = link_mainland * exp(-2*distM_mainland),
			 int_range_intra="nn",int_range_inter="nn")

Z_simu = IsingSamplerCpp(1, MRF$A, MRF$thr, 1, 30, c(-1,1), T,NA+MRF$thr) ## take true occupancy

detmat_format = lapply(1:nrep,function(dummy,nperiod,nsite,nspp){
	matrix(-1,nrow = nsite*nspp,ncol = nperiod)
},nperiod,nsite,nspp) # just for formating, hopefully it works

detmat_simu = Sample_Ising_detection_rep(nrep,nperiod,envX,NULL,
										 theta$beta_det,theta$spp_mat_det,Z_simu,
										 detmat_format,nIter=100,n=1, method = "CFTP")

#year1_dethis = detmat_simu[[1]]
Z_absolute = (sapply(detmat_simu,function(detmat_i){rowSums((detmat_i+1)/2)>0})) * 2 - 1


require(ggplot2)
require(ggmap)

Posi = read.csv("./data/APIS/CT_posi_only_island.csv")

map_data = data.frame(Posi,Z1 = Z_simu[1:155],Z2 = Z_simu[156:310])

APIS_map = get_stamenmap(bbox = c(
  left = -91.05, bottom = 46.75, 
  right =-90.35, top = 47.10),zoom = 12)

Z1_map = ggmap(APIS_map) + 
  geom_point(aes(x=Long,y=Lat,color= Z1),data = map_data,size = 1.2)+
  theme(text = element_text(size=15))


Z2_map = ggmap(APIS_map) + 
  geom_point(aes(x=Long,y=Lat,color= Z2),data = map_data,size = 1.2)+
  theme(text = element_text(size=15))




###### Run the Model! ######

vars_prop = list( beta_occu = c(5e-3,5e-3)
                  ,beta_det = rep(1e-3,nspp * ( ncol(envX)) ) # no extra det thing
                  ,eta_intra = rep(1e-3,nspp)
                  ,eta_inter = c(5e-3,5e-3)
                  ,spp_mat = 1e-2
                  ,spp_mat_det = 1e-2)

para_prior = list( beta_occu = rep(1000,2 * ncol(envX))
                   ,beta_det = rep(1000,2 * (ncol(envX)) )
                   ,eta_intra = rep(1000,nspp)
                   ,eta_inter = rep(1000,nspp)
                   ,d_intra=rep(1000,nspp)
                   ,d_inter = rep(1000,nspp)
                   ,spp_mat = 1000
                   ,spp_mat_det = 1000)


kk = IsingOccu.fit.Murray.sampler_Ising_det(X = envX, detmat =  detmat_simu
                                            , detX =  NULL
                                            , mcmc.iter = 1e6, burn.in = 50000
                                            , vars_prop = vars_prop
                                            , para_prior = para_prior
                                            , Zprop_rate = 1
                                            , uni_prior = F
                                            , distM=distM_full,link_map=link_map
                                            , dist_mainland =  distM_mainland , link_mainland =  link_mainland * exp(-2*distM_mainland)
                                            , int_range_intra="nn",int_range_inter="nn"                                          
                                            , seed = 42
                                            , ini = theta,thin.by = 300,report.by = 5000,nIter = 150,method = "CFTP")


      save.image("Test_MH_neutral_1e6_with_intra.RData")






