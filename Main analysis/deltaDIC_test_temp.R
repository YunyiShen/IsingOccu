source("./R/est_logLik.R")
source("./R/misc_island.R")
Rcpp::sourceCpp('./src/IsingCpp_CFTP_sparse.cpp')
require(Matrix)
list_version_mcmc = make_list_version_mcmc(kk$theta.mcmc,theta)

pointest = kk$means
pointest$spp_mat = matrix(pointest$spp_mat,nrow = 2)
pointest$spp_mat_det = matrix(pointest$spp_mat_det,nrow = 2)

#list_version_mcmc_sample = lapply(1:100,function(i,a){a[[i]]},a = list_version_mcmc)

deltaDIC_try = deltaDIC(theta_a_mcmc = list_version_mcmc
                        ,envX_a = envX,distM = distM_full,link_map_a = link_map
                        ,dist_mainland = distM_mainland,link_mainland_a = link_mainland * exp(-2*distM_mainland)
                        ,int_range_intra_a="nn"
                        ,int_range_inter_a="nn"
                        ,kk$Z.mcmc, detX_a = NULL, theta_a_point = pointest
                                   
                        
                        ,theta_mcmc = list_version_mcmc,envX = envX,link_map,link_mainland = link_mainland * exp(-2*distM_mainland)
                        ,int_range_intra="nn",int_range_inter="nn"
                        ,Z_mcmc = kk$Z.mcmc, detX = NULL, theta_point = pointest, detmat = detmat, nrep = 1
                        , nY = 3000,nIter = 50)
