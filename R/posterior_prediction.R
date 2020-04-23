Posterior_predicting <- function(theta_mcmc, distM,dist_mainland,
                                 envX,link_map,
                                 link_mainland,
                                 int_range_intra="nn",int_range_inter="exp",
                                 detX, nperiod = 10 ,nIter = 100,method = "CFTP"){
  nsite <- nrow(envX)
  nspp <- nrow(theta_mcmc[[1]]$spp_mat)

  lapply(1:length(theta_mcmc),
                     function(i,theta,envX,detX
                              ,distM,link_map,dist_mainland
                              ,link_mainland,int_range_intra
                              ,int_range_inter,nperiod,nspp,nIter){
                       
                       MRF = getMRF(theta[[i]],envX,distM,link_map,dist_mainland,link_mainland,
                                    int_range_intra,int_range_inter)
                       
                       Z_simu = IsingSamplerCpp(1, MRF$A, MRF$thr, 1, nIter, c(-1,1), T,NA+MRF$thr) ## take true occupancy
                       if(is.na(sum(Z_simu))) Z_simu = IsingSamplerCpp(1, MRF$A, MRF$thr, 1, nIter+100, c(-1,1), F,NA+MRF$thr)
                       detmat_format = lapply(1,function(dummy,nperiod,nsite,nspp){
                         matrix(-1,nrow = nsite*nspp,ncol = nperiod)
                       },nperiod,nsite,nspp) # just for formating, hopefully it works
                       
                       detmat_simu = Sample_Ising_detection_rep(1,nperiod,envX,detX,
                                                                theta[[i]]$beta_det,theta[[i]]$spp_mat_det,Z_simu,
                                                                detmat_format,nIter=nIter,n=1, method = "CFTP")
                       
                       return(detmat_simu[[1]])
                       
                        },
                     theta_mcmc,envX,detX,distM,link_map,dist_mainland
                     ,link_mainland,int_range_intra
                     ,int_range_inter,nperiod,nspp,nIter)
  
}

