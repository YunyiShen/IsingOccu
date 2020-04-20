source("./R/est_logLik.R")
source("./R/misc.R")
Rcpp::sourceCpp('./src/IsingCpp_CFTP_sparse.cpp')
require(Matrix)
require(jsonlite)
FM_MI = read_json("./Main analysis/Results/imperfect_obs/Gibbs_norm/CF/CF_MI_60d_0.1.json",T)
FM_SS = read_json("./Main analysis/Results/imperfect_obs/Gibbs_norm/CF/CF_SS_60d_0.1.json",T)

FM_MI$means$spp_mat = matrix(FM_MI$means$spp_mat,2,2)
FM_MI$means$spp_mat_det = matrix(FM_MI$means$spp_mat_det,2,2)

FM_SS$means$spp_mat = matrix(FM_SS$means$spp_mat,2,2)
FM_SS$means$spp_mat_det = matrix(FM_SS$means$spp_mat_det,2,2)


FM_MI$linkmap <- lapply(FM_MI$linkmap,function(w){Matrix(w,sparse = T)})
FM_SS$linkmap <- lapply(FM_SS$linkmap,function(w){Matrix(w,sparse = T)})



logBF_MIminusSS = logBF(theta_a_mcmc = make_list_version_mcmc( FM_MI$theta.mcmc,theta)
                               ,envX_a=FM_MI$envX
                               ,distM = FM_MI$distM
                               ,link_map_a = FM_MI$linkmap 
                               ,dist_mainland = FM_MI$dist_mainland
                               ,link_mainland_a = exp(-2*FM_MI$dist_mainland)
                               ,int_range_intra_a="nn"
                               ,int_range_inter_a="nn"
                               ,Z_a_mcmc = FM_MI$Z.mcmc
                               , detX_a = NULL
                               , theta_a_point = FM_MI$means
                               ,theta_mcmc = make_list_version_mcmc( FM_SS$theta.mcmc,theta)
                               ,envX=FM_SS$envX
                               ,link_map = FM_SS$linkmap
                               ,link_mainland = FM_SS$link_mainland * exp(-2* FM_SS$dist_mainland)
                               ,int_range_intra="nn"
                               ,int_range_inter="nn"
                               ,Z_mcmc = FM_SS$Z.mcmc
                               , detX = NULL
                               , theta_point = FM_SS$means
                               , detmat=detmat, nrep=1
                               , nY = 15000,nIter = 50,method = "CFTP")
