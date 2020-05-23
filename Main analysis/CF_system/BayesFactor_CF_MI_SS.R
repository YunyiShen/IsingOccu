source("./R/est_logLik.R")
source("./R/misc.R")
Rcpp::sourceCpp('./src/IsingCpp_CFTP_sparse.cpp')
require(Matrix)
require(jsonlite)
CF_MI = read_json("./Main analysis/Results/imperfect_obs/Gibbs_norm/CF/CF_MI_60d_0.1.json",T)
CF_SS = read_json("./Main analysis/Results/imperfect_obs/Gibbs_norm/CF/CF_SS_60d_0.1.json",T)

CF_MI$means$spp_mat = matrix(CF_MI$means$spp_mat,2,2)
CF_MI$means$spp_mat_det = matrix(CF_MI$means$spp_mat_det,2,2)

CF_SS$means$spp_mat = matrix(CF_SS$means$spp_mat,2,2)
CF_SS$means$spp_mat_det = matrix(CF_SS$means$spp_mat_det,2,2)


CF_MI$linkmap <- lapply(CF_MI$linkmap,function(w){Matrix(w,sparse = T)})
CF_SS$linkmap <- lapply(CF_SS$linkmap,function(w){Matrix(w,sparse = T)})



logBF_MIminusSS = logBF(theta_a_mcmc = make_list_version_mcmc( CF_MI$theta.mcmc,theta)
                               ,envX_a=CF_MI$envX
                               ,distM = CF_MI$distM
                               ,link_map_a = CF_MI$linkmap 
                               ,dist_mainland = CF_MI$dist_mainland
                               ,link_mainland_a = exp(-2*CF_MI$dist_mainland)
                               ,int_range_intra_a="nn"
                               ,int_range_inter_a="nn"
                               ,Z_a_mcmc = CF_MI$Z.mcmc
                               , detX_a = NULL
                               , theta_a_point = CF_MI$means
                               ,theta_mcmc = make_list_version_mcmc( CF_SS$theta.mcmc,theta)
                               ,envX=CF_SS$envX
                               ,link_map = CF_SS$linkmap
                               ,link_mainland = CF_SS$link_mainland * exp(-2* CF_SS$dist_mainland)
                               ,int_range_intra="nn"
                               ,int_range_inter="nn"
                               ,Z_mcmc = CF_SS$Z.mcmc
                               , detX = NULL
                               , theta_point = CF_SS$means
                               , detmat=detmat, nrep=1
                               , nY = 15000,nIter = 150,method = "MH")
