source("./R/est_logLik.R")
source("./R/misc.R")
Rcpp::sourceCpp('./src/IsingCpp_CFTP_sparse.cpp')
require(Matrix)
FM_MI = read_json("FM_MI_50K.json",T)
FM_SS = read_json("FM_SS_50K.json",T)

FM_MI$means$spp_mat = matrix(FM_MI$means$spp_mat,2,2)
FM_MI$means$spp_mat_det = matrix(FM_MI$means$spp_mat_det,2,2)

FM_SS$means$spp_mat = matrix(FM_SS$means$spp_mat,2,2)
FM_SS$means$spp_mat_det = matrix(FM_SS$means$spp_mat_det,2,2)

deltaDIC_MIminusSS = deltaDIC(theta_a_mcmc = make_list_version_mcmc( FM_MI$theta.mcmc,theta)
                               ,envX_a=FM_MI$envX
                               ,distM = distM_full
                               ,link_map_a = list(inter = 0*link_outer,intra=link_inner)
                               ,dist_mainland = distM_mainland
                               ,link_mainland_a = exp(-2*distM_mainland)
                               ,int_range_intra_a="nn"
                               ,int_range_inter_a="nn"
                               ,Z_a_mcmc = FM_MI$Z.mcmc
                               , detX_a = NULL
                               , theta_a_point = FM_MI$means
                               ,theta_mcmc = make_list_version_mcmc( FM_SS$theta.mcmc,theta)
                               ,envX=FM_MI$envX
                               ,link_map = link_map
                               ,link_mainland = link_mainland * exp(-2*distM_mainland)
                               ,int_range_intra="nn"
                               ,int_range_inter="nn"
                               ,Z_mcmc = FM_SS$Z.mcmc
                               , detX = NULL
                               , theta_point = FM_SS$means
                               , detmat=detmat, nrep=1
                               , nY = 3000,nIter = 100,method = "CFTP")
