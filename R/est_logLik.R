est_logLik_occu = function(theta_0, Y, theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp",Z_vec){
	# Descombes 1999 paper, calculate Z_theta/Z_theta0, using samples Y from theta_0, I will use posterior mean.
	n_theta = length(theta_0)
	theta_diff = lapply(1:n_theta,function(i,theta_0,theta){theta_0[[i]]-theta[[i]]},theta_0,theta)
	partition_ratio = mean( 
		exp(-Hamiltonian(theta_diff,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Y))
	)# ratio of partitioning function: Z_theta/Z_theta0
	
	# this is the estimation of a single theta, we need to loop through every sample
	-Hamiltonian(theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_vec)-log(partition_ratio)
}

make_list_version_mcmc = function(mcmc_list,theta_0){ # gonna return a list, with single element is a theta list
	n_mcmc = nrow(mcmc_list[[1]])
	nspp = sqrt( ncol(mcmc_list[[length(mcmc_list)-1]]))
	nparas = length(mcmc_list)
	lapply(1:n_mcmc,function(i,mcmc_list,theta_0,nspp,nparas){
	    temp = lapply(1:(length(theta_0)-2),function(k,mcmc_list,i){
		    mcmc_list[[k]][i,]
		},mcmc_list,i)
		names(temp) = names(theta_0)[1:(length(theta_0)-2)]
		
		temp$spp_mat = matrix(mcmc_list[[nparas-1]][i,],nspp,nspp) # sppmat, need to be formatted
		temp$spp_mat_det = matrix(mcmc_list[[nparas]][i,],nspp,nspp)
	
	},mcmc_list,theta_0,nspp,nparas)
}



