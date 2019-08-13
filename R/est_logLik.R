est_logLik_occu = function(theta_0,envX_0,link_map_0,link_mainland_0,int_range_intra_0="nn",int_range_inter_0="exp",Y # all settings for theta_0 model
						   ,theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp",Z_vec){
	# Descombes 1999 paper, calculate Z_theta/Z_theta0, using samples Y from theta_0, I will use posterior mean.
	nsite = nrow(envX)
	beta_occu = theta$beta_occu
	eta_intra = theta$eta_intra # intra spp, intra island if apply
	d_intra = theta$d_intra
	spp_mat = theta$spp_mat
	nspp = nrow(spp_mat)
	nrep = ncol(Z_vec)
	A_in = getintralayerGraph(distM,link_map$intra,eta_intra,d_intra,int_range = int_range_intra,spp_mat)
	eta_inter = theta$eta_inter # assume there is one 
	d_inter = theta$d_inter
	A_ex = getintralayerGraph(distM,link_map$inter,eta_inter,d_inter,int_range = int_range_inter,spp_mat) # graph among islands, if apply, distM should only contain graph 
	A=getfullGraph(A_ex,A_in,spp_mat)
	rm(A_ex,A_in)
	thr = lapply(1:nspp,
	  function(i,envX,beta_occu,dist_mainland,link_mainland,eta_inter,d_inter,int_range_inter){
	    envX %*% beta_occu[1:ncol(envX)+(i-1)*ncol(envX)] + 
		      mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
	    },envX,beta_occu,dist_mainland,link_mainland,eta_inter,d_inter,int_range_inter)
	thr = Reduce(rbind,thr)
	
	beta_occu_0 = theta_0$beta_occu
	eta_intra_0 = theta_0$eta_intra # intra spp, intra island if apply
	d_intra_0 = theta_0$d_intra
	spp_mat_0 = theta_0$spp_mat
	
	A_in_0 = getintralayerGraph(distM,link_map_0$intra,eta_intra_0,d_intra_0,int_range = int_range_intra_0,spp_mat_0)
	eta_inter_0 = theta_0$eta_inter # assume there is a 
	d_inter_0 = theta_0$d_inter
	A_ex_0 = getintralayerGraph(distM,link_map_0$inter,eta_inter_0,d_inter_0,int_range = int_range_inter_0,spp_mat_0) # graph among islands, if apply, distM should only contain graph 
	A_0=getfullGraph(A_ex_0,A_in_0,spp_mat_0)
	
	rm(A_in_0,A_ex_0)
	
	thr_0 = lapply(1:nspp,
	  function(i,envX,beta_occu,dist_mainland,link_mainland,eta_inter,d_inter,int_range_inter){
	    envX %*% beta_occu[1:ncol(envX)+(i-1)*ncol(envX)] + 
		      mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
	    },envX_0,beta_occu_0,dist_mainland,link_mainland_0,eta_inter_0,d_inter_0,int_range_inter_0)
	thr_0 = Reduce(rbind,thr_0)
	
	A_dif = as(A_0-A,'dsCMatrix')
	rm(A_0,A)
	
	thr_diff = thr_0-thr
	rm(thr,thr_0)
	
	
	partition_ratio = mean(
		apply(Y,2,function(X,graph,s){exp(-H(graph,X,s))},graph = A_dif,s = thr_dif)
	)# ratio of partitioning function: Z_theta/Z_theta_0
	
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
		
		temp$spp_mat = Matrix(mcmc_list[[nparas-1]][i,],nspp,nspp,sparse = T) # sppmat, need to be formatted
		temp$spp_mat_det = Matrix(mcmc_list[[nparas]][i,],nspp,nspp,sparse = T)
	  return(temp)
	},mcmc_list,theta_0,nspp,nparas)
}
# passed Aug 12 2019


