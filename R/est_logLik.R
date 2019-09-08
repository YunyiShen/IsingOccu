est_logLik_occu = function(A_0,thr_0,Y # all settings for theta_0 model
						   ,A,thr,Z){
	# Descombes 1999 paper, calculate Z_theta/Z_theta0, using samples Y from theta_0, I will use posterior mean.
	A_dif = as(A_0-A,'dsCMatrix')
	rm(A_0)
	thr_diff = thr_0-thr	
	Ham = apply(Y,1,function(X,graph,s){(H(graph,X,s))},graph = A_dif,s = thr_diff) # Hamiltonian with parameters as difference and sample of theta_0
	robust_cri = (max(Ham)-min(Ham))/max(Ham)
	partition_ratio = mean(
		exp(-Ham)
	)# ratio of partitioning function: Z_theta/Z_theta_0
	
	# this is the estimation of a single theta, we need to loop through every sample
	logLik_est = -H(A,Z,thr)-log(partition_ratio)
	return(list(logLik_est=logLik_est,robust_cri=robust_cri))
}

logLik_full = function(theta,envX,detX,distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp",Z, detmat, A_0,thr_0,Y){
    ## logLik for detection part
	nspp = nrow(theta$spp_mat)
    nperiod = ncol(detmat[[1]])
    beta_det = theta$beta_det
    logLdata = Pdet_Ising_rep(nrep,nperiod,envX,detX,theta$beta_det,theta$spp_mat_det,as.matrix(Z),detmat)
	
	## now estimate the logLik for occu
	
	nsite = nrow(envX)
	beta_occu = theta$beta_occu
	eta_intra = theta$eta_intra # intra spp, intra island if apply
	d_intra = theta$d_intra
	spp_mat = theta$spp_mat
	nspp = nrow(spp_mat)
	nrep = ncol(Z)
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
    
	
	logLoccu = est_logLik_occu(A_0,thr_0,Y,A,thr,Z)
	# sum them up
	return(list(logLik = logLoccu$logLik_est + logLdata,robust = logLoccu$robust_cri))

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

# we calculate the deltaDIC: DIC(theta_a)-DIC(theta), luckly delta DIC is summable
deltaDIC = function(theta_a_mcmc,envX_a,distM,link_map_a,dist_mainland,link_mainland_a,int_range_intra_a="nn",int_range_inter_a="exp",Z_a_mcmc, detX_a, theta_a_point
					,theta_mcmc,envX,link_map,link_mainland,int_range_intra="nn",int_range_inter="exp",Z_mcmc, detX, theta_point, detmat, nrep, nY = 3000,nIter = 100,method = "CFTP"){
    nsite = nrow(envX)
	beta_occu_point = theta_point$beta_occu
	eta_intra_point = theta_point$eta_intra # intra spp, intra island if apply
	d_intra_point = theta_point$d_intra
	spp_mat_point = theta_point$spp_mat
	nspp = nrow(spp_mat_point)
	A_in_point = getintralayerGraph(distM,link_map$intra,eta_intra_point,d_intra_point,int_range = int_range_intra,spp_mat_point)
	eta_inter_point = theta_point$eta_inter # assume there is one 
	d_inter_point = theta_point$d_inter
	A_ex_point = getintralayerGraph(distM,link_map$inter,eta_inter_point,d_inter_point,int_range = int_range_inter,spp_mat_point) # graph among islands, if apply, distM should only contain graph 
	A_point=getfullGraph(A_ex_point,A_in_point,spp_mat_point)
	rm(A_ex_point,A_in_point)
	thr_point = lapply(1:nspp,
	  function(i,envX,beta_occu,dist_mainland,link_mainland,eta_inter,d_inter,int_range_inter){
	    envX %*% beta_occu[1:ncol(envX)+(i-1)*ncol(envX)] + 
		      mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
	    },envX,beta_occu_point,dist_mainland,link_mainland,eta_inter_point,d_inter_point,int_range_inter)
	thr_point = Reduce(rbind,thr_point)

    beta_occu_0_point = theta_a_point$beta_occu
	eta_intra_0_point = theta_a_point$eta_intra # intra spp, intra island if apply
	d_intra_0_point = theta_a_point$d_intra
	spp_mat_0_point = theta_a_point$spp_mat
	
	A_in_0_point = getintralayerGraph(distM,link_map_a$intra,eta_intra_0_point,d_intra_0_point,int_range = int_range_intra_a,spp_mat_0_point)
	eta_inter_0_point = theta_a_point$eta_inter # assume there is a 
	d_inter_0_point = theta_a_point$d_inter
	A_ex_0_point = getintralayerGraph(distM,link_map_a$inter,eta_inter_0_point,d_inter_0_point,int_range = int_range_inter_a,spp_mat_0_point) # graph among islands, if apply, distM should only contain graph 
	A_0_point=getfullGraph(A_ex_0_point,A_in_0_point,spp_mat_0_point)
	
	rm(A_in_0_point,A_ex_0_point)
	
	thr_0_point = lapply(1:nspp,
	  function(i,envX,beta_occu,dist_mainland,link_mainland,eta_inter,d_inter,int_range_inter){
	    envX %*% beta_occu[1:ncol(envX)+(i-1)*ncol(envX)] + 
		      mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
	    },envX_a,beta_occu_0_point,dist_mainland,link_mainland_a,eta_inter_0_point,d_inter_0_point,int_range_inter_a)
	thr_0_point = Reduce(rbind,thr_0_point)

    A_mean = as(0.5*(A_point + A_0_point),'dsCMatrix') # theta in 1999 paper, as the constant, we estimate ratio of partition function of theta_a/theta_mean and theta/theta_mean
    thr_mean = 0.5*(thr_point + thr_0_point)
    
	Y = IsingSamplerCpp(n=nY,graph = A_mean, thresholds = thr_mean, beta=1, responses = c(-1L, 1L),nIter = nIter,exact = (method=="CFTP"),constrain = NA+thr_mean)
    
	logLik_a = lapply(1:length(theta_a_mcmc),
					  function(i,theta,envX,detX
							   ,distM,link_map,dist_mainland
							   ,link_mainland,int_range_intra
							   ,int_range_inter,Z, detmat, A_0,thr_0,Y){
	                               logLik_full(theta[[i]],envX,detX
											   ,distM,link_map,dist_mainland
											   ,link_mainland,int_range_intra,int_range_inter
											   ,matrix(Z[i,],ncol = nrep), detmat, A_0,thr_0,Y)
	},theta_a_mcmc,envX_a,detX_a,distM,link_map_a,dist_mainland
	  ,link_mainland_a,int_range_intra_a
	  ,int_range_inter_a,Z_a_mcmc, detmat, A_mean,thr_mean,Y)
    
    
	logLik_theta = lapply(1:length(theta_mcmc),
					  function(i,theta,envX,detX
							   ,distM,link_map,dist_mainland
							   ,link_mainland,int_range_intra
							   ,int_range_inter,Z, detmat, A_0,thr_0,Y){
	                               logLik_full(theta[[i]],envX,detX
											   ,distM,link_map,dist_mainland
											   ,link_mainland,int_range_intra,int_range_inter
											   ,matrix(Z[i,],ncol = nrep), detmat, A_0,thr_0,Y)
	},theta_mcmc,envX,detX,distM,link_map,dist_mainland
	  ,link_mainland,int_range_intra
	  ,int_range_inter,Z_mcmc, detmat, A_mean,thr_mean,Y)
	
	rm(Y,A_mean,thr_mean)
	
	pure_logLik_a = lapply(logLik_a,function(a){a[[1]]})
	robust_cri_a = lapply(logLik_a,function(a){a[[2]]})

	pure_logLik_theta = lapply(logLik_theta,function(a){a[[1]]})
	robust_cri_theta = lapply(logLik_theta,function(a){a[[2]]})	
	
	logLik_a = Reduce(rbind,pure_logLik_a)
	logLik_theta = Reduce(rbind,pure_logLik_theta)
	
	pD_Gelman04_a = .5 * var(logLik_a)
	pD_Gelman04_t = .5 * var(logLik_theta)
    
	deltaDIC = - mean(logLik_a) - 2 * (pD_Gelman04_a) +  mean(logLik_theta) + 2 * (pD_Gelman04_t)
	
	robust_cri = list(robust_a = Reduce(mean,robust_cri_a),robust_theta = Reduce(mean,robust_cri_theta))
	
	return(list(deltaDIC = deltaDIC,robust_criterion = robust_cri))
}

