IsingOccu.fit.Moller.sampler = function(X,detmat,detX
										,mcmc.save = 10000, burn.in = 10 
										, vars_prior = 
											list(
												beta_occu = rep(10000,ncol(X))
												,beta_det = rep(1,ncol(detX))
												,eta_intra = rep(1,nspp)
												,eta_inter=rep(1,nspp*(nspp-1)/2
												,d_intra=rep(1,nspp)
												,d_inter = rep(1,nspp)
												,spp_mat = 1e-5
												)
												)
										,vars_prop = 2000
										,Zprop_rate = .1
										,distM,link_map
										,dist_mainland , link_mainland
										,int_range_intra="nn",int_range_inter="exp"
										,Z,seed = 42,ini){ # ini has same formate of theta
	require(coda)
	source("misc_island.R")
	cat("Initializing...\n\n")
	set.seed(seed)
	nsite = (nrow(distM))
	
	theta = ini
	
	beta_occu = theta$beta_occu # this will be a matrix for cols are species
	beta_det = theta$beta_det
	eta_intra = theta$eta_intra # intra spp, intra island if apply
	d_intra = theta$d_intra
	eta_inter = theta$eta_inter
	d_inter = theta$d_inter # inter island scaling factor
	spp_mat = theta$spp_mat
	
	spp_neig = 1 *( spp_mat>0 )
	
	nspp = nrow(spp_mat)
	nrep = ncol(Z_vec)
	
    theta_tuta=ini
	theta_tuta = lapply(theta_tuta,as.matrix)
	#names(theta_tuta) = ("beta_occu","beta_det","eta_intra","d","eta_inter")
	theta.mcmc = list( # for mcmc results
		beta_occu.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(beta_occu))),
		beta_det.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(beta_det))),
		eta_intra.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(eta_intra))),
		d_intra.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(d_intra))),
		eta_inter.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(eta_inter))),
		d_intra.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(d_intra))),
		spp_mat.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(spp_mat)))
	)
	
	
	Z.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = nrow(detmat)))
	Z_absolute = (as.numeric(rowSums(detmat)>0)) * 2 - 1
	Z_tuta = Z_absolute
	
	
	Z_curr = Z_tuta
	# x_curr = detmat
	Z_temp_curr = Z_tuta
	
	theta_curr = theta_tuta
	cat("Burn in...\n")
	accept_Z = 0
	accept_theta_occu = 0
	accept_theta_det = 0
	low_acc_Z = 0
	low_acc_theta_occu = 0
	low_acc_theta_det = 0
	propose_Z = 0
	timing = proc.time()
	for(i in 1:burn.in){# to burn in 
		#propose theta 
		theta_prop = list()
		for(j in c(1:length(theta_curr))[-2]){ # no detection proposing
			theta_prop[[j]] = matrix( rnorm(length(theta_curr[[j]]),mean = 0,sd = sqrt(vars_prop[[j]])),nrow(theta_curr[[j]]),ncol(theta_curr[[j]]) )+ theta_curr[[j]]
		}
		theta_prop$spp_mat=theta_prop$spp_mat * spp_neig
		#propose Z from uniform distribution 
		Z_temp_prop = rIsingOccu_multi = function(theta_prop,envX,distM,link_map,dist_mainland , link_mainland,int_range_intra,int_range_inter,n=nrep,method = "CFTP",nIter = 100)
		# propose x, from the likelihood
		# x_prop = IsingOccu_multispp_sample.detection(theta_prop, X, Z_temp_prop ,detmat, detX,nspp)
		# MH ratio
		Moller_ratio=Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp_curr, Z_temp_prop
						,detmat
						,vars_prior
						,theta_tuta
						,envX, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(Moller_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			#Z_curr = Z_prop
			# x_curr = x_prop
			Z_temp_curr = Z_temp_prop
			accept_theta_occu = accept_theta_occu + 1
		}
		
		theta_prop[[2]] = matrix( rnorm(length(theta_curr[[2]]),mean = 0,sd = sqrt(vars_prop[[2]])),nrow(theta_curr[[2]]),ncol(theta_curr[[2]]) )+ theta_curr[[2]]
		
		Moller_ratio=Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp_curr, Z_temp_curr
						,detmat
						,vars_prior
						,theta_tuta
						,envX, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(Moller_ratio<exp(-10)) low_acc_theta_det = low_acc_theta_det + 1
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			#Z_curr = Z_prop
			# x_curr = x_prop
			Z_temp_curr = Z_temp_prop
			accept_theta_det = accept_theta_det + 1
		}
		
		Z_prop = Z_curr
		flip = sample(which(Z_absolute==-1),1)
		
		if(runif(1)<Zprop_rate) {
			Z_prop[flip]=-Z_prop[flip]
			propose_Z = propose_Z + 1
		}
		
		
		Z_prop = (Z_absolute==1) + (Z_absolute==-1) * ((runif(length(Z_absolute))>=0.5) * 2 - 1)
		Moller_ratio=Moller.ratio(theta_curr ,theta_curr
						,Z_curr ,Z_prop
						,Z_temp_curr, Z_temp_prop
						,detmat
						,vars_prior
						,theta_tuta
						,envX, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(Moller_ratio<exp(-10)) low_acc_Z = low_acc_Z + 1
		if(r<=Moller_ratio){
			#theta_curr=theta_prop
			Z_curr = Z_prop
			accept_Z = accept_Z + 1
		}
		if(i%%100 == 0) {
		  
		  cat("Burn in iteration: #",i,"\n")
		  cat("# of Z proposed: ",propose_Z,"\n")
		  cat("# of Z acceptance: " , accept_Z-(100-propose_Z),"\n")
		  cat("# of Z acceptance ratio <exp(-10): ",low_acc_Z,"\n\n")
		  cat("# of occupancy theta acceptance: " , accept_theta_occu,"\n")
		  cat("# of occupancy acceptance ratio <exp(-10): ",low_acc_theta_occu,"\n\n")
		  #if(accept_theta_occu==0) cat(theta_curr[c(1:ncov,1:5+(ncov+ncov_det))],"\n\n")
		  cat("# of detection theta acceptance:" , accept_theta_det,"\n")
		  cat("# of detection acceptance ratio <exp(-10): ",low_acc_theta_det,"\n\n")
		  #if(accept_theta_det==0) cat(theta_curr[1:ncov_det + ncov],"\n\n")
		  timing = proc.time()- timing
		  cat("Time used in this 100:",timing[1],"s\n")
		  cat("\n\n")
		  timing = proc.time()
		  accept_Z = 0
		  low_acc_Z = 0
		  accept_theta_occu = 0
		  low_acc_theta_occu = 0
		  accept_theta_det = 0
		  low_acc_theta_det = 0
		  propose_Z = 0
		  }
	}
	cat("Start sampling...\n")
	accept_Z = 0
	low_acc_Z = 0
	accept_theta_occu = 0
	low_acc_theta_occu = 0
	accept_theta_det = 0
	low_acc_theta_det = 0
	timing = proc.time()
	
	for(i in 1:(mcmc.save)){ # for to save 
		#propose theta 
		#theta_prop = rnorm(length(theta_curr),mean = theta_curr,sd = sqrt(vars_prop))
		for(j in c(1:length(theta_prop))[-2]){
			theta_prop[[j]] = matrix( rnorm(length(theta_curr[[j]]),mean = 0,sd = sqrt(vars_prop[[j]])),nrow(theta_curr[[j]]),ncol(theta_curr[[j]]) )+ theta_curr[[j]]
		}
		theta_prop$spp_mat=theta_prop$spp_mat * spp_neig
		#propose Z from uniform distribution 
		Z_temp_prop = rIsingOccu_multi = function(theta_prop,envX,distM,link_map,dist_mainland , link_mainland,int_range_intra,int_range_inter,n=nrep,method = "CFTP",nIter = 100)
		# propose x, from the likelihood
		# x_prop = IsingOccu_sample.detection(theta_prop, X, Z_temp_prop ,detmat, detX)
		# MH ratio
		
		Moller_ratio=Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp_curr, Z_temp_prop
						,detmat
						,vars_prior
						,theta_tuta
						,envX, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(Moller_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			#Z_curr = Z_prop
			# x_curr = x_prop
			Z_temp_curr = Z_temp_prop
			accept_theta_occu = accept_theta_occu + 1
		}
		
		
		theta_prop[[2]] = matrix( rnorm(length(theta_curr[[2]]),mean = 0,sd = sqrt(vars_prop[[2]])),nrow(theta_curr[[2]]),ncol(theta_curr[[2]]) )+ theta_curr[[2]]
		
		Moller_ratio=Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp_curr, Z_temp_curr
						,detmat
						,vars_prior
						,theta_tuta
						,envX, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(Moller_ratio<exp(-10)) low_acc_theta_det = low_acc_theta_det + 1
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			#Z_curr = Z_prop
			# x_curr = x_prop
			Z_temp_curr = Z_temp_prop
			accept_theta_det = accept_theta_det + 1
		}
		
		
		
		for(j in 1:length(theta_curr)){
			theta_mcmc[[j]][i,] =as.vector( theta_curr[[j]])
		}
		
		
		
		Z_prop = (Z_absolute==1) + (Z_absolute==-1) * ((runif(length(Z_absolute))>=0.5) * 2 - 1)
		Moller_ratio=Moller.ratio(theta_curr ,theta_curr
						,Z_curr ,Z_prop
						,Z_temp_curr, Z_temp_prop
						,detmat
						,vars_prior
						,theta_tuta
						,envX, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(Moller_ratio<exp(-10)) low_acc_Z = low_acc_Z + 1
		if(r<=Moller_ratio){
			#theta_curr=theta_prop
			Z_curr = Z_prop
			# x_curr = x_prop
			#Z_temp_curr = Z_temp_prop
		}
		
		
		Z.mcmc[i,]=Z_curr
		if(i%%100 == 0) { # reporting
		  cat("Sampling iteration: #",i,"\n")
		  cat("# of Z proposed: ",propose_Z,"\n")
		  cat("# of Z acceptance: " , accept_Z-(100-propose_Z),"\n")
		  cat("# of Z acceptance ratio <exp(-10): ",low_acc_Z,"\n\n")
		  cat("# of occupancy theta acceptance: " , accept_theta_occu,"\n")
		  cat("# of occupancy acceptance ratio <exp(-10): ",low_acc_theta_occu,"\n\n")
		  #if(accept_theta_occu==0) cat(theta_curr[c(1:ncov,1:5+(ncov+ncov_det))],"\n\n")
		  cat("# of detection theta acceptance: " , accept_theta_det,"\n")
		  cat("# of detection acceptance ratio <exp(-10): ",low_acc_theta_det,"\n\n")
		  #if(accept_theta_det==0) cat(theta_curr[1:ncov_det + ncov],"\n\n")
		  timing = proc.time()-timing
		  cat("Time used in this 100:",timing[1],"s\n")
		  cat("\n\n")
		  timing = proc.time()
		  accept_Z = 0
		  low_acc_Z = 0
		  accept_theta_occu = 0
		  low_acc_theta_occu = 0
		  accept_theta_det = 0
		  low_acc_theta_det = 0
		  propose_Z=0
		  }
	}

	
	res = list(theta.mcmc = theta.mcmc,theta.mean = apply(theta.mcmc,1,mean),vars=vars, interaction.range = int_range, graph = graph, envX=X)
	class(res)="IsingOccu_multispp.Moller"
	return(res)
}

