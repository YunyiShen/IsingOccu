IsingOccu.fit.Murray.sampler_Ising_det = function(X,detmat,detX
										,mcmc.iter = 10000, burn.in = 10 
										, vars_prop = list( beta_occu = rep(1e-5,2 * ncol(X))
										                    ,beta_det = rep(1e-5,2 * (ncol(detX[[1]][[1]]) + ncol(X)) )
										                    ,eta_intra = rep(1e-5,nspp)
										                    ,eta_inter = rep(1e-5,nspp*(nspp-1)/2)
										                    ,d_intra=rep(1e-5,nspp)
										                    ,d_inter = rep(1e-5,nspp)
										                    ,spp_mat = 1e-5
															,spp_mat_det = 1e-5)
										,vars_prior = 2000
										,Zprop_rate = .1
										#,Zprop_rate_missing_obs = .5
										,distM,link_map
										,dist_mainland , link_mainland
										,int_range_intra="nn",int_range_inter="exp"
										#,Z
										,seed = 42,ini,thin.by = 100,report.by=100,nIter=100){ # ini has same formate of theta
	
	cat("Initializing...\n\n")
	require(coda)
	require(IsingSampler)
	source("misc_island.R")
	set.seed(seed)
	
	#missing_obs = is.na(sum(Z))
	
	cat("Setting for imperfect observation and missing sites:\n")
	if(Zprop_rate==0) cat("    Perfect observation, given by Z\n")
	else cat("    Imperfect observation, Z propose with rate",Zprop_rate,"\n")
	#if(!missing_obs) cat("    No missing observation site\n\n")
	#else cat("    Missing some sites, Z propose with rate",Zprop_rate_missing_obs,"\n\n")
	
	cat("MCMC reported every",report.by,"iterations and thinned by",thin.by,"\n\n")
	
	nsite = (nrow(distM))
	
	theta = ini
	spp_neig = 1 *( spp_mat!=0 )
	
	nspp = nrow(spp_mat)
	nrep = length(detmat)
	
    #theta_tuta=ini
	ini = lapply(ini,as.matrix)
	theta.mcmc = list()
	for(i in 1:length(ini)){
	  theta.mcmc[[i]] = mcmc(matrix(nrow = floor(mcmc.iter/thin.by),ncol = length(ini[[i]])),thin = thin.by)
	  
	}
	
	names(theta.mcmc) = names(ini)
	
	
	Z.mcmc = mcmc(matrix(nrow = floor(mcmc.iter/thin.by),ncol = nsite*nspp*nrep),thin = thin.by)
	#Z_absolute = (sapply(detmat,rowSums)>0) * 2 - 1
	detmat_nona = lapply(detmat,function(mat){
	  mat[is.na(mat)]=-1
	  return(mat)
	})
	Z_absolute = (sapply(detmat_nona,function(detmat_i){rowSums((detmat_i+1)/2)>0})) * 2 - 1
	rm(detmat_nona)
	
	Z_curr = Z_absolute
	Z_temp = Z_absolute
	
	theta_curr = ini
	
	cat("Burn in...\n")
	
	accept_Z = 0
	accept_Z_missing_obs = 0
	accept_theta_occu = 0
	accept_theta_det = 0
	low_acc_Z = 0
	low_acc_Z_missing_obs = 0
	low_acc_theta_occu = 0
	low_acc_theta_det = 0
	propose_Z = 0
	propose_Z_missing_obs = 0
	timing = proc.time()
	n_para_group = length(theta_curr)
	for(i in 1:burn.in){# to burn in 
		#propose theta 
	  theta_prop = theta_curr
	  for(j in c(1:length(theta_curr))[-c(2,n_para_group)]){ # no detection proposing
	    theta_prop[[j]] = matrix( rnorm(length(theta_curr[[j]]),mean = 0,sd = sqrt(vars_prop[[j]])),nrow(theta_curr[[j]]),ncol(theta_curr[[j]]) )+ theta_curr[[j]]
	  }
	  
	  theta_prop$spp_mat=theta_prop$spp_mat * spp_neig	#	theta_prop$spp_mat=theta_prop$spp_mat * spp_neig
	  theta_prop$spp_mat = .5*(theta_prop$spp_mat + t( theta_prop$spp_mat)) # must be sym
	  Z_temp = rIsingOccu_multi(theta_prop,X,distM,link_map,dist_mainland , link_mainland,int_range_intra,int_range_inter,n=nrep,method = "CFTP",nIter = nIter)
	  # MH ratio
	  
	  Murray_ratio=Murray.ratio.Ising_det(theta_curr ,theta_prop
	                            ,Z_curr ,Z_curr
	                            ,Z_temp
	                            ,detmat
	                            ,vars_prior
	                             
	                            ,X, detX
	                            ,distM,link_map
	                            ,dist_mainland , link_mainland
	                            ,int_range_intra,int_range_inter)
	  r = runif(1)
	  if(is.na(Murray_ratio)){
	    next
	  }
	  if(Murray_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
	  if(r<=Murray_ratio){
	    theta_curr=theta_prop
	    #Z_temp = Z_temp
	    accept_theta_occu = accept_theta_occu + 1
	  }
	    
		
		theta_prop = theta_curr
		theta_prop[[2]] = matrix( rnorm(length(theta_curr[[2]]),mean = 0,sd = sqrt(vars_prop[[2]])),nrow(theta_curr[[2]]),ncol(theta_curr[[2]]) )+ theta_curr[[2]]
		theta_prop[[n_para_group]] = matrix( rnorm(length(theta_curr[[n_para_group]]),mean = 0,sd = sqrt(vars_prop[[n_para_group]])),nrow(theta_curr[[n_para_group]]),ncol(theta_curr[[n_para_group]]) )+ theta_curr[[n_para_group]]
		theta_prop$spp_mat_det=theta_prop$spp_mat_det * spp_neig
		theta_prop$spp_mat_det = .5*(theta_prop$spp_mat_det + t( theta_prop$spp_mat_det)) # must be sym
		
		Murray_ratio=Murray.ratio.Ising_det(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp
						,detmat
						,vars_prior
						 
						,X, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(is.na(Murray_ratio)) Murray_ratio=0
		if(Murray_ratio<exp(-10)) low_acc_theta_det = low_acc_theta_det + 1
		if(r<=Murray_ratio){
			theta_curr=theta_prop
			accept_theta_det = accept_theta_det + 1
		}
		
		Z_prop = Z_curr
		
		
		if(runif(1)<Zprop_rate) {
		  absence = which(Z_absolute==-1)
		  #abs_obs = absence[!absence%in%no_obs]
		  flip = sample(absence,1)
			Z_prop[flip]=-Z_prop[flip]
			propose_Z = propose_Z + 1
			
		}
		
		
		Murray_ratio=Murray.ratio.Ising_det(theta_curr ,theta_curr
						,Z_curr ,Z_prop
						,Z_temp
						,detmat
						,vars_prior
						 
						,X, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(is.na(Murray_ratio)) Murray_ratio = 0
		if(Murray_ratio<exp(-10)) low_acc_Z = low_acc_Z + 1
		if(r<=Murray_ratio){
			Z_curr = Z_prop
			accept_Z = accept_Z + 1
		}
		
		# # for missing site
		# Z_prop = Z_curr
		# if(missing_obs ){
		# 	if(runif(1)<Zprop_rate_missing_obs){
		# 		
		# 		propose_Z_missing_obs = propose_Z_missing_obs + 1
		# 		flip = sample((no_obs),1)
		# 		Z_prop[flip]=-Z_prop[flip]
		# 	}
		# 	Murray_ratio=Murray.ratio.Ising_det(theta_curr ,theta_curr
		# 				,Z_curr , Z_prop
		# 				,Z_temp
		# 				,detmat
		# 				,vars_prior
		# 				,X, detX
		# 				,distM,link_map
		# 				,dist_mainland , link_mainland
		# 				,int_range_intra,int_range_inter)
		# 	r = runif(1)
		# 	if(Murray_ratio<exp(-10)) low_acc_Z_missing_obs = low_acc_Z_missing_obs + 1
		# 	if(r<=Murray_ratio){
		# 		Z_curr = Z_prop
		# 		accept_Z_missing_obs = accept_Z_missing_obs + 1
		# 	}
		# }
		# 
		# 
		
		if(i%%report.by == 0) {
		  
		  cat("Burn in iteration",i-report.by+1,"to",i,":\n\n")
		  cat("    # of Z proposed for imperfect detection: ",propose_Z,"\n")
		  cat("    # of Z acceptance for imperfect detection: " , accept_Z-(report.by-propose_Z),"\n")
		  cat("    # of Z acceptance ratio <exp(-10): ",low_acc_Z,"\n\n")
		  
		#   #cat("# of Z proposed for missing sites: ",propose_Z,"\n")
		#   if(missing_obs){
		#     cat("    # of Z proposed for missing sites: ",propose_Z_missing_obs,"\n")
		# 	cat("    # of Z acceptance for missing sites: " , accept_Z_missing_obs-(report.by-propose_Z_missing_obs),"\n")
		# 	cat("    # of Z acceptance ratio <exp(-10): ",low_acc_Z_missing_obs,"\n\n")
		#   }
		  
		  cat("    # of occupancy theta acceptance: " , accept_theta_occu,"\n")
		  cat("    # of occupancy acceptance ratio <exp(-10): ",low_acc_theta_occu,"\n\n")
		  cat("    # of detection theta acceptance:" , accept_theta_det,"\n")
		  cat("    # of detection acceptance ratio <exp(-10): ",low_acc_theta_det,"\n\n")
		  timing = proc.time()- timing
		  cat("Time used in this" ,report.by,":",timing[1],"s\n")
		  cat("\n\n")
		  accept_Z = 0
		  accept_Z_missing_obs = 0
		  accept_theta_occu = 0
		  accept_theta_det = 0
		  low_acc_Z = 0
		  low_acc_Z_missing_obs = 0
		  low_acc_theta_occu = 0
		  low_acc_theta_det = 0
		  propose_Z = 0
		  propose_Z_missing_obs = 0
		  timing = proc.time()
		  }
	}
	cat("Start sampling...\n")
	accept_Z = 0
	accept_Z_missing_obs = 0
	accept_theta_occu = 0
	accept_theta_det = 0
	low_acc_Z = 0
	low_acc_Z_missing_obs = 0
	low_acc_theta_occu = 0
	low_acc_theta_det = 0
	propose_Z = 0
	propose_Z_missing_obs = 0
	timing = proc.time()	
	for(i in 1:(mcmc.iter)){ # for to save 
		#propose theta 
	  theta_prop = theta_curr
	  for(j in c(1:length(theta_curr))[-c(2,n_para_group)]){ # no detection proposing
	    theta_prop[[j]] = matrix( rnorm(length(theta_curr[[j]]),mean = 0,sd = sqrt(vars_prop[[j]])),nrow(theta_curr[[j]]),ncol(theta_curr[[j]]) )+ theta_curr[[j]]
	  }
	  
	  theta_prop$spp_mat=theta_prop$spp_mat * spp_neig	#	theta_prop$spp_mat=theta_prop$spp_mat * spp_neig
	  theta_prop$spp_mat = .5*(theta_prop$spp_mat + t( theta_prop$spp_mat)) # must be sym
		Z_temp = rIsingOccu_multi(theta_prop,X,distM,link_map,dist_mainland , link_mainland,int_range_intra,int_range_inter,n=nrep,method = "CFTP",nIter = nIter)

		
		# MH ratio
		Murray_ratio=Murray.ratio.Ising_det(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp
						,detmat
						,vars_prior
						 
						,X, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(is.na(Murray_ratio)){
		  if(i %% thin.by==0){
		    for(j in 1:length(theta_curr)){
		      theta.mcmc[[j]][i/thin.by,] =as.vector( theta_curr[[j]])
		    } # saving the results
		  }
		  next
		}
		if(Murray_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
		if(r<=Murray_ratio){
			theta_curr=theta_prop
			#Z_temp = Z_temp
			accept_theta_occu = accept_theta_occu + 1
		}
		
		
		theta_prop = theta_curr
		theta_prop[[2]] = matrix( rnorm(length(theta_curr[[2]]),mean = 0,sd = sqrt(vars_prop[[2]])),nrow(theta_curr[[2]]),ncol(theta_curr[[2]]) )+ theta_curr[[2]]
		theta_prop[[n_para_group]] = matrix( rnorm(length(theta_curr[[n_para_group]]),mean = 0,sd = sqrt(vars_prop[[n_para_group]])),nrow(theta_curr[[n_para_group]]),ncol(theta_curr[[n_para_group]]) )+ theta_curr[[n_para_group]]
		theta_prop$spp_mat_det=theta_prop$spp_mat_det * spp_neig
		theta_prop$spp_mat_det = .5*(theta_prop$spp_mat_det + t( theta_prop$spp_mat_det)) # must be sym

		
		Murray_ratio=Murray.ratio.Ising_det(theta_curr ,theta_prop
		                          ,Z_curr ,Z_curr
		                          ,Z_temp
		                          ,detmat
		                          ,vars_prior
		                           
		                          ,X, detX
		                          ,distM,link_map
		                          ,dist_mainland , link_mainland
		                          ,int_range_intra,int_range_inter)
		r = runif(1)
		if(is.na(Murray_ratio)) Murray_ratio = 0
		if(Murray_ratio<exp(-10)) low_acc_theta_det = low_acc_theta_det + 1
		if(r<=Murray_ratio){
		  theta_curr=theta_prop
		  accept_theta_det = accept_theta_det + 1
		}
		
		if(i %% thin.by==0){
	  	for(j in 1:length(theta_curr)){
		  	theta.mcmc[[j]][i/thin.by,] =as.vector( theta_curr[[j]])
	  	} # saving the results
		}
		
		
		
		Z_prop = Z_curr
		if(runif(1)<Zprop_rate) {
		  absence = which(Z_absolute==-1)
		  #abs_obs = absence[!absence%in%no_obs]
		  flip = sample(absence,1)
		  Z_prop[flip]=-Z_prop[flip]
		  propose_Z = propose_Z + 1
		  
		}
		
		
		Murray_ratio=Murray.ratio.Ising_det(theta_curr ,theta_curr
		                          ,Z_curr ,Z_prop
		                          ,Z_temp
		                          ,detmat
		                          ,vars_prior
		                           
		                          ,X, detX
		                          ,distM,link_map
		                          ,dist_mainland , link_mainland
		                          ,int_range_intra,int_range_inter)
		r = runif(1)
		if(is.na(Murray_ratio)) Murray_ratio = 0
		if(Murray_ratio<exp(-10)) low_acc_Z = low_acc_Z + 1
		if(r<=Murray_ratio){
		  #theta_curr=theta_prop
		  Z_curr = Z_prop
		  accept_Z = accept_Z + 1
		}
		
		
		# Z_prop = Z_curr
		# if(missing_obs ){
		# 	if(runif(1)<Zprop_rate_missing_obs){
		# 		
		# 		propose_Z_missing_obs = propose_Z_missing_obs + 1
		# 		flip = sample((no_obs),1)
		# 		Z_prop[flip]=-Z_prop[flip]
		# 	}
		# 	Murray_ratio=Murray.ratio.Ising_det(theta_curr ,theta_curr
		# 				,Z_curr ,Z_prop
		# 				,Z_temp
		# 				,detmat
		# 				,vars_prior
		# 				,X, detX
		# 				,distM,link_map
		# 				,dist_mainland , link_mainland
		# 				,int_range_intra,int_range_inter)
		# 	r = runif(1)
		# 	if(Murray_ratio<exp(-10)) low_acc_Z_missing_obs = low_acc_Z_missing_obs + 1
		# 	if(r<=Murray_ratio){
		# 		Z_curr = Z_prop
		# 		accept_Z_missing_obs = accept_Z_missing_obs + 1
		# 	}
		# }		
		
		if(i %% thin.by==0) Z.mcmc[i/thin.by,]=Z_curr
		if(i%%report.by == 0) { # reporting
		  cat("Sampling iteration",i-report.by+1,"to",i,":\n\n")
		  cat("    # of Z proposed: ",propose_Z,"\n")
		  cat("    # of Z acceptance: " , accept_Z-(report.by-propose_Z),"\n")
		  cat("    # of Z acceptance ratio <exp(-10): ",low_acc_Z,"\n\n")
		  
		#   if(missing_obs){
		#     cat("    # of Z proposed for missing sites: ",propose_Z_missing_obs,"\n")
		# 	cat("    # of Z acceptance for missing sites: " , accept_Z_missing_obs-(report.by-propose_Z_missing_obs),"\n")
		# 	cat("    # of Z acceptance ratio <exp(-10): ",low_acc_Z_missing_obs,"\n\n")
		#   }
		  
		  cat("    # of occupancy theta acceptance: " , accept_theta_occu,"\n")
		  cat("    # of occupancy acceptance ratio <exp(-10): ",low_acc_theta_occu,"\n\n")
		  #if(accept_theta_occu==0) cat(theta_curr[c(1:ncov,1:5+(ncov+ncov_det))],"\n\n")
		  cat("    # of detection theta acceptance: " , accept_theta_det,"\n")
		  cat("    # of detection acceptance ratio <exp(-10): ",low_acc_theta_det,"\n\n")
		  #if(accept_theta_det==0) cat(theta_curr[1:ncov_det + ncov],"\n\n")
		  timing = proc.time()-timing
		  cat("Time used in this" ,report.by,":",timing[1],"s\n")
		  cat("\n\n")
		  accept_Z = 0
		  accept_Z_missing_obs = 0
		  accept_theta_occu = 0
		  accept_theta_det = 0
		  low_acc_Z = 0
		  low_acc_Z_missing_obs = 0
		  low_acc_theta_occu = 0
		  low_acc_theta_det = 0
		  propose_Z = 0
		  propose_Z_missing_obs = 0
		  timing = proc.time()
		  }
	}
  
	theta.mean =lapply(theta.mcmc,function(thetaa){ apply(thetaa,2,mean,na.rm=T)})
  
	res = list(theta.mcmc = theta.mcmc,means = theta.mean,Z.mcmc = Z.mcmc,vars=vars_prop, interaction.range =list( int_range_inter,int_range_intra), envX=X)
	class(res)="IsingOccu_island_Ising_det"
	return(res)
}

