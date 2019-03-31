IsingOccu.fit.Murray.sampler = function(X,detmat,detX
										,mcmc.iter = 10000, burn.in = 10 
										, vars_prop = list( beta_occu = rep(1e-5,2 * ncol(X))
										                    ,beta_det = rep(1e-5,2 * (ncol(detX[[1]][[1]]) + ncol(X)) )
										                    ,eta_intra = rep(1e-5,nspp)
										                    ,eta_inter = rep(1e-5,nspp*(nspp-1)/2)
										                    ,d_intra=rep(1e-5,nspp)
										                    ,d_inter = rep(1e-5,nspp)
										                    ,spp_mat = 1e-5)
										,vars_prior = 2000
										,Zprop_rate = .1
										,distM,link_map
										,dist_mainland , link_mainland
										,int_range_intra="nn",int_range_inter="exp"
										,Z,seed = 42,ini,thin.by = 100,report.by=100){ # ini has same formate of theta
	require(coda)
	source("misc_island.R")
	cat("Initializing...\n\n")
	set.seed(seed)
	nsite = (nrow(distM))
	
	theta = ini
	spp_neig = 1 *( spp_mat!=0 )
	
	nspp = nrow(spp_mat)
	nrep = ncol(Z)
	
    #theta_tuta=ini
	ini = lapply(ini,as.matrix)
	theta.mcmc = list()
	for(i in 1:length(ini)){
	  theta.mcmc[[i]] = mcmc(matrix(nrow = floor(mcmc.iter/thin.by),ncol = length(ini[[i]])),thin = thin.by)
	  
	}
	
	names(theta.mcmc) = names(ini)
	
	
	Z.mcmc = mcmc(matrix(nrow = floor(mcmc.iter/thin.by),ncol = nrow(Z)*nrep),thin = thin.by)
	Z_absolute = (sapply(detmat,rowSums)>0) * 2 - 1
	
	
	Z_curr = Z
	Z_temp = Z
	
	theta_curr = ini
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
	  theta_prop = theta_curr
	  for(j in c(1:length(theta_curr))[-2]){ # no detection proposing
	    theta_prop[[j]] = matrix( rnorm(length(theta_curr[[j]]),mean = 0,sd = sqrt(vars_prop[[j]])),nrow(theta_curr[[j]]),ncol(theta_curr[[j]]) )+ theta_curr[[j]]
	  }
	  
	  theta_prop$spp_mat=theta_prop$spp_mat * spp_neig	#	theta_prop$spp_mat=theta_prop$spp_mat * spp_neig
	  theta_prop$spp_mat = .5*(theta_prop$spp_mat + t( theta_prop$spp_mat)) # must be sym
	  Z_temp = rIsingOccu_multi(theta_prop,X,distM,link_map,dist_mainland , link_mainland,int_range_intra,int_range_inter,n=nrep,method = "CFTP",nIter = 100)
	  # MH ratio
	  
	  Murray_ratio=Murray.ratio(theta_curr ,theta_prop
	                            ,Z_curr ,Z_curr
	                            ,Z_temp
	                            ,detmat
	                            ,vars_prior
	                             
	                            ,X, detX
	                            ,distM,link_map
	                            ,dist_mainland , link_mainland
	                            ,int_range_intra,int_range_inter)
	  r = runif(1)
	  if(Murray_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
	  if(r<=Murray_ratio){
	    theta_curr=theta_prop
	    #Z_temp = Z_temp
	    accept_theta_occu = accept_theta_occu + 1
	  }
	    
		
		theta_prop = theta_curr
		theta_prop[[2]] = matrix( rnorm(length(theta_curr[[2]]),mean = 0,sd = sqrt(vars_prop[[2]])),nrow(theta_curr[[2]]),ncol(theta_curr[[2]]) )+ theta_curr[[2]]
		
		Murray_ratio=Murray.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp
						,detmat
						,vars_prior
						 
						,X, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(Murray_ratio<exp(-10)) low_acc_theta_det = low_acc_theta_det + 1
		if(r<=Murray_ratio){
			theta_curr=theta_prop
			accept_theta_det = accept_theta_det + 1
		}
		
		Z_prop = Z_curr
		
		
		if(runif(1)<Zprop_rate) {
		  flip = sample(which(Z_absolute==-1),1)
			Z_prop[flip]=-Z_prop[flip]
			propose_Z = propose_Z + 1
			
		}
		
		
		Murray_ratio=Murray.ratio(theta_curr ,theta_curr
						,Z_curr ,Z_prop
						,Z_temp
						,detmat
						,vars_prior
						 
						,X, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(Murray_ratio<exp(-10)) low_acc_Z = low_acc_Z + 1
		if(r<=Murray_ratio){
			Z_curr = Z_prop
			accept_Z = accept_Z + 1
		}
		if(i%%report.by == 0) {
		  
		  cat("Burn in iteration: #",i,"\n")
		  cat("# of Z proposed: ",propose_Z,"\n")
		  cat("# of Z acceptance: " , accept_Z-(100-propose_Z),"\n")
		  cat("# of Z acceptance ratio <exp(-10): ",low_acc_Z,"\n\n")
		  cat("# of occupancy theta acceptance: " , accept_theta_occu,"\n")
		  cat("# of occupancy acceptance ratio <exp(-10): ",low_acc_theta_occu,"\n\n")
		  cat("# of detection theta acceptance:" , accept_theta_det,"\n")
		  cat("# of detection acceptance ratio <exp(-10): ",low_acc_theta_det,"\n\n")
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
	
	for(i in 1:(mcmc.iter)){ # for to save 
		#propose theta 
	  theta_prop = theta_curr
	  for(j in c(1:length(theta_curr))[-2]){ # no detection proposing
	    theta_prop[[j]] = matrix( rnorm(length(theta_curr[[j]]),mean = 0,sd = sqrt(vars_prop[[j]])),nrow(theta_curr[[j]]),ncol(theta_curr[[j]]) )+ theta_curr[[j]]
	  }
	  
	  theta_prop$spp_mat=theta_prop$spp_mat * spp_neig	#	theta_prop$spp_mat=theta_prop$spp_mat * spp_neig
	  theta_prop$spp_mat = .5*(theta_prop$spp_mat + t( theta_prop$spp_mat)) # must be sym
		Z_temp = rIsingOccu_multi(theta_prop,X,distM,link_map,dist_mainland , link_mainland,int_range_intra,int_range_inter,n=nrep,method = "CFTP",nIter = 100)

		
		# MH ratio
		Murray_ratio=Murray.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp
						,detmat
						,vars_prior
						 
						,X, detX
						,distM,link_map
						,dist_mainland , link_mainland
						,int_range_intra,int_range_inter)
		r = runif(1)
		if(Murray_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
		if(r<=Murray_ratio){
			theta_curr=theta_prop
			#Z_temp = Z_temp
			accept_theta_occu = accept_theta_occu + 1
		}
		
		
		theta_prop = theta_curr
		theta_prop[[2]] = matrix( rnorm(length(theta_curr[[2]]),mean = 0,sd = sqrt(vars_prop[[2]])),nrow(theta_curr[[2]]),ncol(theta_curr[[2]]) )+ theta_curr[[2]]
		
		Murray_ratio=Murray.ratio(theta_curr ,theta_prop
		                          ,Z_curr ,Z_curr
		                          ,Z_temp
		                          ,detmat
		                          ,vars_prior
		                           
		                          ,X, detX
		                          ,distM,link_map
		                          ,dist_mainland , link_mainland
		                          ,int_range_intra,int_range_inter)
		r = runif(1)
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
		  flip = sample(which(Z_absolute==-1),1)
		  Z_prop[flip]=-Z_prop[flip]
		  propose_Z = propose_Z + 1
		  
		}
		
		
		Murray_ratio=Murray.ratio(theta_curr ,theta_curr
		                          ,Z_curr ,Z_prop
		                          ,Z_temp
		                          ,detmat
		                          ,vars_prior
		                           
		                          ,X, detX
		                          ,distM,link_map
		                          ,dist_mainland , link_mainland
		                          ,int_range_intra,int_range_inter)
		r = runif(1)
		if(Murray_ratio<exp(-10)) low_acc_Z = low_acc_Z + 1
		if(r<=Murray_ratio){
		  #theta_curr=theta_prop
		  Z_curr = Z_prop
		  accept_Z = accept_Z + 1
		}
		
		if(i %% thin.by==0) Z.mcmc[i/thin.by,]=Z_curr
		if(i%%report.by == 0) { # reporting
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
  
	theta.mean =lapply(theta.mcmc,function(thetaa){ apply(thetaa,2,mean)})
  
	res = list(theta.mcmc = theta.mcmc,means = theta.mean,Z.mcmc = Z.mcmc,vars=vars_prop, interaction.range =list( int_range_inter,int_range_intra), envX=X)
	class(res)="IsingOccu_island"
	return(res)
}

