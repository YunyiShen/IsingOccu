# THIS is the IsingOccu fitting function using Moller et al. 2006 sampler (if we can only use MCEM to do MPLE, then Bayesian is much faster)
# remember, X contains 1 col while detX doesn't because the design matrix of det is actually cbind(X,detX)
# detX should be a list, with every element is the design matrix WITHOUT 1s.
# bug here, never accept??, may need to try another way of propose change...
IsingOccu.fit.Moller.sampler = function(X,distM, detmat, detX, mcmc.save = 10000, burn.in = 10 , vars_prior = rep(1,4*ncol(X)+2*ncol(detX[[1]])+9),vars_prop = 2,int_range = "exp",seed = 12345,init,thin.by = 1){
	require(coda)
	require(IsingSampler)
	source("misc.R")
    cat("Initializing...\n\n")
	set.seed(seed)
	nsite = nrow(detmat)/2
	abs_pres = ((rowSums(detmat)>0)) # absolute present
	#datatemp = data.frame(r = rowSums(detmat)>0,rbind(X,X))
	if(missing(init)){
		theta_curr = Initial_MPLE(detmat,envX=X,detX,distM,int_range)
		cat("Initial theta:\n")
		cat(theta_curr,"\n\n")
		
	}
	else theta_curr = init

	ncov = 2*ncol(X)
	ncov_det = 2 * ncol(detX[[1]]) + ncov
	
	
	theta_tuta = theta_curr # improvement possible here, make theta_tuta a better approximation of theta_post
	theta.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(theta_curr)),thin = thin.by)
	p = length(theta_tuta)
	colnames(theta.mcmc)[(p - 4):p] = c("eta_spatial_spc1","d_spatial_spc1","eta_spatial_spc2","d_spatial_spc2","eta_interspecies") # eta spatial
	Z.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = nrow(detmat)),thin = thin.by)
	Z_absolute = as.numeric(abs_pres) * 2 - 1
	# Z_absolute = Z_absolute
	Z_curr = Z_absolute
	Z_temp_curr = rIsingOccu(X,distM,theta = theta_curr,method = "CFTP",nIter = 100,int_range = int_range)

	# try do Z and theta separtly!
	cat("Burn in...\n")
	accept_Z = 0
	accept_theta_occu = 0
	accept_theta_det = 0
	low_acc_Z = 0
	low_acc_theta_occu = 0
	low_acc_theta_det = 0
	timing = proc.time()
	for(i in 1:burn.in){# to burn in
		#propose theta
		#occu theta first
		walk = rnorm(length(theta_curr),0,sd = sqrt(vars_prop))
		theta_prop = theta_curr
		theta_prop[c(1:ncov,1:5+(ncov+ncov_det))] = walk [c(1:ncov,1:5+(ncov+ncov_det))] + theta_prop[c(1:ncov,1:5+(ncov+ncov_det))]
		
		# Aux Z proposed from theta_prop and underlaying graph model to cancel out the normalizing constant
		Z_temp_prop = rIsingOccu(X,distM,theta = theta_prop,method = "CFTP",nIter = 100,int_range = int_range)
		# MH ratio
		Moller_ratio = Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						, Z_temp_curr, Z_temp_prop						
						,detmat
						,vars_prior
						,theta_tuta # this is important, which control the importance sampling 
						,envX=X, detX, distM,int_range)
		if(Moller_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
		r = runif(1)
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			Z_temp_curr = Z_temp_prop
			accept_theta_occu = accept_theta_occu + 1
		}
		
		# propose detection theta 
		theta_prop = theta_curr
		theta_prop[1:ncov_det + ncov] = walk [1:ncov_det + ncov] + theta_prop[1:ncov_det + ncov]
		
		# MH ratio
		Moller_ratio = Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp_curr, Z_temp_curr						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range)
		if(Moller_ratio<exp(-10)) low_acc_theta_det = low_acc_theta_det + 1
		r = runif(1)
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			accept_theta_det = accept_theta_det + 1

		}
				
		# propose Z by single flip, This deals with the unobserved data
		Z_prop = Z_curr
		flip = sample(which(Z_absolute==-1),1)
		Z_prop[flip]=-Z_prop[flip]
		
		Moller_ratio = Moller.ratio(theta_curr ,theta_curr
						,Z_curr ,Z_prop
						,Z_temp_curr, Z_temp_curr						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range)
		if(Moller_ratio<exp(-10)) low_acc_Z = low_acc_Z + 1
		r = runif(1)
		if(r<=Moller_ratio){
			Z_curr = Z_prop
			accept_Z = accept_Z + 1
		}
		
		if(i%%100 == 0) {
		  
		  cat("Burn in iteration: #",i,"\n")
		  cat("# of Z acceptance: " , accept_Z,"\n")
		  cat("# of Z acceptance ratio <exp(-10): ",low_acc_Z,"\n\n")
		  cat("# of occupancy theta acceptance: " , accept_theta_occu,"\n")
		  cat("# of occupancy acceptance ratio <exp(-10): ",low_acc_theta_occu,"\n\n")
		  if(accept_theta_occu==0) cat(theta_curr[c(1:ncov,1:5+(ncov+ncov_det))],"\n\n")
		  cat("# of detection theta acceptance:" , accept_theta_det,"\n")
		  cat("# of detection acceptance ratio <exp(-10): ",low_acc_theta_det,"\n\n")
		  if(accept_theta_det==0) cat(theta_curr[1:ncov_det + ncov],"\n\n")
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
	for(i in 1:(mcmc.save)){ # to save
		#propose theta
		# theta_occu first
		walk = rnorm(length(theta_curr),0,sd = sqrt(vars_prop))
		theta_prop = theta_curr
		theta_prop[c(1:ncov,1:5+(ncov+ncov_det))] = walk [c(1:ncov,1:5+(ncov+ncov_det))] + theta_prop[c(1:ncov,1:5+(ncov+ncov_det))]
		
		# Aux Z proposed from the proposal theta and underlaying graph model
		Z_temp_prop = rIsingOccu(X,distM,theta = theta_prop,method = "CFTP",nIter = 100,int_range = int_range)
		# MH ratio
		Moller_ratio = Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp_curr, Z_temp_prop						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range)
		r = runif(1)
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			Z_temp_curr = Z_temp_prop
			accept_theta_occu = accept_theta_occu + 1
		}
		if(Moller_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
		theta_prop = theta_curr
		theta_prop[1:ncov_det + ncov] = walk [1:ncov_det + ncov] + theta_prop[1:ncov_det + ncov]
		
		# MH ratio
		Moller_ratio = Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp_curr, Z_temp_curr						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range)
		r = runif(1)
		if(Moller_ratio<exp(-10)) low_acc_theta_det = low_acc_theta_det + 1
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			accept_theta_det = accept_theta_det + 1
		}

		# propose Z by single flip, may try Gibbs, this deals with missing obervation of Z
		Z_prop = Z_curr
		flip = sample(which(Z_absolute==-1),1)
		Z_prop[flip]=-Z_prop[flip]
		
		Moller_ratio = Moller.ratio(theta_curr ,theta_curr
						,Z_curr ,Z_prop
						, Z_temp_curr, Z_temp_curr						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range)
		if(Moller_ratio<exp(-10)) low_acc_Z = low_acc_Z + 1
		r = runif(1)
		if(r<=Moller_ratio){
			Z_curr = Z_prop
			accept_Z = accept_Z + 1
		}
		
		theta.mcmc[i,]=theta_curr
		Z.mcmc[i,]=Z_curr
		
		if(i%%100 == 0) { # reporting
		  cat("Sampling iteration: #",i,"\n")
		  cat("# of Z acceptance: " , accept_Z,"\n")
		  cat("# of Z acceptance ratio <exp(-10): ",low_acc_Z,"\n\n")
		  cat("# of occupancy theta acceptance: " , accept_theta_occu,"\n")
		  cat("# of occupancy acceptance ratio <exp(-10): ",low_acc_theta_occu,"\n\n")
		  if(accept_theta_occu==0) cat(theta_curr[c(1:ncov,1:5+(ncov+ncov_det))],"\n\n")
		  cat("# of detection theta acceptance: " , accept_theta_det,"\n")
		  cat("# of detection acceptance ratio <exp(-10): ",low_acc_theta_det,"\n\n")
		  if(accept_theta_det==0) cat(theta_curr[1:ncov_det + ncov],"\n\n")
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
		  }
	}
    theta.mcmc[,c(p-3,p-1)] = abs(theta.mcmc[,c(p-3,p-1)])
	res = list(theta.mcmc = theta.mcmc
	           ,theta.mean = apply(theta.mcmc,2,mean)
	           ,vars_prior=vars_prior
	           ,interaction.range = int_range
	           ,graph = getGraph(distM,apply(theta.mcmc,2,mean),int_range = int_range,full=FALSE), envX=X)
	class(res)="IsingOccu.Moller"
	return(res)
}




