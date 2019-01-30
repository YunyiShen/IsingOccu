
getGraph = function(distM,theta,int_range = "exp",full=TRUE)
{
	p = length(theta)
	sites = nrow(distM)
	eta01 = theta[p-4]
	d01 = theta[p-3]
	eta02 = theta[p-2]
	d02 = theta[p-1]
	eta1 = theta[p]
	if(int_range=="arth"){
		D1 = eta01*as.matrix(1/((distM)^(2+d01))) # long range interactions dimension constant d01 and d02 (sigma in some physics papers regarding Ising model)
		D2 = eta02*as.matrix(1/((distM)^(2+d02)))
	}
	else{
	if(int_range=="exp"){
		D1 = eta01*as.matrix(exp(-abs(d01)*distM))
		D2 = eta02*as.matrix(exp(-abs(d02)*distM))
	}
	else{
	print("int_range must be exp or arth, will assume exp")
		D1 = eta01*as.matrix(exp(-abs(d01)*distM))
		D2 = eta02*as.matrix(exp(-abs(d02)*distM))
		}
	}
	if(full){
		I = matrix(0,nrow=sites,ncol=sites)
		diag(I) = eta1
		A = cbind(D1,I)
		B = cbind(I,D2)
		A = rbind(A,B)
		rm(B)
		diag(A)=0
		row.names(A)=colnames(A)
	}
	else{
		diag(D1)=0
		diag(D2)=0
		A=list(D1=D1,D2=D2,eta1=eta1)
	}
	return(A)
}

rIsingOccu = function(envX, distM,theta,method = "MH",nIter=nIter,n=1,int_range = "exp") #distM is a distance matrix
{
	require(IsingSampler)
	p = length(theta)
	sites = nrow(distM)
	ncov = ncol(envX)
	zeros = matrix(0,nrow=sites,ncol=ncov)
	beta1 = as.numeric( matrix(c(theta[1:(ncol(envX))])))
	beta2 = as.numeric(matrix(theta[1:ncol(envX) + ncol(envX)]))
	#Xfull = cbind(rbind(envX,zeros),rbind(zeros,envX))
	thr = rbind(envX %*% beta1,envX %*% beta2)
	#rm(Xfull)
	A = getGraph(distM,theta,int_range = int_range)
	as.matrix(IsingSampler(n=n,graph = A,thresholds=thr, responses = c(-1L, 1L),nIter=nIter,method=method)) # use IsingSampler
}

# sample detection given Z
IsingOccu_sample.detection = function(theta, envX, Z ,detmat, detX){
	p = length(theta)
	RN = matrix( runif(nrow(detmat)*ncol(detmat)),nrow = nrow(detmat),ncol = ncol(detmat) )
	beta_det = theta[(2*ncol(envX) + 1):(p-5)] # length(beta_det) = 2 * ncol(detX[[1]]) + 2 * ncol(X)  # beta for detections
	P_det = Pdet(envX, detmat, detX, beta_det)
	# occustatus = matrix(rep((Z+1)/2,nperiod),nrow = length(Z),ncol=nperiod)
	P_det_occu = apply(P_det,2,function(X,Y){X*Y},Y = (Z+1)/2)
	# P_det_occu = P_det * occustatus

	sample_detHistory = 1 * (RN<P_det_occu)
	sample_detHistory
}


Pdet = function(envX, detmat, detX, beta_det) # likelihood given detections if occupy
{
	nperiod = ncol(detmat) # detmat is the data of 0 and 1 for detections
	#beta_det = theta[(2*ncol(envX) + 1):(p-5)] # length(beta_det) = 2 * ncol(detX[[1]]) + 2 * ncol(X)  # beta for detections
	if(is.null(detX)){
	  detDesign = envX
	  P_det1 = envX %*% beta_det[1:npardet]
	  P_det2 = envX %*% beta_det[1:npardet + npardet]
	  P_det = rbind(P_det1,P_det2)
	  P_det = matrix(P_det,nrow = nrow(envX),ncol = nperiod)
	  
	  }
	else{
	  detDesign = lapply(detX,function(x,y){cbind(y,x)},y = envX) # This is the full design matrix list of detection probability p at time
	  npardet = ncol(detDesign[[1]])
	  Xbeta_det1 = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = beta_det[1:npardet]) # Xbeta for detections
	  Xbeta_det2 = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = beta_det[1:npardet + npardet])
	  P_det1 = lapply(Xbeta_det1,function(W){exp(W) / (1 + exp(W))}) # just a logistic regression for detection
	  rm(Xbeta_det1)
	  P_det2 = lapply(Xbeta_det2,function(W){exp(W) / (1 + exp(W))})
	  rm(Xbeta_det2)
	  P_det1 = (matrix(unlist(P_det1),nrow = nrow(envX),ncol = nperiod)) # detection probability, row is site i col is period j
	  P_det2 = (matrix(unlist(P_det2),nrow = nrow(envX),ncol = nperiod))
	  P_det = rbind(P_det1,P_det2)
	 }
	return(P_det)
}
# This has no problem 20190118


Hamiltonian = function(theta, envX, distM, Z ,int_range = "exp"){
	Z=matrix(Z,nrow=length(Z),ncol=1) # make Z col vector
	p = length(theta)
	nsite = nrow(distM)
	ncov = ncol(envX)
	
	# zeros = matrix(0,nrow=nsite,ncol=ncov)
	beta1 = as.numeric( matrix(c(theta[1:(ncov)])))
	beta2 = as.numeric( matrix(c(theta[1:(ncov) + ncov])))
	# Xfull = cbind(rbind(envX,zeros),rbind(zeros,envX))
	thr1 = envX%*%beta1
	thr2 = envX%*%beta2
	# rm(Xfull)
	A = getGraph(distM,theta,int_range = int_range,full=FALSE)
	negPot = t(thr1)%*%Z[1:nsite] + t(thr2)%*%Z[1:nsite+nsite] + 0.5*t(Z[1:nsite])%*%A$D1%*%Z[1:nsite] + 0.5*t(Z[1:nsite+nsite])%*%A$D2%*%Z[1:nsite+nsite] + A$eta1*t(Z[1:nsite+nsite])%*%Z[1:nsite]

	return(negPot)


}

IsingOccu.logL.innorm = function(theta, envX, distM, Z ,detmat, detX, int_range = "exp"){ # the in-normalized log likelihood of IsingOccu Model
    # require(IsingSampler)
    
	p = length(theta)
	H = Hamiltonian(theta, envX, distM, Z ,int_range)

	beta_det = theta[(2*ncol(envX) + 1):(p-5)] # length(beta_det) = 2 * ncol(detX[[1]]) + 2 * ncol(X)  # beta for detections
	P_det = Pdet(envX, detmat, detX, beta_det)
	LP_Z1 = as.matrix(rowSums(detmat * log(P_det) + (1-detmat) * log(1-P_det)))
	LP_Z0 = as.matrix(log(1*(rowSums(detmat)==0) + 1e-13 * (1-(rowSums(detmat)==0)))) # I(data = 0), do not want err for those have detections
	logLdata = sum(as.numeric((Z+1)/2) * LP_Z1 + as.numeric(1-(Z+1)/2) * LP_Z0)

	return(H+logLdata)
}

# Moller MH ratio
# assume to be good 20190119
Moller.ratio = function(theta_curr ,theta_prop
						,Z_curr ,Z_prop
						,Z_temp_curr, Z_temp_prop
						,detmat
						,vars_prior
						,theta_tuta
						,envX, detX, distM,int_range ){
	log_H_theta_tuta_Z_temp_prop = Hamiltonian(theta_tuta, envX, distM, Z_temp_prop , int_range = int_range)
	# then auxiliented variable x_prop is same to detmat, together with Z_temp_prop from underlaying Isingmodel. It was proposed using likelihood function with parameter theta_prop and in the main sampler, which is important in canceling out the normalizing constant.
	log_pi_theta_prop =log(dnorm(theta_prop,0,sd=sqrt(vars_prior)))
	log_pi_theta_prop = sum(log_pi_theta_prop)
	#prior of proposed theta
	log_q_theta_Z_prop_detmat = IsingOccu.logL.innorm(theta_prop, envX, distM, Z_prop ,detmat = detmat, detX, int_range = int_range)
	# theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
	log_H_theta_Z_temp_curr = Hamiltonian(theta_curr, envX, distM, Z_temp_curr , int_range = int_range)

	#### end of the upper part, start the lower

	log_H_theta_tuta_Z_temp_curr = Hamiltonian(theta_tuta, envX, distM, Z_temp_curr , int_range = int_range)
	log_pi_theta_curr =log(dnorm(theta_curr,0,sd=sqrt(vars_prior)))
	log_pi_theta_curr = sum(log_pi_theta_curr)
	log_q_theta_Z_curr_detmat = IsingOccu.logL.innorm(theta_curr, envX, distM, Z_curr ,detmat = detmat, detX, int_range = int_range)
	log_H_theta_Z_temp_prop = Hamiltonian(theta_prop, envX, distM, Z_temp_prop , int_range = int_range)

	log_MH_ratio = (log_H_theta_tuta_Z_temp_prop + log_pi_theta_prop + log_q_theta_Z_prop_detmat + log_H_theta_Z_temp_curr)-
				   (log_H_theta_tuta_Z_temp_curr + log_pi_theta_curr + log_q_theta_Z_curr_detmat + log_H_theta_Z_temp_prop)

	return(min(1,exp(log_MH_ratio)))
}



# THIS is the IsingOccu fitting function using Moller et al. 2006 sampler (if we can only use MCEM to do MPLE, then Bayesian is much faster)
# remember, X contains 1 col while detX doesn't because the design matrix of det is actually cbind(X,detX)
# detX should be a list, with every element is the design matrix WITHOUT 1s.
# bug here, never accept??, may need to try another way of propose change...
IsingOccu.fit.Moller.sampler = function(X,distM, detmat, detX, mcmc.save = 10000, burn.in = 10 , vars_prior = rep(1,4*ncol(X)+2*ncol(detX[[1]])+9),vars_prop = 2,int_range = "exp",seed = 12345,init){
	require(coda)
  cat("Initializing...\n\n")
	set.seed(seed)
	nsite = nrow(detmat)/2
	abs_pres = ((rowSums(detmat)>0)) # absolute present
	#datatemp = data.frame(r = rowSums(detmat)>0,rbind(X,X))
	if(missing(init)){
	start = c(glm.fit(X,abs_pres[1:nsite],family = binomial())$coef
	          ,glm.fit(X,abs_pres[1:nsite+nsite],family = binomial())$coef)
	cat("Initial occupancy theta:\n")
	cat(start,"\n\n")
	start_det = c( glm.fit(cbind(X,detX[[1]])[which(abs_pres[1:nsite]),] , detmat[which(abs_pres[1:nsite]),1],family = binomial())$coef
	               , glm.fit(cbind(X,detX[[1]])[which(abs_pres[1:nsite+nsite]),] , detmat[which(abs_pres[1:nsite+nsite]),1],family = binomial())$coef)
	cat("Initial detection theta:\n")
	cat(start_det,"\n\n")
	ncov = length(start)
	ncov_det = length(start_det)
	theta_curr = matrix( c(start,start_det,1,1,1,1,-1))
	cat("Initial interactions:\n")
	cat("eta_sigma:1 eta_tau:1 d_sigma:1 d_tau:1 eta:-1\n\n\n")
	}
	else theta_curr = init
	# rm(datatemp)
	
	
	
	theta_tuta = theta_curr
	theta.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(theta_curr)))
	p = length(theta_tuta)
	colnames(theta.mcmc)[(p - 4):p] = c("eta_spatial_spc1","d_spatial_spc1","eta_spatial_spc2","d_spatial_spc2","eta_interspecies") # eta spatial
	Z.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = nrow(detmat)))
	Z_absolute = as.numeric(abs_pres) * 2 - 1
	# Z_absolute = Z_absolute
	Z_curr = Z_absolute
	Z_temp_curr = rIsingOccu(X,distM,theta = theta_curr,method = "MH",nIter = 300,int_range = int_range)
		# propose x, from the likelihood using proposed theta and sampled Z
	#x_curr = IsingOccu_sample.detection(theta_curr, X, Z_temp_curr ,detmat, detX)
	# try do Z and theta separtly!
	cat("Burn in...\n")
	accept_Z = 0
	accept_theta_occu = 0
	accept_theta_det = 0
	for(i in 1:burn.in){# to burn in
		#propose theta
		# theta first
		walk = rnorm(length(theta_curr),0,sd = sqrt(vars_prop))
		theta_prop = theta_curr
		theta_prop[c(1:ncov,1:5+(ncov+ncov_det))] = walk [c(1:ncov,1:5+(ncov+ncov_det))] + theta_prop[c(1:ncov,1:5+(ncov+ncov_det))]
		# theta_prop = rnorm(length(theta_curr),mean = theta_curr,sd = sqrt(vars_prop))
		
		# Aux Z
		Z_temp_prop = rIsingOccu(X,distM,theta = theta_prop,method = "MH",nIter = 300,int_range = int_range)
		# propose x, from the likelihood using proposed theta and sampled Z
		# x_prop = IsingOccu_sample.detection(theta_prop, X, Z_temp_prop ,detmat, detX)
		# MH ratio
		Moller_ratio = Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						, Z_temp_curr, Z_temp_prop						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range)
		r = runif(1)
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			#Z_curr = Z_prop
			Z_temp_curr = Z_temp_prop
			#x_curr = x_prop
			accept_theta_occu = accept_theta_occu + 1
			#cat("accepted!\n")
		}
		
		theta_prop = theta_curr
		theta_prop[1:ncov_det + ncov] = walk [1:ncov_det + ncov] + theta_prop[1:ncov_det + ncov]
		# theta_prop = rnorm(length(theta_curr),mean = theta_curr,sd = sqrt(vars_prop))
		
		# Aux Z
		# Z_temp_prop = rIsingOccu(X,distM,theta = theta_prop,method = "MH",nIter = 300,int_range = int_range)
		# propose x, from the likelihood using proposed theta and sampled Z
		# x_prop = IsingOccu_sample.detection(theta_prop, X, Z_temp_prop ,detmat, detX)
		# MH ratio
		Moller_ratio = Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp_curr, Z_temp_curr						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range)
		r = runif(1)
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			#Z_curr = Z_prop
			#Z_temp_curr = Z_temp_prop
			#x_curr = x_prop
			accept_theta_det = accept_theta_det + 1
			#cat("accepted!\n")
		}
		
		
		
		
		
		# propose Z by single flip, may try Gibbs
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
		r = runif(1)
		if(r<=Moller_ratio){
			#theta_curr=theta_prop
			Z_curr = Z_prop
			#Z_temp_curr = Z_temp_prop
			#x_curr = x_prop
			accept_Z = accept_Z + 1
			#cat("accepted!\n")
		}
		
		if(i%%100 == 0) {
		  cat("Burn in iteration: #",i,"\n")
		  cat("# of Z acceptance:" , accept_Z,"\n")
		  cat("# of occupancy theta acceptance:" , accept_theta_occu,"\n")
		  if(accept_theta_occu==0) cat(theta_curr[c(1:ncov,1:5+(ncov+ncov_det))])
		  cat("# of detection theta acceptance:" , accept_theta_det,"\n")
		  if(accept_theta_det==0) cat(theta_curr[1:ncov_det + ncov])
		  
		  accept_Z = 0
		  accept_theta_occu = 0
		  accept_theta_det = 0
		  }
	}
	cat("Start sampling...\n")
	accept_Z = 0
	accept_theta_occu = 0
	accept_theta_det = 0
	for(i in 1:(mcmc.save)){ # to save
		#propose theta
		# theta first
		walk = rnorm(length(theta_curr),0,sd = sqrt(vars_prop))
		theta_prop = theta_curr
		theta_prop[c(1:ncov,1:5+(ncov+ncov_det))] = walk [c(1:ncov,1:5+(ncov+ncov_det))] + theta_prop[c(1:ncov,1:5+(ncov+ncov_det))]
		# theta_prop = rnorm(length(theta_curr),mean = theta_curr,sd = sqrt(vars_prop))
		
		# Aux Z
		Z_temp_prop = rIsingOccu(X,distM,theta = theta_prop,method = "MH",nIter = 300,int_range = int_range)
		# propose x, from the likelihood using proposed theta and sampled Z
		# x_prop = IsingOccu_sample.detection(theta_prop, X, Z_temp_prop ,detmat, detX)
		# MH ratio
		Moller_ratio = Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						, Z_temp_curr, Z_temp_prop						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range)
		r = runif(1)
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			#Z_curr = Z_prop
			Z_temp_curr = Z_temp_prop
			#x_curr = x_prop
			accept_theta_occu = accept_theta_occu + 1
			#cat("accepted!\n")
		}
		
		theta_prop = theta_curr
		theta_prop[1:ncov_det + ncov] = walk [1:ncov_det + ncov] + theta_prop[1:ncov_det + ncov]
		# theta_prop = rnorm(length(theta_curr),mean = theta_curr,sd = sqrt(vars_prop))
		
		# Aux Z
		# Z_temp_prop = rIsingOccu(X,distM,theta = theta_prop,method = "MH",nIter = 300,int_range = int_range)
		# propose x, from the likelihood using proposed theta and sampled Z
		# x_prop = IsingOccu_sample.detection(theta_prop, X, Z_temp_prop ,detmat, detX)
		# MH ratio
		Moller_ratio = Moller.ratio(theta_curr ,theta_prop
						,Z_curr ,Z_curr
						,Z_temp_curr, Z_temp_curr						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range)
		r = runif(1)
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			#Z_curr = Z_prop
			#Z_temp_curr = Z_temp_prop
			#x_curr = x_prop
			accept_theta_det = accept_theta_det + 1
			#cat("accepted!\n")
		}

		# propose Z by single flip, may try Gibbs
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
		r = runif(1)
		if(r<=Moller_ratio){
			#theta_curr=theta_prop
			Z_curr = Z_prop
			#Z_temp_curr = Z_temp_prop
			#x_curr = x_prop
			accept_Z = accept_Z + 1
			#cat("accepted!\n")
		}
		
		theta.mcmc[i,]=theta_curr
		Z.mcmc[i,]=Z_curr
		
		if(i%%100 == 0) {
		  cat("Sampling iteration: #",i,"\n")
		  cat("# of Z acceptance:" , accept_Z,"\n")
		  cat("# of occupancy theta acceptance:" , accept_theta_occu,"\n")
		  if(accept_theta_occu==0) cat(theta_curr[c(1:ncov,1:5+(ncov+ncov_det))])
		  cat("# of detection theta acceptance:" , accept_theta_det,"\n")
		  if(accept_theta_det==0) cat(theta_curr[1:ncov_det + ncov])
		  
		  accept_Z = 0
		  accept_theta_occu = 0
		  accept_theta_det = 0
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




