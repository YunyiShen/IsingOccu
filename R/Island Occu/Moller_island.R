IsingOccu.fit.Moller.sampler = function(X,detmat,detX,spp_neig,nspp=nrow(spp_mat)
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
										,vars_prop = 2,
										distM,distM_island
										,int_range_intra="nn",int_range_inter="exp"
										,Z,seed = 42,ini){ # ini has same formate of theta
	require(coda)
	set.seed(seed)
	nsite = (nrow(distM))
	nspp = nrow(spp_neig)
	
	beta_occu = theta$beta_occu # this will be a matrix for cols are species
	beta_det = theta$beta_det
	eta_intra = theta$eta_intra # intra spp, intra island if apply
	d_intra = theta$d_intra
	eta_inter = theta$eta_inter
	d_inter = theta$d_inter # inter island scaling factor
	spp_mat = theta$spp_mat
	
	nrep = ncol(Z_vec)
	
    theta_tuta=ini
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
	for(i in 1:burn.in){# to burn in 
		#propose theta 
		#theta_prop = rnorm(length(theta_curr),mean = theta_curr,sd = sqrt(vars_prop))
		for(j in 1:5){
			theta_prop[[j]] = rnorm(length(theta_curr[[j]]),mean = theta_curr[[j]],sd = sqrt(vars_prop[[j]]))
		}
		#propose Z from uniform distribution 
		Z_temp_prop = rIsing_multispp(theta_prop,X,distM,int_range,nspp,iter = 300,1)
		# propose x, from the likelihood
		# x_prop = IsingOccu_multispp_sample.detection(theta_prop, X, Z_temp_prop ,detmat, detX,nspp)
		# MH ratio
		Moller_ratio = Moller.ratio_multi(theta_curr ,theta_prop
						,Z_curr ,Z_curr, Z_temp_curr,Z_temp_prop						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range,nspp = nspp)
		r = runif(1)
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			#Z_curr = Z_prop
			# x_curr = x_prop
			Z_temp_curr = Z_temp_prop
		}
		
		
		Z_prop = (Z_absolute==1) + (Z_absolute==-1) * ((runif(length(Z_absolute))>=0.5) * 2 - 1)
		Moller_ratio = Moller.ratio_multi(theta_curr ,theta_curr
						,Z_curr ,Z_prop, Z_temp_curr,Z_temp_curr						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range,nspp = nspp)
		r = runif(1)
		if(r<=Moller_ratio){
			#theta_curr=theta_prop
			Z_curr = Z_prop
			# x_curr = x_prop
			#Z_temp_curr = Z_temp_prop
		}
	}
	
	for(i in 1:(mcmc.save)){ # for to save 
		#propose theta 
		#theta_prop = rnorm(length(theta_curr),mean = theta_curr,sd = sqrt(vars_prop))
		for(j in 1:5){
			theta_prop[[j]] = rnorm(length(theta_curr[[j]]),mean = theta_curr[[j]],sd = sqrt(vars_prop[[j]]))
		}
		
		#propose Z from uniform distribution 
		
		#propose Z_temp from theta_prop's Ising
		Z_temp_prop = rIsing_multispp(theta_prop,X,distM,int_range,nspp,iter = 300,1)
		# propose x, from the likelihood
		# x_prop = IsingOccu_sample.detection(theta_prop, X, Z_temp_prop ,detmat, detX)
		# MH ratio
		Moller_ratio = Moller.ratio_multi(theta_curr ,theta_prop
						,Z_curr ,Z_prop, Z_temp_curr,Z_temp_prop
						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range, nspp)
		r = runif(1)
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			# Z_curr = Z_prop
			# x_curr = x_prop
			Z_temp_curr = Z_temp_prop
		}
		for(j in 1:5){
			theta_mcmc[[j]][i,] =as.vector( theta_curr[[j]])
		}
		Z_prop = (Z_absolute==1) + (Z_absolute==-1) * ((runif(length(Z_absolute))>=0.5) * 2 - 1)
		Moller_ratio = Moller.ratio_multi(theta_curr ,theta_prop
						,Z_curr ,Z_prop, Z_temp_curr,Z_temp_prop
						
						,detmat
						,vars_prior
						,theta_tuta
						,envX=X, detX, distM,int_range, nspp)
		r = runif(1)
		if(r<=Moller_ratio){
			# theta_curr=theta_prop
			Z_curr = Z_prop
			# x_curr = x_prop
			Z_temp_curr = Z_temp_prop
		}
		
		
		
		Z.mcmc[i,]=Z_curr
	}
	
	res = list(theta.mcmc = theta.mcmc,theta.mean = apply(theta.mcmc,1,mean),vars=vars, interaction.range = int_range, graph = graph, envX=X)
	class(res)="IsingOccu_multispp.Moller"
	return(res)
}

