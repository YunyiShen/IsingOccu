getintralayerGraph = function(distM,eta,d,int_range = "exp",nspp)
{
	#eta = eta[1:nspp]
	A = list()
	if(int_range=="arth"){
		for(i in 1:nspp){
			A[i] = eta[i]*as.matrix(1/((distM)^(2+d[i])))
		}
	}
	else{
		if(int_range=="exp"){
			for(i in 1:nspp){
				A[i] = eta[i]*as.matrix(exp(-abs(d[i])*distM))
			}
		}
		else{
			print("int_range must be exp or arth, will assume exp")
			for(i in 1:nspp){
				A[i] = eta[i]*as.matrix(exp(-abs(d[i])*distM))
			}
		}
	}
	return(A)
}

IsingOccu_multispp.logL.innorm = function(theta, envX, distM, Z ,detmat, detX,int_range = "exp",nspp){ # the in-normalized log likelihood of IsingOccu Model beta is matrix here

	beta_occu = theta$beta_occu
	beta_det = theta$beta_det
	eta_intra = theta$eta_intra
	d = theta$d
	eta_inter = theta$eta_inter
# eta_inter should be a symmetric matrix, eta_intra should be vector
# detmat here should be a matrix, with 1:nspp row for first spp, detmat = rbind(detmat1,detmat2,...,etc.)
# beta_det is vector
    #require(IsingSampler)
	#nspp = ncol(Z)
    Z_vec=matrix(Z,nrow=length(Z),ncol=1) # Z can be a vector with 1:nsite be first species.
	
	sites = nrow(distM)
	ncov = ncol(envX)
	zeros = matrix(0,nrow=sites,ncol=ncov)
	#beta1 = as.numeric( matrix(c(theta[1:(2*ncol(envX))])))
	#Xfull = cbind(rbind(envX,zeros),rbind(zeros,envX))
	thr = X%*%beta # a matrix
	#rm(Xfull)
	A = getintralayerGraph(distM,eta_intra,d,int_range = int_range)
	negPot = 0
	for(i in 1:nspp){ # intralayer terms:
		negPot = negPot + t(thr[,i])%*%Z_vec[1:nsite + (i-1) * nsite] + 0.5*t(Z_vec[1:nsite + (i-1) * nsite])%*%A[[i]]%*%Z_vec[1:nsite + (i-1) * nsite]
	}
	for(i in 2:nspp-1){
		for (j in (i+1):nspp){
			negPot = negPot + eta_inter[i,j] * t(Z_vec[1:nsite + (i-1) * nsite])%*%(Z_vec[1:nsite + (j-1) * nsite])
		}
	}
	
	#beta_det = matrix(detbeta,nrow = length(detbeta),ncol = 1)# 
	P_det = Pdet_multi(envX, detmat, detX, beta_det,nspp)
	LP_Z1 = as.matrix(rowSums(detmat * log(P_det) + (1-detmat) * log(1-P_det)))
	LP_Z0 = as.matrix(log(1*(rowSums(detmat)==0) + 1e-13 * (1-(rowSums(detmat)==0)))) # I(data = 0), do not want err for those have detections
	logLdata = sum(as.numeric((Z_vec+1)/2) * LP_Z1 + as.numeric(1-((Z_vec+1)/2)) * LP_Z0)

	return(negPot+logLdata)
}

Pdet_multi = function(envX, detmat, detX, beta_det, nspp) # likelihood given Z and detections
{

	# this is still 2 spp case, need to change to multi case
	nperiod = ncol(detmat) # detmat is the data of 0 and 1 for detections
	# length(beta_det) = 2 * ncol(detX[[1]]) + 2 * ncol(X)  # beta for detections
	detDesign = lapply(detX,function(x,y){cbind(y,x)},y = envX) # This is the full design matrix list of detection probability p at time
	npardet = ncol(detDesign[[1]])
	P_det = matrix(0,nrow = 1,ncol = nperiod)
	P_det = P_det[,-1]
	for(i in 1:nspp){
		Xbeta_temp = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = matrix( beta_det[1:ncol(detX) + (i-1) * ncol(detX)]))
		P_det_temp = lapply(Xbeta_temp,function(W){exp(W) / (1 + exp(W))}) # just a logistic regression for detection
		P_det_temp = (matrix(unlist(P_det_temp),nrow = nrow(X),ncol = nperiod)) # detection probability, row is site i col is period j
		P_det = rbind(P_det,P_det_temp)
	}
	return(P_det)
}


# Moller MH ratio
Moller.ratio = function(theta_curr ,theta_prop
						,Z_curr ,Z_prop
						,x_curr,x_prop
						,detmat
						,vars_prior
						,theta_tuta,Z_tuta
						,envX, detX, distM,int_range ){
	log_q_theta_Z_tuta_x_prop = IsingOccu_multispp.logL.innorm(theta_tuta, envX, distM, Z_tuta ,detmat = x_prop, detX, int_range = int_range) 
	# then auxiliented variable x_prop is same to detmat, and proposed using likelihood function in the main sampler
	log_pi_theta_prop =log(dnorm(theta_prop,0,sd=sqrt(vars_prior)))
	log_pi_theta_prop = sum(log_pi_theta_prop)
	#prior of proposed theta
	log_q_theta_Z_prop_detmat = IsingOccu_multispp.logL.innorm(theta_prop, envX, distM, Z_prop ,detmat = detmat, detX, int_range = int_range)
	# theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
	log_q_theta_Z_curr_x_curr = IsingOccu_multispp.logL.innorm(theta_curr, envX, distM, Z_curr ,detmat = x_curr, detX, int_range = int_range)
	
	#### end of the upper part, start the lower
	
	log_q_theta_Z_tuta_x_curr = IsingOccu_multispp.logL.innorm(theta_tuta, envX, distM, Z_tuta ,detmat = x_curr, detX, int_range = int_range)
	log_pi_theta_curr =log(dnorm(theta_curr,0,sd=sqrt(vars_prior)))
	log_pi_theta_curr = sum(log_pi_theta_curr)
	log_q_theta_Z_curr_detmat = IsingOccu_multispp.logL.innorm(theta_curr, envX, distM, Z_curr ,detmat = detmat, detX, int_range = int_range)
	log_q_theta_Z_prop_x_prop = IsingOccu_multispp.logL.innorm(theta_prop, envX, distM, Z_prop ,detmat = x_prop, detX, int_range = int_range)
	
	log_MH_ratio = (log_q_theta_Z_tuta_x_prop + log_pi_theta_prop + log_q_theta_Z_prop_detmat + log_q_theta_Z_curr_x_curr)-
				   (log_q_theta_Z_tuta_x_curr + log_pi_theta_curr + log_q_theta_Z_curr_detmat + log_q_theta_Z_prop_x_prop)

	return(min(1,exp(log_MH_ratio)))
}

IsingOccu_multispp_sample.detection = function(theta, X, Z ,detmat, detX,nspp){
	Pdet = Pdet_multi(envX, detmat, detX, theta$beta_det, nspp)
	occustatus = matrix(rep((Z+1)/2,nperiod),nrow = length(Z),ncol=nperiod)
	P_det_occu = P_det * occustatus
	RN = matrix( runif(nrow(detmat)*ncol(detmat)),nrow = nrow(detmat),ncol = ncol(detmat) )
	sample_detHistory = 1 * (RN<P_det_occu)
	return(sample_detHistory)
}


# THIS is the IsingOccu fitting function using Moller et al. 2006 sampler (if we can only use MCEM to do MPLE, then Bayesian is much faster)
# remember, X contains 1 col while detX doesn't because the design matrix of det is actually cbind(X,detX)
# detX should be a list, with every element is the design matrix WITHOUT 1s.
IsingOccu.fit.Moller.sampler = function(X,distM, detmat, detX, nspp,mcmc.save = 10000, burn.in = 10 , vars_prior = list(beta_occu = rep(1,ncol(X)),beta_det = rep(1,ncol(detX)),eta_intra = rep(1,nspp),d=rep(1,nspp),eta_inter=rep(1,nspp*(nspp-1)/2)),vars_prop = 2,int_range = "exp",seed = 42){
	require(coda)
	set.seed(seed)
	# starting value has problem, should be multi species
	beta_occu = glm(as.numeric(rowSums(detmat)>0) ~ X - 1, family = binomial)$coef
	beta_det = glm(detmat[,1] ~ cbind(X,detX[[1]]) - 1,family = binomial)$coef
	eta_intra = matrix(1,1,nspp)
	d = matrix(1,1,nspp)
	eta_inter = matrix(1,1,nspp*(nspp-1)/2)
	theta_tuta = list(beta_occu,beta_det,eta_intra,d,eta_inter)
	names(theta_tuta) = ("beta_occu","beta_det","eta_intra","d","eta_inter")
	theta.mcmc = list(
		beta_occu.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(beta_occu))),
		beta_det.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(beta_det))),
		eta_intra.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(eta_intra))),
		d.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(d))),
		eta_inter.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(eta_inter)))
	)
	
	
	Z.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = nrow(detmat)))
	Z_absolute = (as.numeric(rowSums(detmat)>0)) * 2 - 1
	Z_tuta = Z_absolute
	
	
	Z_curr = Z_tuta
	x_curr = detmat
	
	theta_curr = theta_tuta
	for(i in 1:burn.in){# to burn in 
		#propose theta 
		#theta_prop = rnorm(length(theta_curr),mean = theta_curr,sd = sqrt(vars_prop))
		for(j in 1:5){
			theta_prop[[j]] = rnorm(length(theta_curr[[j]]),mean = theta_curr[[j]],sd = sqrt(vars_prop[[j]]))
		}
		#propose Z from uniform distribution 
		Z_prop = (Z_absolute==1) + (Z_absolute==-1) * ((runif(length(Z_absolute))>=0.5) * 2 - 1)
		# propose x, from the likelihood
		x_prop = IsingOccu_multispp_sample.detection(theta_curr, X, Z_curr ,detmat, detX,nspp)
		# MH ratio
		Moller_ratio = Moller.ratio_multi(theta_curr ,theta_prop
						,Z_curr ,Z_prop
						,x_curr,x_prop
						,detmat
						,vars_prior
						,theta_tuta,Z_tuta
						,envX=X, detX, distM,int_range)
		r = runif(1)
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			Z_curr = Z_prop
			x_curr = x_prop
		}
	}
	for(i in 1:(mcmc.save)){ # for to save 
		#propose theta 
		#theta_prop = rnorm(length(theta_curr),mean = theta_curr,sd = sqrt(vars_prop))
		for(j in 1:5){
			theta_prop[[j]] = rnorm(length(theta_curr[[j]]),mean = theta_curr[[j]],sd = sqrt(vars_prop[[j]]))
		}
		#propose Z from uniform distribution 
		Z_prop = (Z_absolute==1) + (Z_absolute==-1) * ((runif(length(Z_absolute))>=0.5) * 2 - 1)
		# propose x, from the likelihood
		x_prop = IsingOccu_sample.detection(theta_curr, X, Z_curr ,detmat, detX)
		# MH ratio
		Moller_ratio = Moller.ratio_multi(theta_curr ,theta_prop
						,Z_curr ,Z_prop
						,x_curr,x_prop
						,detmat
						,vars_prior
						,theta_tuta,Z_tuta
						,envX=X, detX, distM,int_range)
		r = runif(1)
		if(r<=Moller_ratio){
			theta_curr=theta_prop
			Z_curr = Z_prop
			x_curr = x_prop
		}
		for(j in 1:5){
			theta_mcmc[[j]][i,] = theta_curr[[j]]
		}
		Z.mcmc[i,]=Z_curr
	}
	
	res = list(theta.mcmc = theta.mcmc,theta.mean = apply(theta.mcmc,1,mean),vars=vars, interaction.range = int_range, graph = graph, envX=X)
	class(res)="IsingOccu_multispp.Moller"
	return(res)
}

