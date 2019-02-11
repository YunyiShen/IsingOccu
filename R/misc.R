########## helper functions ##########

## get graph of the underlaying graph model
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
	if(int_range=="nn"){
		#cat("assume distM is the graph")
		D1=eta01 * distM
		D2=eta02 * distM
	}
	else{
		print("int_range must be exp or arth, will assume exp")
		D1 = eta01*as.matrix(exp(-abs(d01)*distM))
		D2 = eta02*as.matrix(exp(-abs(d02)*distM))
	
	}
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

## detection probability
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

## Hamiltonian of the underlaying graph model
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
	negPot = t(thr1) %*% Z[1:nsite] + t(thr2) %*%Z [1:nsite+nsite] + .5*t(Z[1:nsite])%*%A$D1%*%Z[1:nsite] + .5*t(Z[1:nsite+nsite])%*%A$D2%*%Z[1:nsite+nsite] + A$eta1* t(Z[1:nsite+nsite]) %*% Z[1:nsite]

	return(negPot)
}

## Total non normalized log likelihood of the graph-detection model
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

## Moller 2006 MH ratio
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

######### data sampling functions ##########

## Take samples from the underlaying graph model
rIsingOccu = function(envX, distM,theta,method = "CFTP",nIter,n=1,int_range = "exp") #distM is a distance matrix
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

## sample detection given Z
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

## helper function for initialization
Init_det = function(beta_det,detmat,envX,detX){
	ncov = ncol(envX)
	ncov_det = ncol(detX[[1]])
	nsite = nrow(envX)
	Z_abs = apply(detmat,1,max)
	#detmat_abs = detmat[Z_abs==1,]
	#detmat_abs2 = detmat[Z_abs[1:nsite + nsite]==1,]
	#n_Z_abs = sum(Z_abs)
	
	P_det = Pdet(envX, detmat, detX, beta_det)
	loglik = as.matrix(rowSums(detmat * log(P_det) + (1-detmat) * log(1-P_det)))
	-sum(loglik[Z_abs==1])
}

# PSEUDO LIKELIHOOD of AUTOOCCU MODEL
# Pseudo Likelihood here, first version with non perfect detection is done Sept/22/2018 and checked NOT CHECKED YET
IsingOccu.logPL = function(theta, envX, distM, Z ,detmat, detX, int_range = "exp") # assume known Z, data should form as nrow(data)==nrow(X),ncol(data)=# of periods. detX is design matrix of detections, WITHOUT 1s, should be a list. length(detX) = ncol(data) = # of periods # A1 is dist matrix for species 2, all about first spc1 should be 0
{
	# deal with the underlaying Ising model
	# REMEMBER TO ADD 1s to envX, detX will automaticlly include 1s because envX in cbinded
	# nsite = nrow(X)/2	
	require(IsingSampler)
	p = length(theta)
	sites = nrow(distM)
	ncov = ncol(envX)
	ncov_det = ncol(detX[[1]])
	#zeros = matrix(0,nrow=sites,ncol=ncov)
	beta1 = as.numeric( matrix(c(theta[1:(2*ncov)])))
	thr = rbind(envX %*% beta1[1:ncov],envX %*% beta1[1:ncov+ncov])
	#rm(Xfull)
	#A = getGraph(distM,theta,int_range = int_range)
	log_PL = logPL(theta,Z,envX,distM,int_range)
	#rm(A)
	
	P_det = Pdet(envX, detmat, detX, theta[1:(2*ncov_det + 2*ncov)+2*ncov])
	LP_Z1 = as.matrix(rowSums(detmat * log(P_det) + (1-detmat) * log(1-P_det)))
	LP_Z0 = as.matrix(log(1*(rowSums(detmat)==0) + 1e-13 * (1-(rowSums(detmat)==0)))) # I(data = 0), do not want err for those have detections
	logLdata = sum(as.numeric((Z+1)/2) * LP_Z1 + as.numeric(1-((Z+1)/2)) * LP_Z0) # likelihood of data 
	
	# total neg log likelihood
	log_PL-logLdata 
}


logPL = function(theta,Z,envX,distM,int_range){
	Z = matrix(Z)
	nsite = nrow(envX)
	ncov = ncol(envX)
	A = getGraph(distM,theta,int_range = int_range,full=FALSE)
	Xbeta1 = envX %*% matrix( theta[1:ncov])
	Xbeta2 = envX %*% matrix(theta[1:ncov + ncov])
	
	logPL1 = Xbeta1 + A$D1%*%Z[1:nsite] + A$eta1*Z[1:nsite + nsite]
	logPL1 = t(Z[1:nsite]) %*% logPL1 -sum(log(exp(-logPL1)+exp(logPL1)))
	logPL2 = Xbeta2 + A$D2%*%Z[1:nsite+nsite] + A$eta1*Z[1:nsite]
	logPL2 = t(Z[1:nsite + nsite]) %*% logPL2 -sum(log(exp(-logPL2)+exp(logPL2)))
	log_PL = -logPL1-logPL2
	return(log_PL)



}

Initial_MPLE = function(detmat,envX,detX,distM,int_range){
	ncov = ncol(envX)
	ncov_det = ncol(detX[[1]])
	beta_det_ini = runif(2*(ncov + ncov_det))
	beta_det_ini = optim(beta_det_ini,Init_det,detmat = detmat,envX = envX,detX = detX,method = "SANN")$par
	P_det = Pdet(envX, detmat, detX, beta_det_ini)
	no_det = rowSums(log(1-P_det))
	Z_abs = apply(detmat,1,max)
	thr = min(no_det[Z_abs==1])
	Z_abs[Z_abs==0] = no_det[Z_abs==0]<=thr
	Z_abs = Z_abs * 2 -1
	
	theta_occu_ini = rnorm(4*ncov + 2*ncov_det+5)
	#theta_occu_ini = optim(theta_occu_ini,IsingOccu.logPL,Z=Z_abs,detmat = detmat,envX = envX,distM=distM,int_range = int_range)
	theta_occu_ini = optim(par=(theta_occu_ini),fn=IsingOccu.logPL,envX=envX,distM=distM,Z=Z_abs,detmat=detmat,detX=detX,int_range = int_range,method = "SANN")

  
	return(theta_occu_ini$par)

}