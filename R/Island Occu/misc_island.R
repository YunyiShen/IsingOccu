getintralayerGraph = function(distM,eta,d,int_range = "exp",spp_mat) #it can be used multiple times for interislan and intra-island
{
	#eta = eta[1:nspp]
	nspp = nrow(spp_mat) # which is the interspecific neighborhood matrix
	A = list() # intralayer graphs are passed using lists
	if(int_range=="arth"){
		for(i in 1:nspp){
			A[[i]] = eta[i]*as.matrix(1/((distM)^(2+d[i])))
		}
	}
	else{
		if(int_range=="exp"){
			for(i in 1:nspp){
				A[[i]] = eta[i]*as.matrix(exp(-abs(d[i])*distM))
			}
		}
		else{
			if(int_range=="nn"){
			for(i in 1:nspp){
				A[[i]] = eta[i]*as.matrix((*distM))
			}
			}
			else{
			#print("int_range must be exp or arth, will assume exp")
			for(i in 1:nspp){
				A[[i]] = eta[i]*as.matrix(exp(-abs(d[i])*distM))
			}
			}
		}
	}
	return(A)
}

getfullGraph = function(A_ex,A_in,spp_mat){
	nspp = nrow(spp_mat)
	nsite = nrow(A_ex[[1]])
	A = matrix(0,nspp*nsite,nspp*nsite)
	for(i in 2:nspp-1){
		A[1:nsite + (i-1)*nsite,1:nsite + (i-1)*nsite]=A_ex[[i]] + A_in[[i]] # diagonal part
		for(j in (i+1):nspp){
			
			diag(A[1:nsite + (i-1)*nsite,1:nsite + (j-1)*nsite])=spp_mat[i,j]
			diag(A[1:nsite + (j-1)*nsite,1:nsite + (i-1)*nsite])=spp_mat[j,i]
			
		}
	}
	i=nspp
	A[1:nsite + (i-1)*nsite,1:nsite + (i-1)*nsite]=A_ex[[i]] + A_in[[i]]
	return(A)
}

Hamiltonian = function(theta,envX,distM,distM_island,int_range_intra="nn",int_range_inter="exp",Z_vec){
	beta_occu = theta$beta_occu # this will be a matrix for cols are species
	beta_det = theta$beta_det
	eta_intra = theta$eta_intra # intra spp, intra island if apply
	d_intra = theta$d_intra
	spp_mat = theta$spp_mat
	eta_inter = theta$eta_inter
	nspp = nrow(spp_mat)
	d_inter = theta$d_inter # inter island scaling factor
	nsites = nrow(distM)
	ncov = ncol(envX) # number of covs
	nrep = ncol(Z_vec)
	#zeros = matrix(0,nrow=nsites,ncol=ncov)
	#beta1 = as.numeric( matrix(c(theta[1:(2*ncol(envX))])))
	#Xfull = cbind(rbind(envX,zeros),rbind(zeros,envX))
	thr = X%*%beta_occu # a matrix
	#rm(Xfull)
	A = getintralayerGraph(distM,eta_intra,d,int_range = int_range_intra)
	negPot = matrix(0,1,nrep)
	for(i in 1:nspp){ # intralayer terms:
		negPot = negPot + t(thr[,i])%*%Z_vec[1:nsite + (i-1) * nsite,] + 0.5*t(Z_vec[1:nsite + (i-1) * nsite,])%*%A[[i]]%*%Z_vec[1:nsite + (i-1) * nsite,]
	}
	for(i in 2:nspp-1){
		for (j in (i+1):nspp){
			negPot = negPot + spp_mat[i,j] * t(Z_vec[1:nsite + (i-1) * nsite,])%*%(Z_vec[1:nsite + (j-1) * nsite,])
		}
	}
	if(!is.null(distM_island) & !is.null(theta$eta_inter) & !is.null(int_range_inter) & !is.null(theta$d_inter)){
		eta_inter = theta$eta_inter # assume there is a 
		d_inter = theta$d_inter
		A_inter = getintralayerGraph(distM,eta_inter,d_inter,int_range = int_range_inter) # graph among islands, if apply, distM should only contain graph among different islands, here will be exp for between two island
		for(i in 1:nspp){ # intralayer terms:
			negPot = negPot  + 0.5*t(Z_vec[1:nsite + (i-1) * nsite,])%*%A_inter[[i]]%*%Z_vec[1:nsite + (i-1) * nsite,]
		}
	
	}
	
	return(sum(negPot)) # if we have repeat, just make Z_vec has two cols 
	
}

rIsingOccu_multi = function(theta,envX,distM,distM_island,int_range_intra="nn",int_range_inter="exp",n=1,method = "CFTP",nIter = 100){
	require(IsingSampler)
	beta_occu = theta$beta_occu
	eta_intra = theta$eta_intra # intra spp, intra island if apply
	d_intra = theta$d_intra
	#eta_inter = theta$eta_inter
	spp_mat = theta$spp_mat
	nspp = nrow(spp_mat)
	A_in = getintralayerGraph(distM,eta_intra,d_intra,int_range = int_range_intra)
	if(!is.null(distM_island) & !is.null(theta$eta_inter) & !is.null(int_range_inter) & !is.null(theta$d_inter)){
		eta_inter = theta$eta_inter # assume there is a 
		d_inter = theta$d_inter
		A_ex = getintralayerGraph(distM,eta_inter,d_inter,int_range = int_range_inter) # graph among islands, if apply, distM should only contain graph 
		}
	A=getfullGraph(A_ex,A_in,spp_mat)
	
	thr = envX %*% beta_occu
	thr = matrix(thr,length(thr),1)
	
	Z = IsingSampler(n=n,graph = A,thresholds=thr, responses = c(-1L, 1L),nIter=nIter,method=method)
	return(Z)
}

Pdet_multi = function(detmat, envX,detX, beta_det, nspp){ # likelihood given Z and detections If have repeat, use this multiple times.


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

IsingOccu_multi.logL.innorm = function(theta, envX, distM,distM_island,int_range_intra="nn",int_range_inter="exp", Z ,detmat, detX){ # the in-normalized log likelihood of IsingOccu Model beta is matrix here detX should be a list of list detmat should be a list, they should have the same length
	nspp = nrow(theta$spp_mat)
	beta_det = theta$beta_det
	negPot = Hamiltonian(theta,envX,distM,distM_island,int_range,int_range_island,Z)
	nrep = ncol(Z)
	#beta_det = matrix(detbeta,nrow = length(detbeta),ncol = 1)#
	logLdata=0
	for(i in 1:nrep){
		P_det = Pdet_multi(detmat[[i]], detX[[i]], beta_det,nspp) # detX should be a list of list 
		LP_Z1 = as.matrix(rowSums(detmat[[i]] * log(P_det) + (1-detmat[[i]]) * log(1-P_det)))
		LP_Z0 = as.matrix(log(1*(rowSums(detmat[[i]])==0) + 1e-13 * (1-(rowSums(detmat[[i]])==0)))) # I(data = 0), do not want err for those have detections
		logLdata = logLdata +  sum(as.numeric((Z[,i]+1)/2) * LP_Z1 + as.numeric(1-((Z[,i]+1)/2)) * LP_Z0)
	}
	return(negPot+logLdata)
}

Moller.ratio = function(theta_curr ,theta_prop
						,Z_curr ,Z_prop
						,Z_temp_curr, Z_temp_prop
						,detmat
						,vars_prior
						,theta_tuta
						,envX, detX
						,distM,distM_island
						,int_range_intra="nn",int_range_inter="exp"){
	log_H_theta_tuta_Z_temp_prop = Hamiltonian(theta_tuta,envX,distM,distM_island,int_range_intra,int_range_inter,Z_temp_prop)
	# then auxiliented variable x_prop is same to detmat, together with Z_temp_prop from underlaying Isingmodel. It was proposed using likelihood function with parameter theta_prop and in the main sampler, which is important in canceling out the normalizing constant.
	log_pi_theta_prop =log(dnorm(theta_prop,0,sd=sqrt(vars_prior)))
	log_pi_theta_prop = sum(log_pi_theta_prop)
	#prior of proposed theta
	log_q_theta_Z_prop_detmat = IsingOccu_multi.logL.innorm(theta_prop, envX, distM, distM_island,int_range_intra,int_range_inter,Z_prop ,detmat = detmat, detX)
	# theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
	log_H_theta_Z_temp_curr = Hamiltonian(theta_curr,envX,distM,distM_island,int_range_intra,int_range_inter,Z_temp_curr)

	#### end of the upper part, start the lower

	log_H_theta_tuta_Z_temp_curr = Hamiltonian(theta_tuta,envX,distM,distM_island,int_range_intra,int_range_inter,Z_temp_curr)
	log_pi_theta_curr =log(dnorm(theta_curr,0,sd=sqrt(vars_prior)))
	log_pi_theta_curr = sum(log_pi_theta_curr)
	log_q_theta_Z_curr_detmat = IsingOccu_multi.logL.innorm(theta_curr, envX, distM, distM_island,int_range_intra,int_range_inter,Z_curr ,detmat = detmat, detX)
	log_H_theta_Z_temp_prop = Hamiltonian(theta_prop,envX,distM,distM_island,int_range_intra,int_range_inter,Z_temp_prop)

	log_MH_ratio = (log_H_theta_tuta_Z_temp_prop + log_pi_theta_prop + log_q_theta_Z_prop_detmat + log_H_theta_Z_temp_curr)-
				   (log_H_theta_tuta_Z_temp_curr + log_pi_theta_curr + log_q_theta_Z_curr_detmat + log_H_theta_Z_temp_prop)

	return(min(1,exp(log_MH_ratio)))
}

