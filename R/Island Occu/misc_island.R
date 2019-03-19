getintralayerGraph = function(distM,link_map,eta,d,int_range = "exp",spp_mat) #it can be used multiple times for interislan and intra-island
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
        A[[i]] = eta[i]*as.matrix(exp(-exp(d[i])*distM)) * (link_map)
        diag(A[[i]])=0
      }
    }
    else{
      if(int_range=="nn"){
        for(i in 1:nspp){
          A[[i]] = eta[i]*as.matrix((link_map))
        }
      }
      else{
        #print("int_range must be exp or arth, will assume exp")
        for(i in 1:nspp){
          A[[i]] = eta[i]*as.matrix(exp(-exp(d[i])*distM)) * (link_map)
          diag(A[[i]])=0
        }
      }
    }
  }
  return(A)
} 
  # passed 2019/3/18

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
  # passed 2019/3/18

mainland_thr = function(dist_mainland,link_mainland,eta,d,int_range_inter="exp"){
	A = 0*dist_mainland
	if(int_range_inter=="arth"){
			A = eta*as.matrix(1/((dist_mainland)^(2+d)))
	}
	else{
		if(int_range_inter=="exp"){
			A = eta*as.matrix(exp(-exp(d)*dist_mainland)) * (link_mainland)
		}
	}
	return(A)
	# test for 2spp passed 3/18/2019
}
  # passed 2019/3/18

Hamiltonian = function(theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp",Z_vec){
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
	#thr = envX%*%beta_occu # a matrix
	thr = apply(matrix(1:nspp),1, function(k,beta_occu,envX){ envX %*% beta_occu[1:ncol(envX)+(k-1)*ncol(envX)]},beta_occu,envX)
	#rm(Xfull)
	thr_mainland = 0*thr
	for(i in 1:nspp){
	  thr_mainland[,i ] = mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
	}
	A = getintralayerGraph(distM,link_map$intra,eta_intra,d,int_range = int_range_intra,spp_mat)
	negPot = matrix(0,1,nrep)
	for(i in 1:nspp){ # intralayer terms:
		negPot = negPot + t(as.matrix(thr[,i] ))%*%Z_vec[1:nsites+ (i-1) * nsites,] + 
			apply(Z_vec[1:nsites + (i-1) * nsites,],2,function(Z,A){.5*t(Z)%*%A%*%(Z)},A=A[[i]])
	}
	for(i in 2:nspp-1){
		for (j in (i+1):nspp){
			negPot = negPot + spp_mat[i,j] * diag(t(Z_vec[1:nsites + (i-1) * nsites,])%*%(Z_vec[1:nsites + (j-1) * nsites,]))
		}
	}
	#if(!is.null(link_map$inter) & !is.null(theta$eta_inter) & !is.null(int_range_inter) & !is.null(theta$d_inter)){
	eta_inter = theta$eta_inter # assume there is a 
	d_inter = theta$d_inter
	A_inter = getintralayerGraph(distM,link_map$inter,eta_inter,d_inter,int_range = int_range_inter,spp_mat) # graph among islands, if apply, distM should only contain graph among different islands, here will be exp for between two island
	for(i in 1:nspp){ # intralayer terms:
			thr_mainland = mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
			negPot = negPot  + t(as.matrix(thr_mainland))%*%Z_vec[1:nsites + (i-1) * nsites,] + #mainland part
				apply(Z_vec[1:nsites + (i-1) * nsites,],2,function(Z,A){.5*t(Z)%*%A%*%Z},A=A_inter[[i]])  
				#0.5*t(Z_vec[1:nsite + (i-1) * nsite,])%*%A_inter[[i]]%*%Z_vec[1:nsite + (i-1) * nsite,]
	#	}
	
	}
	
	return(sum(negPot)) # if we have repeat, just make Z_vec has two cols 
	
}
  # passed 2019/3/19

rIsingOccu_multi = function(theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp",n=1,method = "CFTP",nIter = 100){
	require(IsingSampler)
	nsite = nrow(envX)
	beta_occu = theta$beta_occu
	eta_intra = theta$eta_intra # intra spp, intra island if apply
	d_intra = theta$d_intra
	#eta_inter = theta$eta_inter
	spp_mat = theta$spp_mat
	nspp = nrow(spp_mat)
	A_in = getintralayerGraph(distM,link_map$intra,eta_intra,d_intra,int_range = int_range_intra,spp_mat)
	#if(!is.null(link_map$inter) & !is.null(theta$eta_inter) & !is.null(int_range_inter) & !is.null(theta$d_inter)){
		eta_inter = theta$eta_inter # assume there is a 
		d_inter = theta$d_inter
		A_ex = getintralayerGraph(distM,link_map$inter,eta_inter,d_inter,int_range = int_range_inter,spp_mat) # graph among islands, if apply, distM should only contain graph 
	#	}
	#else{
	#	A_ex=0*A_in
		
		
	#	}
	A=getfullGraph(A_ex,A_in,spp_mat)
	
	thr = apply(matrix(1:nspp),2, function(k,beta_occu,envX){ envX %*% beta_occu[1:ncol(envX)+(k-1)*ncol(envX)]},beta_occu,envX)
	#thr = matrix(thr,length(thr),1)
	thr_mainland = 0*thr
	for(i in 1:nspp){
		thr_mainland[1:nsite + (i-1)*nsite] = mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
	}
	
	Z = IsingSampler(n=n,graph = A,thresholds=thr + thr_mainland, responses = c(-1L, 1L),nIter=nIter,method=method)
	return(t(Z))
	# test for 2spp case, passed 3/18/2019
}
  # passed 2019/3/18

Pdet_multi = function(nperiod, envX,detX, beta_det, nspp){ # likelihood given Z and detections If have repeat, use this multiple times.
	# this is still 2 spp case, need to change to multi case
	#nperiod = ncol(detmat) # detmat is the data of 0 and 1 for detections
	# length(beta_det) = 2 * ncol(detX[[1]]) + 2 * ncol(X)  # beta for detections
	detDesign = lapply(detX,function(x,y){ cbind(y,x)},y = envX) # This is the full design matrix list of detection probability p at time
	npardet = ncol(detDesign[[1]])
	#Pdet = list()
	#nrep = length(detX)
	#for(j in 1:nrep){
	  P_det = matrix(0,nrow = 1,ncol = nperiod)
	  P_det = P_det[-1,]
	  for(i in 1:nspp){
		  Xbeta_temp = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = matrix( beta_det[1:npardet + (i-1) * npardet]))
		  P_det_temp = lapply(Xbeta_temp,function(W){exp(W) / (1 + exp(W))}) # just a logistic regression for detection
		  P_det_temp = (matrix(unlist(P_det_temp),nrow = nrow(envX),ncol = nperiod)) # detection probability, row is site i col is period j
		  P_det = rbind(P_det,P_det_temp)
	  }
	  #Pdet[[j]] = P_det
	#}
	return(P_det)
}
  # passed, will return a matrix, with nrow = nspp*nsite, ncol = nperiod,
  #  1:nsite for species 1 and 1:nsite+(i-1)*nsite rows for species i.

Sample_detection = function(nrep,nperiod,envX,detX,beta_det,nspp,Z){
  detmat = list()
  nsite = nrow(envX)
  for(i in 1:nrep){
    r = matrix( runif(nperiod * nspp * nsite) , nspp * nsite,nperiod )
    Pdet = Pdet_multi(nperiod, envX,detX[[i]], beta_det, nspp)
    detmat[[i]] = apply(  1.0 * (r<Pdet),2,function(det,Z){det*Z},Z=(Z[,i]==1) )  
  }
  return(detmat)
}
  # passed 2019/3/19

IsingOccu_multi.logL.innorm = function(theta, envX, distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp", Z ,detmat, detX){ # the in-normalized log likelihood of IsingOccu Model beta is matrix here detX should be a list of list detmat should be a list, they should have the same length
	nspp = nrow(theta$spp_mat)
	beta_det = theta$beta_det
	negPot = Hamiltonian(theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z)
	nrep = ncol(Z)
	#beta_det = matrix(detbeta,nrow = length(detbeta),ncol = 1)#
	logLdata=0
	for(i in 1:nrep){
		P_det = Pdet_multi(ncol(detmat[[i]]), envX,detX[[i]], beta_det,nspp) # detX should be a list of list 
		LP_Z1 = as.matrix(rowSums(detmat[[i]] * log(P_det) + (1-detmat[[i]]) * log(1-P_det)))
		LP_Z0 = as.matrix(log(1*(rowSums(detmat[[i]])==0) + 1e-13 * (1-(rowSums(detmat[[i]])==0)))) # I(data = 0), do not want err for those have detections
		logLdata = logLdata +  sum(as.numeric((Z[,i]+1)/2) * LP_Z1 + as.numeric(1-((Z[,i]+1)/2)) * LP_Z0)
	}
	return(negPot+logLdata)
}
  # passed 2019/3/18

Moller.ratio = function(theta_curr ,theta_prop
						,Z_curr ,Z_prop
						,Z_temp_curr, Z_temp_prop
						,detmat
						,vars_prior
						,theta_tuta
						,envX, detX
						,distM,link_map
						,dist_mainland,link_mainland
						,int_range_intra="nn",int_range_inter="exp"){
	log_H_theta_tuta_Z_temp_prop = Hamiltonian(theta_tuta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp_prop)
	# then auxiliented variable x_prop is same to detmat, together with Z_temp_prop from underlaying Isingmodel. It was proposed using likelihood function with parameter theta_prop and in the main sampler, which is important in canceling out the normalizing constant.
	log_pi_theta_prop = lapply(theta_prop,function(theta_temp,vars_prior){ sum(log(dnorm(theta_temp,0,sd=sqrt(vars_prior))))},vars_prior)
	log_pi_theta_prop = sum(unlist(log_pi_theta_prop))
	#log_pi_theta_prop = sum(log_pi_theta_prop)
	#prior of proposed theta
	log_q_theta_Z_prop_detmat = IsingOccu_multi.logL.innorm(theta_prop, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_prop ,detmat = detmat, detX)
	# theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
	log_H_theta_Z_temp_curr = Hamiltonian(theta_curr,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp_curr)

	#### end of the upper part, start the lower

	log_H_theta_tuta_Z_temp_curr = Hamiltonian(theta_tuta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp_curr)
	log_pi_theta_curr = lapply(theta_curr,function(theta_temp,vars_prior){ sum(log(dnorm(theta_temp,0,sd=sqrt(vars_prior))))},vars_prior)
	log_pi_theta_curr = sum(unlist(log_pi_theta_curr))
	log_q_theta_Z_curr_detmat = IsingOccu_multi.logL.innorm(theta_curr, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_curr ,detmat = detmat, detX)
	log_H_theta_Z_temp_prop = Hamiltonian(theta_prop,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp_prop)

	log_MH_ratio = (log_H_theta_tuta_Z_temp_prop + log_pi_theta_prop + log_q_theta_Z_prop_detmat + log_H_theta_Z_temp_curr)-
				   (log_H_theta_tuta_Z_temp_curr + log_pi_theta_curr + log_q_theta_Z_curr_detmat + log_H_theta_Z_temp_prop)

	return(min(1,exp(log_MH_ratio)))
}
  # passed 2019/3/18

