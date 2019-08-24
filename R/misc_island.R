getintralayerGraph = function(distM,link_map,eta,d,int_range = "exp",spp_mat) #it can be used multiple times for interislan and intra-island
{
  # pass all graphs as sparse matrix in package Matrix
  nspp = nrow(spp_mat) # which is the interspecific neighborhood matrix
  A = list() # intralayer graphs are passed using lists
  link_map = as((link_map),"symmetricMatrix")
  if(int_range=="arth"){
    A = lapply(1:nspp,function(i,eta,distM,d){
      eta[i]*as.matrix(1/((distM)^(2+d[i])))
    },eta,distM,d)
  }
  else{
    if(int_range=="exp"){
	  A = lapply(1:nspp,function(i,eta,d,distM,link_map){
		At = eta[i]*as.matrix(exp(-exp(d[i])*distM)) * (link_map)
	    diag(At)=0
	    return(At)
	  },eta,d,distM,link_map)
    }
    else{
      if(int_range=="nn"){
	  A = lapply(1:nspp, function(i,eta,link_map){
		    eta[i]*((link_map))
	    },eta,link_map)
      }
      else{
        #print("int_range must be exp or arth, will assume exp")
		A = lapply(1:nspp,function(i,eta,d,distM,link_map){
		  At = eta[i]*as.matrix(exp(-exp(d[i])*distM)) * (link_map)
	      diag(At)=0
	      return(At)
	    },eta,d,distM,link_map)
      }
    }
  }
  return(A) # if link map is sparse, then A is sparse
} 
  # passed 2019/3/18

getfullGraph = function(A_ex,A_in,spp_mat){
  nspp = nrow(spp_mat)
  nsite = nrow(A_ex[[1]])
  A = Matrix(0,nspp*nsite,nspp*nsite,sparse = T)
  for(i in 2:nspp-1){
    A[1:nsite + (i-1)*nsite,1:nsite + (i-1)*nsite]=A_ex[[1]] + A_in[[1]] # diagonal part
    A_ex[[1]] = NULL # recycle
    A_in[[1]] = NULL
    for(j in (i+1):nspp){
      
      diag(A[1:nsite + (i-1)*nsite,1:nsite + (j-1)*nsite])=spp_mat[i,j]
      diag(A[1:nsite + (j-1)*nsite,1:nsite + (i-1)*nsite])=spp_mat[j,i]
      
    }
  }
  i=nspp
  A[1:nsite + (i-1)*nsite,1:nsite + (i-1)*nsite]=A_ex[[1]] + A_in[[1]]
  rm(A_in,A_ex)
  A = as(A,'symmetricMatrix')
  return(A)
} 
  # passed 2019/3/18

mainland_thr = function(dist_mainland,link_mainland,eta,d,int_range_inter="exp"){
	A = 0*dist_mainland
	link_mainland = (as.matrix(link_mainland))
	if(int_range_inter=="arth"){
			A = eta*as.matrix(1/((dist_mainland)^(2+d)))
	}
	else{
		if(int_range_inter=="exp"){
			A = eta*as.matrix(exp(-exp(d)*dist_mainland)) * (link_mainland)
		}
	  else{
	    if(int_range_inter=="nn")
	    A = eta * (link_mainland)
	  }
	}
	return(A)
	# test for 2spp passed 3/18/2019
}
  # passed 2019/3/18

getMRF = function(theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp"){
  	nsite = nrow(envX)
	beta_occu = theta$beta_occu
	eta_intra = theta$eta_intra # intra spp, intra island if apply
	d_intra = theta$d_intra
	spp_mat = theta$spp_mat
	nspp = nrow(spp_mat)
	#nrep = ncol(Z_vec)
	A_in = getintralayerGraph(distM,link_map$intra,eta_intra,d_intra,int_range = int_range_intra,spp_mat)
	eta_inter = theta$eta_inter # assume there is a 
	d_inter = theta$d_inter
	A_ex = getintralayerGraph(distM,link_map$inter,eta_inter,d_inter,int_range = int_range_inter,spp_mat) # graph among islands, if apply, distM should only contain graph 
	A=getfullGraph(A_ex,A_in,spp_mat)
    rm(A_ex,A_in)
	thr = lapply(1:nspp,
	  function(i,envX,beta_occu,dist_mainland,link_mainland,eta_inter,d_inter,int_range_inter){
	    envX %*% beta_occu[1:ncol(envX)+(i-1)*ncol(envX)] + 
		      mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
	    },envX,beta_occu,dist_mainland,link_mainland,eta_inter,d_inter,int_range_inter)
	thr = Reduce(rbind,thr)
    return(list(A = A,thr = thr))
}

IsingStateProb = function (s, graph, thresholds, beta, responses = c(-1L, 1L)) 
{
  if (!is.list(s)) 
    s <- list(s)
  N <- length(s[[1]])
  Allstates <- do.call(expand.grid, lapply(1:N, function(x) responses))
  Dist <- exp(-beta * apply(Allstates, 1, function(s) H(graph, 
    s, ( thresholds))))
  Z <- sum(Dist)
  sapply(s, function(x) exp(-beta * H(graph, x, ( thresholds)))/Z)
}

Hamiltonian = function(MRF,Z_vec){
	nrep = ncol(Z_vec)
	Ham = lapply(1:nrep,function(i,Z,J,h){H(J,Z[,i],h)},Z=Z_vec,J=MRF$A,h=( MRF$thr))
  
	Ham = Reduce(rbind,Ham)
	return(Ham) # if we have repeat, just make Z_vec has two cols 
	
}

negHamiltonian_posterior = function(theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp",Z_vec){
	H = matrix(0,1,2)
	namestr = H
	beta_occu = theta$beta_occu # this will be a matrix for cols are species
	beta_det = theta$beta_det
	eta_intra = theta$eta_intra # intra spp, intra island if apply
	d_intra = theta$d_intra
	eta_inter = theta$eta_inter # assume there is a 
	d_inter = theta$d_inter
	spp_mat = theta$spp_mat
	eta_inter = theta$eta_inter
	nspp = sqrt(length(spp_mat))
	spp_mat = matrix(spp_mat,nspp,nspp)
	ncov = ncol(envX)
	nsites = nrow(envX)
	
	A = getintralayerGraph(distM,link_map$intra,eta_intra,d_intra,int_range = int_range_intra,spp_mat)
	B = getintralayerGraph(distM,link_map$inter,eta_inter,d_inter,int_range = int_range_inter,spp_mat)
	k = 1
	for(i in 1:nspp){ # intralayer terms:
		#H[[i]]=list()
		for(j in 1:ncov){ # species i, envj
			H[k] = sum(beta_occu[j+(i-1)*ncov]*envX[,j]*Z_vec[1:nsites + (i-1) * nsites,])
			
			namestr[k]=paste0("spp_",i,"_beta_",j)
			k = k + 1
		}
		H[k] = .5 * t(as.matrix( Z_vec[1:nsites + (i-1) * nsites,]))%*%A[[i]]%*%as.matrix( Z_vec[1:nsites + (i-1) * nsites,]) # intra-island 
		
		namestr[k]=paste0("spp_",i,"_Intra_Island")
		k = k + 1
		H[k] = .5 * t(as.matrix( Z_vec[1:nsites + (i-1) * nsites,]))%*%B[[i]]%*%as.matrix( Z_vec[1:nsites + (i-1) * nsites,]) # inter-island
		
		namestr[k]=paste0("spp_",i,"_Inter_Island")
		k = k + 1
		thr_mainland = mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
		H[k] = sum(thr_mainland*Z_vec[1:nsites + (i-1) * nsites,])

		namestr[k]=paste0("spp_",i,"_Mainland")
		k = k + 1
	}
	for(i in 2:nspp-1){
		for(j in (i+1):nspp){
			H[k] = spp_mat[i,j] * (t(Z_vec[1:nsites + (i-1) * nsites,])%*%(Z_vec[1:nsites + (j-1) * nsites,]))
			namestr[k]=paste0("Cor_spp_",i,"_spp_",j)
			k = k + 1
			
		}
	}
	names(H)=namestr
	return(H)
}
  # passed 2019/3/19

rIsingOccu_multi = function(MRF,n=1,method = "CFTP",nIter = 100){
	Z = IsingSamplerCpp(n=n,graph = MRF$A,thresholds=MRF$thr, responses = matrix( c(-1L, 1L),2,1),beta = 1,nIter=nIter,exact = (method=="CFTP"),constrain = NA + MRF$thr)
  return(t(Z))
	# test for 2spp case, passed 3/18/2019
}
  # passed 2019/3/18

Pdet_multi = function(nperiod, envX,detX, beta_det, nspp){ # likelihood given Z and detections If have repeat, use this multiple times.
  if(is.null(detX)) {
    detDesign = lapply(1:nperiod,function(dummy,envX){envX},envX)
    
  }
  else detDesign = lapply(detX,function(x,y){ as.matrix( cbind(y,x))},y = envX)
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
  
Pdet_Ising_single_site = function(thr, Z, dethis, sppmat_det){
	spp_exist = Z==1
  sppmat_det = as(sppmat_det,"dgCMatrix")
	if(sum(spp_exist)==0 | sum(!is.na(dethis))==0){return(0)} # no species there, probability one to be no detection, or no observation here
	if(prod(spp_exist)==0){
	  thr_exis = as.matrix( thr[,spp_exist])
	  # thr_abs = - apply(matrix(sppmat_det[!spp_exist,spp_exist],sum(!spp_exist),sum(spp_exist)),2,sum) # condition on some species not exist here thus never be detected 
	  # do not include thr_abs since if two species never coexist we cannot infer "what if they coexist", i.e. thr_exis will be total colinear with thr_exis
	  thr = thr_exis
	  #thr = apply(matrix(1:ncol(thr_exis)),1,function(k,ww,kk){ww[,k]+kk[k]},thr_exis,( thr_abs))
	}
	graph = sppmat_det[spp_exist,spp_exist]
	has_obs = !(is.na(rowSums(dethis)))
	dethis = dethis[has_obs,spp_exist]# convert it to nrow = nperiod, ncol = nspp for single site, single repeat
	thr = thr[has_obs,]
	
	Pdet_site = apply(matrix(1:sum(has_obs)),1,function(k,dethis,thr,graph){
		IsingStateProb(dethis[k,], graph, thr[k,], beta=1, responses = c(-1L, 1L))
	} ,matrix( dethis,sum(has_obs),sum(spp_exist)), matrix( thr,sum(has_obs),sum(spp_exist)), as( as.matrix(graph),'dsCMatrix'))
	
	return(sum(log(Pdet_site + 1e-15)))
	
}

Sample_Ising_det_single_site = function(thr, Z, dethis, sppmat_det,nIter,n=1, method = "CFTP"){
	spp_exist = Z==1
	dethis[,!spp_exist] = -1# convert it to nrow = nperiod, ncol = nspp for single site, single repeat
	if(sum(spp_exist)==0) return(dethis)
	if(prod(spp_exist)==0){
	  thr_exis = as.matrix( thr[,spp_exist])
	  thr_abs = - apply(matrix(sppmat_det[!spp_exist,spp_exist],sum(!spp_exist),sum(spp_exist)),2,sum) # condition on some species not exist here thus never be detection
	  # check here, may be sth wrong 
	  thr = apply(matrix(1:ncol(thr_exis)),1,function(k,ww,kk){ww[,k]+kk[k]},thr_exis,( thr_abs))
	}
	graph = sppmat_det[spp_exist,spp_exist]
	dethis_exist = dethis[,spp_exist]
	dethis_exist = apply(matrix(1:nrow( as.matrix( dethis))),1,function(k,dethis_exist,thr,graph,nIter,n,method){
		IsingSamplerCpp(n=n,graph = graph, thresholds = thr[k,], beta=1, responses = c(-1L, 1L),nIter = nIter,exact = (method=="CFTP"),constrain = NA+thr[k,])
	},matrix( dethis,sum(has_obs),sum(spp_exist)), matrix( thr,sum(has_obs),sum(spp_exist)), as.matrix( graph),nIter,n,method)
	dethis[,spp_exist] = t(dethis_exist)
	return(dethis)
}

extract_thr = function(i,thr_list){
	nspp = length(thr_list)
	thr = sapply(thr_list,function(thr1,i){t(thr1[i,])},i=i) # thr at site i for all spps, will return a matrix with ncol = nspp, nrow = nperiod
	return(thr)
}

Pdet_Ising = function(nperiod,envX,detX,beta_det,sppmat_det,Z,detmat){
	#require(IsingSamplerCpp)
	 # This is the full design matrix list of detection probability p at time
	if(is.null(detX)) {
	  detDesign = lapply(1:nperiod,function(dummy,envX){envX},envX)
	  
	}
	else detDesign = lapply(detX,function(x,y){ as.matrix( cbind(y,x))},y = envX)
	npardet = ncol(detDesign[[1]])
	nsite = nrow(envX)
	nspp = nrow(sppmat_det)
	thr_list = lapply( 1:nspp, function(i,detDesign,beta_det,naprdet,n_row,nperiod){ 
		temp = lapply(detDesign,function(w,beta1,i){w%*%beta1},beta1 = matrix( beta_det[1:npardet + (i-1) * npardet]),i=i)
		thr1 = (matrix(unlist(temp),nrow = n_row,ncol = nperiod))
		return(thr1) # now here is basically a matrix, for each species at site and period
		},detDesign,beta_det,npardet,nrow(envX),nperiod) # this is gonna be  a list for all species, 
	
	Pdet = lapply(1:nsite,function(i,thr_list,detmat,Z,sppmat_det,nsite,nspp){
		thr1 = extract_thr(i,thr_list)
		rows1 = i + (1:nspp-1)*nsite
		dethis = t(detmat[rows1,])
		Z_site = Z[rows1,]
		Pdet_Ising_single_site(thr1, Z_site, dethis, sppmat_det)
	},thr_list,detmat,as.matrix( Z),sppmat_det,nsite,nspp)# loop over sites
	return(Reduce('+',Pdet))
}

## sampleIsingdet
Sample_Ising_detection = function(nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,nIter=100,n=1, method = "CFTP"){
	#require(IsingSamplerCpp)
  #detDesign = lapply(detX,function(x,y){ as.matrix( cbind(y,x))},y = envX) # This is the full design matrix list of detection probability p at time
  if(is.null(detX)) {
    detDesign = lapply(1:nperiod,function(dummy,envX){envX},envX)
    
  }
  else detDesign = lapply(detX,function(x,y){ as.matrix( cbind(y,x))},y = envX)
  
  npardet = ncol(detDesign[[1]])
  nsite = nrow(envX)
  nspp = nrow(sppmat_det)
  thr_list = lapply( 1:nspp, function(i,detDesign,beta_det,naprdet,n_row,nperiod){ 
    temp = lapply(detDesign,function(w,beta1,i){w%*%beta1},beta1 = matrix( beta_det[1:npardet + (i-1) * npardet]),i=i)
    thr = (matrix(unlist(temp),nrow = n_row,ncol = nperiod))
    return(thr) # now here is basically a matrix, for each species at site and period
  },detDesign,beta_det,npardet,nrow(envX),nperiod) # this is gonna be  a list for all species, 
  
	detmat_list = lapply(1:nsite,function(i,thr_list,detmat,Z,sppmat_det,nsite,nspp,nIter,n, method){
		thr = extract_thr(i,thr_list)
		rows1 = i + (1:nspp-1)*nsite
		dethis = t(detmat[rows1,])
		Z_site = Z[rows1]
		Sample_Ising_det_single_site(thr, Z_site, dethis, sppmat_det,nIter,n, method)
	},thr_list,detmat,Z,sppmat_det,nsite,nspp,nIter,n, method)# loop over sites
	
	det_Ising_spp_list = lapply(1:nspp,function(k,det_list){
	  sapply(det_list,function(sitelist,k){
	    t(sitelist[,k])
	  },k=k)
	},detmat_list)
	detmat = Reduce(cbind,det_Ising_spp_list)
	
	return(t(detmat))
}

Sample_Ising_detection_rep = function(nrep,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,nIter=100,n=1, method = "CFTP"){
  detmat = lapply(1:nrep,function(k,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,nIter,n, method){
    Sample_Ising_detection(nperiod,envX,detX[[k]],beta_det,sppmat_det,Z,detmat[[k]],nIter,n, method)
  },nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,nIter,n, method)
}

Pdet_Ising_rep = function(nrep,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat){
  Pdets = lapply(1:nrep,function(k,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat){
    Pdet_Ising(nperiod,envX,detX[[k]],beta_det,sppmat_det,Z[,k],detmat[[k]])
  },nperiod,envX,detX,beta_det,sppmat_det,Z,detmat)
  return(Reduce('+',Pdets))
}

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

IsingOccu_multi.logL.innorm = function(MRF, Z, beta_det ,detmat, detX,no_obs){ # the in-normalized log likelihood of IsingOccu Model beta is matrix here detX should be a list of list detmat should be a list, they should have the same length
	nspp = nrow(theta$spp_mat)
	beta_det = theta$beta_det
	negPot = Hamiltonian(MRF,Z)
	negpot = -sum(negPot)
	nrep = ncol(Z)
	#beta_det = matrix(detbeta,nrow = length(detbeta),ncol = 1)#
	logLdata=0
	for(i in 1:nrep){
		P_det = Pdet_multi(ncol(detmat[[i]]), envX,detX[[i]], beta_det,nspp) # detX should be a list of list 
		LP_Z1 = as.matrix(rowSums(detmat[[i]] * log(P_det) + (1-detmat[[i]]) * log(1-P_det)))
		LP_Z0 = as.matrix(log(1*(rowSums(detmat[[i]])==0) + 1e-13 * (1-(rowSums(detmat[[i]])==0)))) # I(data = 0), do not want err for those have detections
		temp = as.numeric((Z[,i]+1)/2) * LP_Z1 + as.numeric(1-((Z[,i]+1)/2)) * LP_Z0
		temp[no_obs[,i]] = 0 # not so good, assume no_obs is a matrix with 1 shows there is no camera/camera fail there (no observation)
		logLdata = logLdata +  sum(temp)
	}
	return(negPot+logLdata)
}
  # passed 2019/3/18

IsingOccu_Ising_det_multi_logL_innorm = function(MRF, beta_det,sppmat_det, Z ,detmat, detX){ # the in-normalized log likelihood of IsingOccu Model beta is matrix here detX should be a list of list detmat should be a list, they should have the same length
  nspp = nrow(sppmat_det)
  nperiod = ncol(detmat[[1]])
  negPot = Hamiltonian(MRF,Z)
  negPot = -sum(negPot)
  nrep = ncol(Z)
  logLdata = Pdet_Ising_rep(nrep,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat)
  return(negPot+logLdata)
}

getlogprior = function(theta_prop,theta_curr,vars_prior){
  log_pi_theta_prop = lapply(theta_prop,function(theta_temp,vars_prior){ sum(log(dnorm(as.vector(theta_temp),0,sd=sqrt(vars_prior))))},vars_prior)
  log_pi_theta_prop = sum(unlist(log_pi_theta_prop))
  
  log_pi_theta_curr = lapply(theta_curr,function(theta_temp,vars_prior){ sum(log(dnorm(as.vector( theta_temp),0,sd=sqrt(vars_prior))))},vars_prior)
  log_pi_theta_curr = sum(unlist(log_pi_theta_curr))

  return(list(prop = log_pi_theta_prop,curr = log_pi_theta_curr))

}

Murray.ratio.Ising_det = function(MRF_prop,MRF_curr,log_pi
                        ,Z_curr ,Z_prop
                        ,Z_temp 
                        ,detmat
                        ,beta_det_curr
					    ,beta_det_prop
                        , detX
                        , sppmat_det_curr, sppmat_det_prop){
  log_q_theta_Z_prop_detmat = IsingOccu_Ising_det_multi_logL_innorm(MRF_prop, beta_det_prop, sppmat_det_prop,Z_prop ,detmat, detX)
  # theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
  log_H_theta_curr_Z_temp = -sum(Hamiltonian(MRF_curr,Z_temp))
  
  #### end of the numerator part, start the denominator
  
  log_q_theta_Z_curr_detmat = IsingOccu_Ising_det_multi_logL_innorm(MRF_curr, beta_det_curr, sppmat_det_curr, Z_curr ,detmat = detmat, detX)
  log_H_theta_prop_Z_temp = -sum(Hamiltonian(MRF_prop,Z_temp))
  
  log_MH_ratio = (log_pi$prop + log_q_theta_Z_prop_detmat + log_H_theta_curr_Z_temp)-
    (log_pi$curr + log_q_theta_Z_curr_detmat + log_H_theta_prop_Z_temp)
  
  return(min(1,exp(log_MH_ratio)))
}

write_json.IsingOccu_samples = function(x,path){
  n_sample = nrow(x$Z.mcmc)
  x$theta.mcmc = lapply(x$theta.mcmc,matrix,nrow = n_sample)
  x$Z.mcmc = matrix(x$Z.mcmc,nrow = n_sample)
  class(x) = 'list'
  jsonlite::write_json(x,path)
}