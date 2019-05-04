getintralayerGraph = function(distM,link_map,eta,d,int_range = "exp",spp_mat) #it can be used multiple times for interislan and intra-island
{
  #eta = eta[1:nspp]
  nspp = nrow(spp_mat) # which is the interspecific neighborhood matrix
  A = list() # intralayer graphs are passed using lists
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
		    eta[i]*as.matrix((link_map))
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
	  else{
	    if(int_range_inter=="nn")
	    A = eta * as.matrix(link_mainland)
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
	#thr = apply(matrix(1:nspp),1, function(k,beta_occu,envX){ envX %*% beta_occu[1:ncol(envX)+(k-1)*ncol(envX)]},beta_occu,envX)
	
	#rm(Xfull)
	#thr_mainland = 0*thr
	A = getintralayerGraph(distM,link_map$intra,eta_intra,d_intra,int_range = int_range_intra,spp_mat)
	negPot = matrix(0,1,nrep)
	for(i in 1:nspp){ # intralayer terms:
	  thr = envX %*% beta_occu[1:ncol(envX)+(i-1)*ncol(envX)]
		negPot = negPot + t(as.matrix(thr ))%*%Z_vec[1:nsites + (i-1) * nsites,] + 
			apply(as.matrix(Z_vec[1:nsites + (i-1) * nsites,]),2,function(Z,A){.5*t(Z)%*%A%*%(Z)},A=A[[i]])
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
	for(i in 1:nspp){ # intralayer, inter island terms:
			thr_mainland = mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
			negPot = negPot  + t(as.matrix(thr_mainland))%*%Z_vec[1:nsites + (i-1) * nsites,] + #mainland part
				apply(as.matrix( Z_vec[1:nsites + (i-1) * nsites,]),2,function(Z,A){.5*t(Z)%*%A%*%Z},A=A_inter[[i]])  
				#0.5*t(Z_vec[1:nsite + (i-1) * nsite,])%*%A_inter[[i]]%*%Z_vec[1:nsite + (i-1) * nsite,]
	#	}
	
	}
	
	return(sum(negPot)) # if we have repeat, just make Z_vec has two cols 
	
}

Hamiltonian_posterior = function(theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp",Z_vec){
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
	eta_inter = theta$eta_inter # assume there is a 
	d_inter = theta$d_inter
	A_ex = getintralayerGraph(distM,link_map$inter,eta_inter,d_inter,int_range = int_range_inter,spp_mat) # graph among islands, if apply, distM should only contain graph 
	A=getfullGraph(A_ex,A_in,spp_mat)
	thr = matrix(0,nspp*ncol(envX))
	#thr = apply(matrix(1:nspp),1, function(k,beta_occu,envX){ envX %*% beta_occu[1:ncol(envX)+(k-1)*ncol(envX)]},beta_occu,envX)
	#thr = matrix(thr,length(thr),1)
	thr_mainland = thr
	for(i in 1:nspp){
	  thr[1:nsite + (i-1)*nsite] = envX %*% beta_occu[1:ncol(envX)+(i-1)*ncol(envX)]
		thr_mainland[1:nsite + (i-1)*nsite] = mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
	}
	
	Z = IsingSampler(n=n,graph = A,thresholds=thr + thr_mainland, responses = c(-1L, 1L),nIter=nIter,method=method,CFTPretry = 1)
	return(t(Z))
	# test for 2spp case, passed 3/18/2019
}
  # passed 2019/3/18

Pdet_multi = function(nperiod, envX,detX, beta_det, nspp){ # likelihood given Z and detections If have repeat, use this multiple times.
	# this is still 2 spp case, need to change to multi case
	#nperiod = ncol(detmat) # detmat is the data of 0 and 1 for detections
	# length(beta_det) = 2 * ncol(detX[[1]]) + 2 * ncol(X)  # beta for detections
	#detDesign = lapply(detX,function(x,y){ as.matrix( cbind(y,x))},y = envX) # This is the full design matrix list of detection probability p at time
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
	if(sum(spp_exist)==0){ return(0)} # no species there, probability one to be no detection
	graph = sppmat_det[spp_exist,spp_exist]
	thr_exis = as.matrix( thr[,spp_exist])
	thr_abs = - apply(as.matrix(sppmat_det[!spp_exist,spp_exist]),2,sum) # condition on some species not exist here thus never be detection
	thr = apply(matrix(1:ncol(thr_exis)),1,function(k,ww,kk){ww[,k]+kk[k]},thr_exis,( thr_abs))
	#thr = t(thr)
	dethis = dethis[,spp_exist]# convert it to nrow = nperiod, ncol = nspp for single site, single repeat
	
	Pdet_site = apply(matrix(1:nrow(as.matrix(dethis))),1,function(k,dethis,thr,graph){
		IsingStateProb(dethis[k,], graph, thr[k,], beta=1, responses = c(-1L, 1L))
	} ,as.matrix( dethis), as.matrix( thr), as.matrix( graph))
	
	return(sum(log(Pdet_site)))
	
}

Sample_Ising_det_single_site = function(thr, Z, dethis, sppmat_det,nIter,n=1, method = "CFTP"){
	spp_exist = Z==1
	dethis[,!spp_exist] = -1# convert it to nrow = nperiod, ncol = nspp for single site, single repeat
	if(sum(spp_exist)==0) return(dethis)
	graph = sppmat_det[spp_exist,spp_exist]
	thr_exis = as.matrix( thr[,spp_exist])
	thr_abs = - apply(as.matrix(sppmat_det[!spp_exist,spp_exist]),2,sum) # condition on some species not exist here thus never be detection
	thr = apply(matrix(1:ncol(thr_exis)),1,function(k,ww,kk){ww[,k]+kk[k]},thr_exis,( thr_abs))
	dethis_exist = dethis[,spp_exist]
	dethis_exist = apply(matrix(1:nrow( as.matrix( dethis))),1,function(k,dethis_exist,thr,graph,nIter,n,method){
		IsingSampler(n=n,graph = graph, thresholds = thr[k,], beta=1, responses = c(-1L, 1L),nIter = nIter)
	}, as.matrix( dethis), as.matrix( thr), as.matrix( graph),nIter,n,method)
	dethis[,spp_exist] = dethis_exist
	return(dethis)
}

extract_thr = function(i,thr_list){
	nspp = length(thr_list)
	thr = sapply(thr_list,function(thr,i){t(thr[i,])},i=i) # thr at site i for all spps, will return a matrix with ncol = nspp, nrow = nperiod
	return(thr)
}

Pdet_Ising = function(nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,no_obs){
	require(IsingSampler)
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
		thr = (matrix(unlist(temp),nrow = n_row,ncol = nperiod))
		return(thr) # now here is basically a matrix, for each species at site and period
		},detDesign,beta_det,npardet,nrow(envX),nperiod) # this is gonna be  a list for all species, 
	
	Pdet = lapply(1:nsite,function(i,thr_list,detmat,Z,sppmat_det,nsite,nspp){
		thr = extract_thr(i,thr_list)
		rows1 = i + (1:nspp-1)*nsite
		dethis = t(detmat[rows1,])
		Z_site = Z[rows1,]
		Pdet_Ising_single_site(thr, Z_site, dethis, sppmat_det)
	},thr_list,detmat,as.matrix( Z),sppmat_det,nsite,nspp)# loop over sites
	 Pdet[no_obs]=0
	return(Reduce('+',Pdet))
}

## sampleIsingdet
Sample_Ising_detection = function(nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,nIter=100,n=1, method = "CFTP"){
	require(IsingSampler)
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

Pdet_Ising_rep = function(nrep,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,no_obs=NULL){
  if(!is.null(no_obs)) no_obs = as.matrix(no_obs) 
  Pdets = lapply(1:nrep,function(k,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,no_obs,nIter,n, method){
    Pdet_Ising(nperiod,envX,detX[[k]],beta_det,sppmat_det,Z[,k],detmat[[k]],no_obs[,k])
  },nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,(no_obs),nIter,n, method)
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

IsingOccu_multi.logL.innorm = function(theta, envX, distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp", Z ,detmat, detX,no_obs){ # the in-normalized log likelihood of IsingOccu Model beta is matrix here detX should be a list of list detmat should be a list, they should have the same length
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
		temp = as.numeric((Z[,i]+1)/2) * LP_Z1 + as.numeric(1-((Z[,i]+1)/2)) * LP_Z0
		temp[no_obs[,i]] = 0 # not so good, assume no_obs is a matrix with 1 shows there is no camera/camera fail there (no observation)
		logLdata = logLdata +  sum(temp)
	}
	return(negPot+logLdata)
}
  # passed 2019/3/18

IsingOccu_Ising_det_multi_logL_innorm = function(theta, envX, distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp", Z ,detmat, detX,no_obs){ # the in-normalized log likelihood of IsingOccu Model beta is matrix here detX should be a list of list detmat should be a list, they should have the same length
  nspp = nrow(theta$spp_mat)
  nperiod = ncol(detmat[[1]])
  beta_det = theta$beta_det
  negPot = Hamiltonian(theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z)
  nrep = ncol(Z)
  logLdata = Pdet_Ising_rep(nrep,nperiod,envX,detX,theta$beta_det,theta$spp_mat_det,Z,detmat,no_obs)
  return(negPot+logLdata)
}




Moller.ratio = function(theta_curr ,theta_prop
						,Z_curr ,Z_prop
						,Z_temp_curr, Z_temp_prop
						,detmat,no_obs # give detmat all 0 if no observation there
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
	log_q_theta_Z_prop_detmat = IsingOccu_multi.logL.innorm(theta_prop, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_prop ,detmat = detmat, detX,no_obs)
	# theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
	log_H_theta_Z_temp_curr = Hamiltonian(theta_curr,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp_curr)

	#### end of the upper part, start the lower

	log_H_theta_tuta_Z_temp_curr = Hamiltonian(theta_tuta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp_curr)
	log_pi_theta_curr = lapply(theta_curr,function(theta_temp,vars_prior){ sum(log(dnorm(theta_temp,0,sd=sqrt(vars_prior))))},vars_prior)
	log_pi_theta_curr = sum(unlist(log_pi_theta_curr))
	log_q_theta_Z_curr_detmat = IsingOccu_multi.logL.innorm(theta_curr, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_curr ,detmat = detmat, detX,no_obs)
	log_H_theta_Z_temp_prop = Hamiltonian(theta_prop,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp_prop)

	log_MH_ratio = (log_H_theta_tuta_Z_temp_prop + log_pi_theta_prop + log_q_theta_Z_prop_detmat + log_H_theta_Z_temp_curr)-
				   (log_H_theta_tuta_Z_temp_curr + log_pi_theta_curr + log_q_theta_Z_curr_detmat + log_H_theta_Z_temp_prop)

	return(min(1,exp(log_MH_ratio)))
}
  # passed 2019/3/18

Murray.ratio = function(theta_curr ,theta_prop
						,Z_curr ,Z_prop
						,Z_temp 
						,detmat,no_obs
						,vars_prior
						,envX, detX
						,distM,link_map
						,dist_mainland,link_mainland
						,int_range_intra="nn",int_range_inter="exp"){
	#log_H_theta_tuta_Z_temp_prop = Hamiltonian(theta_tuta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp_prop)
	# then auxiliented variable x_prop is same to detmat, together with Z_temp_prop from underlaying Isingmodel. It was proposed using likelihood function with parameter theta_prop and in the main sampler, which is important in canceling out the normalizing constant.
	log_pi_theta_prop = lapply(theta_prop,function(theta_temp,vars_prior){ sum(log(dnorm(theta_temp,0,sd=sqrt(vars_prior))))},vars_prior)
	log_pi_theta_prop = sum(unlist(log_pi_theta_prop))
	#log_pi_theta_prop = sum(log_pi_theta_prop)
	#prior of proposed theta
	log_q_theta_Z_prop_detmat = IsingOccu_multi.logL.innorm(theta_prop, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_prop ,detmat = detmat, detX, no_obs)
	# theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
	log_H_theta_curr_Z_temp = Hamiltonian(theta_curr,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp)

	#### end of the upper part, start the lower

	#log_H_theta_tuta_Z_temp_curr = Hamiltonian(theta_tuta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp_curr)
	log_pi_theta_curr = lapply(theta_curr,function(theta_temp,vars_prior){ sum(log(dnorm(theta_temp,0,sd=sqrt(vars_prior))))},vars_prior)
	log_pi_theta_curr = sum(unlist(log_pi_theta_curr))
	log_q_theta_Z_curr_detmat = IsingOccu_multi.logL.innorm(theta_curr, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_curr ,detmat = detmat, detX, no_obs=NULL)
	log_H_theta_prop_Z_temp = Hamiltonian(theta_prop,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp)

	log_MH_ratio = (log_pi_theta_prop + log_q_theta_Z_prop_detmat + log_H_theta_curr_Z_temp)-
				   (log_pi_theta_curr + log_q_theta_Z_curr_detmat + log_H_theta_prop_Z_temp)

	return(min(1,exp(log_MH_ratio)))
}

Murray.ratio.Ising_det = function(theta_curr ,theta_prop
                        ,Z_curr ,Z_prop
                        ,Z_temp 
                        ,detmat,no_obs
                        ,vars_prior
                        ,envX, detX
                        ,distM,link_map
                        ,dist_mainland,link_mainland
                        ,int_range_intra="nn",int_range_inter="exp"){
  #log_H_theta_tuta_Z_temp_prop = Hamiltonian(theta_tuta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp_prop)
  # then auxiliented variable x_prop is same to detmat, together with Z_temp_prop from underlaying Isingmodel. It was proposed using likelihood function with parameter theta_prop and in the main sampler, which is important in canceling out the normalizing constant.
  log_pi_theta_prop = lapply(theta_prop,function(theta_temp,vars_prior){ sum(log(dnorm(theta_temp,0,sd=sqrt(vars_prior))))},vars_prior)
  log_pi_theta_prop = sum(unlist(log_pi_theta_prop))
  #log_pi_theta_prop = sum(log_pi_theta_prop)
  #prior of proposed theta
  log_q_theta_Z_prop_detmat = IsingOccu_Ising_det_multi_logL_innorm(theta_prop, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_prop ,detmat = detmat, detX, no_obs)
  # theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
  log_H_theta_curr_Z_temp = Hamiltonian(theta_curr,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp)
  
  #### end of the upper part, start the lower
  
  #log_H_theta_tuta_Z_temp_curr = Hamiltonian(theta_tuta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp_curr)
  log_pi_theta_curr = lapply(theta_curr,function(theta_temp,vars_prior){ sum(log(dnorm(theta_temp,0,sd=sqrt(vars_prior))))},vars_prior)
  log_pi_theta_curr = sum(unlist(log_pi_theta_curr))
  log_q_theta_Z_curr_detmat = IsingOccu_Ising_det_multi_logL_innorm(theta_curr, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_curr ,detmat = detmat, detX, no_obs)
  log_H_theta_prop_Z_temp = Hamiltonian(theta_prop,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp)
  
  log_MH_ratio = (log_pi_theta_prop + log_q_theta_Z_prop_detmat + log_H_theta_curr_Z_temp)-
    (log_pi_theta_curr + log_q_theta_Z_curr_detmat + log_H_theta_prop_Z_temp)
  
  return(min(1,exp(log_MH_ratio)))
}
