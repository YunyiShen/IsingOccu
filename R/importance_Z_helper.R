# sampling Z from detection 
Hamdet_Ising_single_site = function(thr, Z, dethis, sppmat_det){
	spp_exist = Z==1
    sppmat_det = as(sppmat_det,"dgCMatrix")
	if(sum(spp_exist)==0 | sum(!is.na(dethis))==0){return(0)} # no species there, probability one to be no detection, or no observation here
	if(prod(spp_exist)==0){
	  thr_exis = as.matrix( thr[,spp_exist])
	  # do not include thr_abs since if two species never coexist we cannot infer "what if they coexist", i.e. thr_exis will be total colinear with thr_exis
	  thr = thr_exis
	  #thr = apply(matrix(1:ncol(thr_exis)),1,function(k,ww,kk){ww[,k]+kk[k]},thr_exis,( thr_abs))
	}
	graph = sppmat_det[spp_exist,spp_exist]
	has_obs = !(is.na(rowSums(dethis)))
	dethis = dethis[has_obs,spp_exist]# convert it to nrow = nperiod, ncol = nspp for single site, single repeat
	thr = thr[has_obs,]
	
	Hamdet_site = apply(matrix(1:sum(has_obs)),1,function(k,dethis,thr,graph){
		H(graph,(dethis[k,]), matrix( thr[k,]))
	} ,matrix( dethis,sum(has_obs),sum(spp_exist)), matrix( thr,sum(has_obs),sum(spp_exist)), as( as.matrix(graph),'dsCMatrix'))
	
	return(sum((Hamdet_site)))
	
}



Hamdet_Ising = function(nperiod,envX,detX,beta_det,sppmat_det,Z,detmat){
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
	
	Hamdet = lapply(1:nsite,function(i,thr_list,detmat,Z,sppmat_det,nsite,nspp){
		thr1 = extract_thr(i,thr_list)
		rows1 = i + (1:nspp-1)*nsite
		dethis = t(detmat[rows1,])
		Z_site = Z[rows1,]
		Hamdet_Ising_single_site(thr1, Z_site, dethis, sppmat_det)
	},thr_list,detmat,as.matrix( Z),sppmat_det,nsite,nspp)# loop over sites
	return(Reduce(rbind,Hamdet)) # change 25/8/2019
}


Hamdet_Ising_rep = function(nrep,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat){
  Hamdets = lapply(1:nrep,function(k,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat){
    sum(Hamdet_Ising(nperiod,envX,detX[[k]],beta_det,sppmat_det,Z[,k],detmat[[k]]))
  },nperiod,envX,detX,beta_det,sppmat_det,Z,detmat)
  return(Reduce('+',Hamdets))
}


propose_Z = function(theta, envX, detX,detmat,Z_curr,Z_absolute,Zprop_rate){
  nperiod = ncol(detmat)
  sppmat_det = theta$spp_mat_det
  Pdets = Pdet_Ising(nperiod,envX,detX,theta$beta_det,sppmat_det,Z_curr,detmat)
  Pdets_flip = Pdet_Ising(nperiod,envX,detX,theta$beta_det,sppmat_det,-Z_curr,detmat) # if flip will increase the likelihood of detection, flip
  rs = runif(length(Pdets)*ncol(sppmat_det))
  flip_or_not = ( log(rs)-log(Zprop_rate)<=rep(Pdets_flip-Pdets,ncol(sppmat_det)))
  Z_prop = Z_curr
  which_flip = which(Z_absolute==-1 & flip_or_not) 
  Z_prop[which_flip] = -Z_prop[which_flip]
  
  return(Z_prop)
}

propose_Z_rep = function(theta, envX, detX,detmat,Z_curr,Z_absolute,Zprop_rate){
  temp = lapply(1:ncol(Z_curr), function(i,theta, envX, detX,detmat,Z_curr,Z_absolute,Zprop_rate){
    propose_Z(theta, envX, detX[[i]],detmat[[i]],Z_curr[,i],Z_absolute[,i],Zprop_rate) 
  },theta, envX, detX,detmat,Z_curr,Z_absolute,Zprop_rate)
  
  return(as.matrix(Reduce(cbind,temp))) # matrix
}       




propose_rate = function(theta, envX, detX,detmat,Z_prop,Z_curr,Z_absolute,Zprop_rate){
  sppmat_det = theta$spp_mat_det
  nperiod = ncol(detmat)
  nspp = ncol(sppmat_det)
  Z_prop_temp = matrix(Z_prop,ncol = nspp)
  Z_curr_temp = matrix(Z_curr,ncol = nspp)
  Z_abso_temp = matrix(Z_absolute,ncol = nspp)
  
  
  flip_site = Z_prop_temp!=Z_curr_temp # find flipped sites
  
  No_flip = Z_abso_temp==-1 & !flip_site
  
  rm(Z_prop_temp,Z_curr_temp,Z_abso_temp)
  
  P_prop_Z =Zprop_rate * sapply(exp( Hamdet_Ising(nperiod,envX,detX,theta$beta_det,sppmat_det,-Z_prop,detmat) - 
								  Hamdet_Ising(nperiod,envX,detX,theta$beta_det,sppmat_det,Z_prop,detmat)),min,1) # probability of flipping
  P_curr_Z =Zprop_rate * sapply(exp( Hamdet_Ising(nperiod,envX,detX,theta$beta_det,sppmat_det,-Z_curr,detmat) - 
								  Hamdet_Ising(nperiod,envX,detX,theta$beta_det,sppmat_det,Z_curr,detmat)),min,1)
  
  prop2curr = sum(apply(flip_site,2,function(flip,P){sum(log(P[flip]))},P = P_prop_Z)) +  # fliped site, from prop to curr
              sum(apply(No_flip,2,function(No_flip,P){sum(log(1-P[No_flip]))},P = P_prop_Z)) # sites that can flip but did not 
  
  curr2prop = sum(apply(flip_site,2,function(flip,P){sum(log(P[flip]))},P = P_curr_Z)) + # fliped site, from curr to prop
              sum(apply(No_flip,2,function(No_flip,P){sum(log(1-P[No_flip]))},P = P_curr_Z)) 
  

    
  
  return(list(
  prop2curr = prop2curr,
  curr2prop = curr2prop
    )
  )
}

propose_rate_rep = function(theta, envX, detX,detmat,Z_prop,Z_curr,Z_absolute,Zprop_rate){
  nrep = ncol(Z_curr)
  temp = lapply(1:nrep,function(i,theta, envX, detX,detmat,Z_prop,Z_curr,Z_absolute,Zprop_rate){
    propose_rate(theta, envX, detX[[i]],detmat[[i]],Z_prop[,i],Z_curr[,i],Z_absolute[,i],Zprop_rate)
  },theta, envX, detX,detmat,Z_prop,Z_curr,Z_absolute,Zprop_rate)
  
  prop2curr = lapply(temp,function(t){t$prop2curr})
  curr2prop = lapply(temp,function(t){t$curr2prop})
  
  return(
    list(
      prop2curr = Reduce('+',prop2curr),
      curr2prop = Reduce('+',curr2prop)
    )
  )
  
}

MH_ratio_Z = function(theta, Z_curr, Z_prop,Z_absolute,Zprop_rate
                      ,detmat,envX, detX
                      ,distM,link_map
                      ,dist_mainland,link_mainland
                      ,int_range_intra="nn",int_range_inter="exp"){
  propose_rate_Z = propose_rate_rep(theta, envX, detX,detmat,Z_prop,Z_curr,Z_absolute,Zprop_rate)
  # numerator
  log_q_theta_Z_prop_detmat = IsingOccu_Ising_det_multi_logL_innorm(theta, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_prop ,detmat = detmat, detX)
  log_p_theta_Z_curr = propose_rate_Z$prop2curr

  # denominator
  log_q_theta_Z_curr_detmat = IsingOccu_Ising_det_multi_logL_innorm(theta, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_curr ,detmat = detmat, detX)
  log_p_theta_Z_prop = propose_rate_Z$curr2prop
  MH_ratio = log_q_theta_Z_prop_detmat + log_p_theta_Z_curr - log_q_theta_Z_curr_detmat - log_p_theta_Z_prop
  return(min(1,exp(MH_ratio)))
}       


propose_Z_plain = function(Z_curr,Z_absolute,Zprop_rate){
  Z_prop = Z_curr
  if(runif(1)<Zprop_rate){
    absence = which(Z_absolute==-1)
    #abs_obs = absence[!absence%in%no_obs]
    flip = sample(absence,1)
    Z_prop[flip]=-Z_prop[flip]
  }
  return(Z_prop)
}       



MH_ratio_Z_plain = function(theta, Z_curr, Z_prop,Z_absolute,Zprop_rate
                      ,detmat,envX, detX
                      ,distM,link_map
                      ,dist_mainland,link_mainland
                      ,int_range_intra="nn",int_range_inter="exp"){
  # numerator
  log_q_theta_Z_prop_detmat = IsingOccu_Ising_det_multi_logL_innorm(theta, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_prop ,detmat = detmat, detX)

  # denominator
  log_q_theta_Z_curr_detmat = IsingOccu_Ising_det_multi_logL_innorm(theta, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_curr ,detmat = detmat, detX)
  MH_ratio = log_q_theta_Z_prop_detmat - log_q_theta_Z_curr_detmat 
  return(min(1,exp(MH_ratio)))
}

