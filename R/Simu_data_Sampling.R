Sample_Ising_detection_rep = function(nrep,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,nIter=100,n=1, method = "CFTP"){
  detmat = lapply(1:nrep,function(k,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,nIter,n, method){
    Sample_Ising_detection(nperiod,envX,detX[[k]],beta_det,sppmat_det,Z,detmat[[k]],nIter,n, method)
  },nperiod,envX,detX,beta_det,sppmat_det,Z,detmat,nIter,n, method)
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
