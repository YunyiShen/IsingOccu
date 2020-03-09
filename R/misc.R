## misc helper functions
getintralayerGraph = function(distM,link_map,eta,d,int_range = "exp",spp_mat) #it can be used multiple times for interislan and intra-island
{
  # pass all graphs as sparse matrix in package Matrix
  nspp = nrow(spp_mat) # which is the interspecific neighborhood matrix
  A = list() # intralayer graphs are passed using lists
  link_map = as(as.matrix(link_map),"dgCMatrix")
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
  return(A) # if link map is sparse, then A is sparse
} 
  # passed 2019/3/18

getfullGraph = function(A_ex,A_in,spp_mat){
  nspp = nrow(spp_mat)
  nsite = nrow(A_ex[[1]])
  A = Matrix(0,nspp*nsite,nspp*nsite,sparse = T)
  for(i in 2:nspp-1){
    A[1:nsite + (i-1)*nsite,1:nsite + (i-1)*nsite]=A_ex[[i]] + A_in[[i]] # diagonal part
    for(j in (i+1):nspp){
      
      diag(A[1:nsite + (i-1)*nsite,1:nsite + (j-1)*nsite])=spp_mat[i,j]
      diag(A[1:nsite + (j-1)*nsite,1:nsite + (i-1)*nsite])=spp_mat[j,i]
      
    }
  }
  i=nspp
  A[1:nsite + (i-1)*nsite,1:nsite + (i-1)*nsite]=A_ex[[i]] + A_in[[i]]
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
	Z = PartitionCpp(graph,thresholds,beta,responses)
    sapply(s, function(x) exp(-beta * H(graph, x, ( thresholds)))/Z)
}

Hamiltonian = function(MRF,Z_vec){
	nrep = ncol(Z_vec)
	Ham = lapply(1:nrep,function(i,Z,J,h){H(J,Z[,i],h)},Z=Z_vec,J=MRF$A,h=( MRF$thr))
  
	Ham = Reduce(rbind,Ham)
	return(Ham) # if we have repeat, just make Z_vec has two cols 
	
}


rIsingOccu_multi = function(MRF,n=1,method = "CFTP",nIter = 100){
	Z = IsingSamplerCpp(n=n,graph = MRF$A,thresholds=MRF$thr, responses = matrix( c(-1L, 1L),2,1),beta = 1,nIter=nIter,exact = (method=="CFTP"),constrain = NA + MRF$thr)
  return(t(Z))
	# test for 2spp case, passed 3/18/2019
}
  
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
	return(Reduce(rbind,Pdet)) # change 25/8/2019
}


Pdet_Ising_rep = function(nrep,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat){
  Pdets = lapply(1:nrep,function(k,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat){
    sum(Pdet_Ising(nperiod,envX,detX[[k]],beta_det,sppmat_det,Z[,k],detmat[[k]])) # change 25/8/2019
  },nperiod,envX,detX,beta_det,sppmat_det,Z,detmat)
  return(Reduce('+',Pdets))
}



IsingOccu_Ising_det_multi_logL_innorm = function(MRF, beta_det,sppmat_det, Z ,detmat, envX,detX){ # the in-normalized log likelihood of IsingOccu Model beta is matrix here detX should be a list of list detmat should be a list, they should have the same length
  nspp = nrow(sppmat_det)
  nperiod = ncol(detmat[[1]])
  negPot = Hamiltonian(MRF,Z)
  negPot = -sum(negPot)
  nrep = ncol(Z)
  logLdata = Pdet_Ising_rep(nrep,nperiod,envX,detX,beta_det,sppmat_det,Z,detmat)
  return(negPot+logLdata)
}


         
write_json.IsingOccu_samples = function(x,path){
  n_sample = nrow(x$Z.mcmc)
  x$theta.mcmc = lapply(x$theta.mcmc,matrix,nrow = n_sample)
  x$Z.mcmc = matrix(x$Z.mcmc,nrow = n_sample)
  class(x) = 'list'
  jsonlite::write_json(x,path)
}
		   
		   
adjacency.matrix = function(m, n = NULL)
{
    if (missing(n))
    {
        A = Matrix(0, m^2, m^2,sparse = T)
        for (i in 1:m^2)
        {
            up = i - m
            down = i + m
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= m^2)
                A[i, down] = 1
            if (left %% m != 0)
                A[i, left] = 1
            if (i %% m != 0)
                A[i, right] = 1
        }
    }
    else
    {
        A = Matrix(0, m * n, m * n,sparse = T)
        for (i in 1:(m * n))
        {
            up = i - n
            down = i + n
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= (m * n))
                A[i, down] = 1
            if (left %% n != 0)
                A[i, left] = 1
            if (i %% n != 0)
                A[i, right] = 1
        }
    }
    A
}   
		   
		   
		   
		   