propose_Z_plain <- function(Z_curr,Z_absolute,Zprop_rate){
  Z_prop <- Z_curr
  if(runif(1)<Zprop_rate){
    absence <- which(Z_absolute==-1)
    #abs_obs = absence[!absence%in%no_obs]
    flip <- sample(absence,1)
    Z_prop[flip] <- -Z_prop[flip]
  }
  return(Z_prop)
}       



MH_ratio_Z_plain <- function(theta, MRF,Z_curr, Z_prop
                      ,detmat,envX, detX
                      ){
  # numerator
  log_q_theta_Z_prop_detmat <- IsingOccu_Ising_det_multi_logL_innorm(MRF, theta$beta_det,theta$spp_mat_det, Z_prop ,detmat, envX,detX)

  # denominator
  log_q_theta_Z_curr_detmat <- IsingOccu_Ising_det_multi_logL_innorm(MRF, theta$beta_det,theta$spp_mat_det, Z_curr ,detmat, envX,detX)
  MH_ratio <- log_q_theta_Z_prop_detmat - log_q_theta_Z_curr_detmat 
  return(min(1,exp(MH_ratio)))
}



Gibbs_Z <- function(theta, envX, detX,detmat,Z_curr,Z_absolute,MRF){
  Z_curr <- as.matrix(Z_curr)
  nperiod <- ncol(detmat)
  sppmat_det <- theta$spp_mat_det
  nspp <- nrow(theta$spp_mat_det)
  n_site <- nrow(envX)
  minus1s <- which(Z_absolute==-1) # only scan 
  beta_det <- theta$beta_det
  A <- MRF$A
  thr <- MRF$thr
  if(is.null(detX)) {
      detDesign <- lapply(1:nperiod,function(dummy,envX){envX},envX) 
  }
  else detDesign <- lapply(detX,function(x,y){ as.matrix( cbind(y,x))},y = envX)
  npardet <- ncol(detDesign[[1]])
  nsite <- nrow(envX)
  nspp <- nrow(sppmat_det)
  
  det_thr_list <- lapply( 1:nspp, function(i,detDesign,beta_det,naprdet,n_row,nperiod){ 
        temp <- lapply(detDesign,function(w,beta1,i){w%*%beta1},beta1 = matrix( beta_det[1:npardet + (i-1) * npardet]),i=i)
        thr1 <- (matrix(unlist(temp),nrow = n_row,ncol = nperiod))
        return(thr1) # now here is basically a matrix, for each species at site and period
      },detDesign,beta_det,npardet,nrow(envX),nperiod) # this is gonna be  a list for all species, each element was a site, which is helpful
    


  sites_adder <- (1:nspp-1)*nsite
  
  for(scan in minus1s){
    Ham_plus1 <- A[scan,] %*% Z_curr + thr[scan] # hamiltonian when it is +1
    
    which_site <- scan %% nsite # determin which site to work with 
    if(which_site==0) which_site <- nsite
    which_spp <- (scan-which_site)/nsite + 1 # determin which species we were working with 
    which_row <- sites_adder + which_site # the row of all species at such site
    
    Z_temp <- Z_curr[which_row,] # grab all species at such site
    Z_temp[which_spp] <- 1
    
    dethis <- t(detmat[which_row,])
    
    det_thr_temp <- extract_thr(which_site,det_thr_list) # this get the thr of such site
    
    Pplus1 <- Pdet_Ising_single_site(det_thr_temp, Z_temp, dethis, sppmat_det) # this was logged
    Z_temp[which_spp] <- -1
    Pminus1 <- Pdet_Ising_single_site(det_thr_temp, Z_temp, dethis, sppmat_det)
    
    PZiplus1 <- (exp(Ham_plus1+Pplus1))/(exp(-Ham_plus1+Pminus1)+exp(Ham_plus1+Pplus1))
    
    Z_curr[scan,] <- ifelse(runif(1)<PZiplus1,1,-1)
    
  }

  
  return(Z_curr)
}





Gibbs_Z_rep <- function(theta, envX, detX,detmat,Z_curr,Z_absolute,MRF){
  temp <- lapply(1:ncol(Z_curr), function(i,theta, envX, detX,detmat,Z_curr,Z_absolute,MRF){
    Gibbs_Z(theta, envX, detX[[i]],detmat[[i]],Z_curr[,i],Z_absolute[,i],MRF) 
  },theta, envX, detX,detmat,Z_curr,Z_absolute,MRF)
  return(as.matrix(Reduce(cbind,temp))) # matrix
}       
























