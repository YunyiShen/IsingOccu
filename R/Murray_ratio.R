## various M-H ratios and Murray ratios used in sampler

getlogprior_normal = function(theta_prop,theta_curr,vars_prior){
  n_para = length(theta_prop)
  
  log_pi_theta_prop = lapply(1:n_para,function(i,theta_temp,vars_prior){ sum(log(dnorm(as.vector(theta_temp[[i]]),0,sd=sqrt(vars_prior[[i]]))))},theta_prop,vars_prior)
  log_pi_theta_prop = sum(unlist(log_pi_theta_prop))
  
  log_pi_theta_curr = lapply(1:n_para,function(i,theta_temp,vars_prior){ sum(log(dnorm(as.vector( theta_temp[[i]]),0,sd=sqrt(vars_prior[[i]]))))},theta_curr,vars_prior)
  log_pi_theta_curr = sum(unlist(log_pi_theta_curr))

  return(list(prop = log_pi_theta_prop,curr = log_pi_theta_curr))
}

getlogprior_uniform = function(theta_prop,theta_curr,lim_prior){
  log_pi_theta_prop = lapply(1:length(theta_prop),function(i,theta,lim_prior){ 
	  sum(log(dunif(as.vector(theta[[i]]),-as.vector(lim_prior[[i]]),as.vector(lim_prior[[i]]))))}
							 ,theta_prop,lim_prior)
  log_pi_theta_prop = sum(unlist(log_pi_theta_prop))
  
  log_pi_theta_curr = lapply(1:length(theta_curr),function(i,theta,lim_prior){ 
	  sum(log(dunif(as.vector(theta[[i]]),-as.vector(lim_prior[[i]]),as.vector(lim_prior[[i]]))))}
							 ,theta_curr,lim_prior)
  log_pi_theta_curr = sum(unlist(log_pi_theta_prop))
  return(list(prop = log_pi_theta_prop,curr = log_pi_theta_curr))
}

Murray_ratio_occu_theta = function( MRF_curr, MRF_prop, log_pi
                        ,Z
                        ,Z_temp
                        ,vars_prior
                        ,distM,link_map
                        ,dist_mainland,link_mainland
                        ,int_range_intra="nn",int_range_inter="exp"){
  
  #prior of proposed theta
  log_q_theta_prop_Z = -sum(Hamiltonian(MRF_prop,Z ))
  # theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
  log_H_theta_curr_Z_temp = -sum(Hamiltonian(MRF_curr,Z_temp))
  
  #### end of the numerator part, start the denominator
  
  log_q_theta_curr_Z = -sum(Hamiltonian(MRF_curr,Z ))
  log_H_theta_prop_Z_temp = -sum(Hamiltonian(MRF_prop,Z_temp))
  
  log_MH_ratio = (log_pi$prop + log_q_theta_prop_Z + log_H_theta_curr_Z_temp)-
    (log_pi$curr + log_q_theta_curr_Z + log_H_theta_prop_Z_temp)
  
  return(min(1,exp(log_MH_ratio)))
}

MH.ratio.Ising_det = function(theta_curr ,theta_prop,log_pi
                        ,Z
                        ,detmat
                        ,vars_prior
                        ,envX, detX
                        ){
  nrep = ncol(Z)
  nperiod = ncol(detmat[[1]])  
  #prior of proposed theta
  log_q_theta_prop_Z_detmat = Pdet_Ising_rep(nrep,nperiod,envX,detX,theta_prop$beta_det,theta_prop$spp_mat_det,Z,detmat)  
  #### end of the numerator part, start the denominator
  
  log_q_theta_curr_Z_detmat = Pdet_Ising_rep(nrep,nperiod,envX,detX,theta_curr$beta_det,theta_curr$spp_mat_det,Z,detmat)
  
  log_MH_ratio = (log_pi$prop + log_q_theta_prop_Z_detmat)-
    (log_pi$curr + log_q_theta_curr_Z_detmat )
  
  return(min(1,exp(log_MH_ratio)))
}



Murray.ratio.Ising_occudet = function(MRF_curr,MRF_prop,log_pi
                        ,Z_curr ,Z_prop
                        ,Z_temp 
                        ,detmat
                        ,beta_det_curr
					    ,beta_det_prop
                        ,detX
                        ,sppmat_det_curr, sppmat_det_prop){
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
