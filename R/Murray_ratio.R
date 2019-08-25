## various M-H ratios and Murray ratios used in sampler
Murray_ratio_occu_theta = function(theta_curr ,theta_prop
                        ,Z
                        ,Z_temp
                        ,vars_prior
                        ,distM,link_map
                        ,dist_mainland,link_mainland
                        ,int_range_intra="nn",int_range_inter="exp"){
  log_pi_theta_prop = lapply(theta_prop,function(theta_temp,vars_prior){ sum(log(dnorm(as.vector(theta_temp),0,sd=sqrt(vars_prior))))},vars_prior)
  log_pi_theta_prop = sum(unlist(log_pi_theta_prop))
  
  #prior of proposed theta
  log_q_theta_prop_Z = -sum(Hamiltonian(theta_prop, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z ))
  # theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
  log_H_theta_curr_Z_temp = -sum(Hamiltonian(theta_curr,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp))
  
  #### end of the numerator part, start the denominator
  
  log_pi_theta_curr = lapply(theta_curr,function(theta_temp,vars_prior){ sum(log(dnorm(as.vector( theta_temp),0,sd=sqrt(vars_prior))))},vars_prior)
  log_pi_theta_curr = sum(unlist(log_pi_theta_curr))
  log_q_theta_curr_Z = -sum(Hamiltonian(theta_curr, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z ))
  log_H_theta_prop_Z_temp = -sum(Hamiltonian(theta_prop,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp))
  
  log_MH_ratio = (log_pi_theta_prop + log_q_theta_prop_Z + log_H_theta_curr_Z_temp)-
    (log_pi_theta_curr + log_q_theta_curr_Z + log_H_theta_prop_Z_temp)
  
  return(min(1,exp(log_MH_ratio)))
}

MH.ratio.Ising_det = function(theta_curr ,theta_prop
                        ,Z
                        ,detmat
                        ,vars_prior
                        ,envX, detX
                        ){
  nrep = ncol(Z)
  nperiod = ncol(detmat[[1]])
  log_pi_theta_prop = lapply(theta_prop,function(theta_temp,vars_prior){ sum(log(dnorm(as.vector(theta_temp),0,sd=sqrt(vars_prior))))},vars_prior)
  log_pi_theta_prop = sum(unlist(log_pi_theta_prop))
  
  #prior of proposed theta
  log_q_theta_prop_Z_detmat = Pdet_Ising_rep(nrep,nperiod,envX,detX,theta_prop$beta_det,theta_prop$spp_mat_det,Z,detmat)  
  #### end of the numerator part, start the denominator
  
  log_pi_theta_curr = lapply(theta_curr,function(theta_temp,vars_prior){ sum(log(dnorm(as.vector( theta_temp),0,sd=sqrt(vars_prior))))},vars_prior)
  log_pi_theta_curr = sum(unlist(log_pi_theta_curr))
  log_q_theta_curr_Z_detmat = Pdet_Ising_rep(nrep,nperiod,envX,detX,theta_curr$beta_det,theta_curr$spp_mat_det,Z,detmat)
  
  log_MH_ratio = (log_pi_theta_prop + log_q_theta_prop_Z_detmat)-
    (log_pi_theta_curr + log_q_theta_curr_Z_detmat )
  
  return(min(1,exp(log_MH_ratio)))
}



Murray.ratio.Ising_occudet = function(theta_curr ,theta_prop
                        ,Z
                        ,Z_temp
                        ,detmat
                        ,vars_prior
                        ,envX, detX
                        ,distM,link_map
                        ,dist_mainland,link_mainland
                        ,int_range_intra="nn",int_range_inter="exp"){
  log_pi_theta_prop = lapply(theta_prop,function(theta_temp,vars_prior){ sum(log(dnorm(as.vector(theta_temp),0,sd=sqrt(vars_prior))))},vars_prior)
  log_pi_theta_prop = sum(unlist(log_pi_theta_prop))
  
  #prior of proposed theta
  log_q_theta_prop_Z_detmat = IsingOccu_Ising_det_multi_logL_innorm(theta_prop, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z ,detmat = detmat, detX)
  # theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
  log_H_theta_curr_Z_temp = -sum(Hamiltonian(theta_curr,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp))
  
  #### end of the numerator part, start the denominator
  
  log_pi_theta_curr = lapply(theta_curr,function(theta_temp,vars_prior){ sum(log(dnorm(as.vector( theta_temp),0,sd=sqrt(vars_prior))))},vars_prior)
  log_pi_theta_curr = sum(unlist(log_pi_theta_curr))
  log_q_theta_curr_Z_detmat = IsingOccu_Ising_det_multi_logL_innorm(theta_curr, envX, distM, link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z ,detmat = detmat, detX)
  log_H_theta_prop_Z_temp = -sum(Hamiltonian(theta_prop,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_temp))
  
  log_MH_ratio = (log_pi_theta_prop + log_q_theta_prop_Z_detmat + log_H_theta_curr_Z_temp)-
    (log_pi_theta_curr + log_q_theta_curr_Z_detmat + log_H_theta_prop_Z_temp)
  
  return(min(1,exp(log_MH_ratio)))
}


