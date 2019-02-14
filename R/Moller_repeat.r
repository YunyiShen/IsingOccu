Moller.ratio = function(theta_curr ,theta_prop
                        ,Z_curr ,Z_prop
                        ,Z_temp_curr, Z_temp_prop
                        ,vars_prior
                        ,theta_tuta
                        ,envX, distM ){
  log_H_theta_tuta_Z_temp_prop = Hamiltonian(theta_tuta, envX, distM, Z_temp_prop )
  # then auxiliented variable x_prop is same to detmat, together with Z_temp_prop from underlaying Isingmodel. It was proposed using likelihood function with parameter theta_prop and in the main sampler, which is important in canceling out the normalizing constant.
  log_pi_theta_prop =log(dnorm(theta_prop,0,sd=sqrt(vars_prior)))
  log_pi_theta_prop = sum(log_pi_theta_prop)
  #prior of proposed theta
  log_q_theta_Z_prop_detmat = Hamiltonian(theta_prop, envX, distM, Z_prop )
  # theta_prop should be sample from independent Gaussian distribution with mean theta_curr, Z_prop should be directly sample from a uniform configuration (of course where exist detection should be 1 with probability 1, actually sample all 0s, then we can cancel out the proposal probability from the MH ratio)
  log_H_theta_Z_temp_curr = Hamiltonian(theta_curr, envX, distM, Z_temp_curr )
  
  #### end of the upper part, start the lower
  
  log_H_theta_tuta_Z_temp_curr = Hamiltonian(theta_tuta, envX, distM, Z_temp_curr )
  log_pi_theta_curr =log(dnorm(theta_curr,0,sd=sqrt(vars_prior)))
  log_pi_theta_curr = sum(log_pi_theta_curr)
  log_q_theta_Z_curr_detmat = Hamiltonian(theta_curr, envX, distM, Z_curr )
  log_H_theta_Z_temp_prop = Hamiltonian(theta_prop, envX, distM, Z_temp_prop )
  
  log_MH_ratio = (log_H_theta_tuta_Z_temp_prop + log_pi_theta_prop + log_q_theta_Z_prop_detmat + log_H_theta_Z_temp_curr)-
    (log_H_theta_tuta_Z_temp_curr + log_pi_theta_curr + log_q_theta_Z_curr_detmat + log_H_theta_Z_temp_prop)
  
  return(min(1,exp(log_MH_ratio)))
}


Hamiltonian = function(theta, envX, distM, Z){
  Z=matrix(Z,nrow=length(Z),ncol=1) # make Z col vector
  p = length(theta)
  nsite = nrow(distM)
  ncov = ncol(envX)
  
  # zeros = matrix(0,nrow=nsite,ncol=ncov)
  d = theta[p]
  beta1 = theta[-c(p,p-1)]
  eta = theta[p-1]
  # Xfull = cbind(rbind(envX,zeros),rbind(zeros,envX))
  thr1 = envX%*%beta1
  G = eta*(exp(-abs(d)*distM))
  diag(G)=0
  # rm(Xfull)
  #A = getGraph(distM,theta,int_range = int_range,full=FALSE)
  negPot = t(thr1) %*% Z  + .5* eta*t(Z)%*%G%*%Z
  
  return(negPot)
}


rIsing=function(X,distM,theta,method = "CFTP",nIter = 100,n=1){
  require(IsingSampler)
  p = length(theta)
  beta1 = theta[-c(p,p-1)]
  d=theta[p]
  eta = theta[p-1]
  thr = X%*%beta1
  G = eta*(exp(-abs(d)*distM))
  diag(G)=0
  Z = IsingSampler(n=n,graph = G,thresholds = thr,nIter = nIter,responses = c(-1L,1L),method = method)
  return(Z)
}

# THIS is the IsingOccu fitting function using Moller et al. 2006 sampler (if we can only use MCEM to do MPLE, then Bayesian is much faster)
# remember, X contains 1 col while detX doesn't because the design matrix of det is actually cbind(X,detX)
# detX should be a list, with every element is the design matrix WITHOUT 1s.
# bug here, never accept??, may need to try another way of propose change...
Moller.sampler_repeat = function(X,distM, Z,mcmc.save = 10000, burn.in = 10 , vars_prior = rep(1,4*ncol(X)+2*ncol(detX[[1]])+9),vars_prop = 2,int_range = "exp",seed = 12345,init,thin.by = 1){
  require(coda)
  require(IsingSampler)
  #source("misc.R")
  cat("Initializing...\n\n")
  set.seed(seed)
  #nsite = nrow(detmat)/2
  #abs_pres = ((rowSums(detmat)>0)) # absolute present
  #datatemp = data.frame(r = rowSums(detmat)>0,rbind(X,X))
  theta_curr = init
  
  #ncov = ncol(X)
  #ncov_det = 2 * ncol(detX[[1]]) + ncov
  
  
  theta_tuta = theta_curr # improvement possible here, make theta_tuta a better approximation of theta_post
  theta.mcmc = mcmc(matrix(nrow = (mcmc.save),ncol = length(theta_curr)),thin = thin.by)
  p = length(theta_tuta)
  #colnames(theta.mcmc)[(p - 4):p] = c("eta_spatial_spc1","d_spatial_spc1","eta_spatial_spc2","d_spatial_spc2","eta_interspecies") # eta spatial
  
  Z_temp_curr = rIsing(X,distM,theta = theta_curr,method = "CFTP",nIter = 100)
  Z_curr = Z
  # try do Z and theta separtly!
  cat("Burn in...\n")
  accept_Z = 0
  accept_theta_occu = 0
  accept_theta_det = 0
  low_acc_Z = 0
  low_acc_theta_occu = 0
  low_acc_theta_det = 0
  Z_curr = Z
  timing = proc.time()
  for(i in 1:burn.in){# to burn in
    #propose theta
    #occu theta first
    walk = rnorm(length(theta_curr),0,sd = sqrt(vars_prop))
    theta_prop = theta_curr + walk
    #theta_prop[c(1:ncov,1:5+(ncov+ncov_det))] = walk [c(1:ncov,1:5+(ncov+ncov_det))] + theta_prop[c(1:ncov,1:5+(ncov+ncov_det))]
    
    # Aux Z proposed from theta_prop and underlaying graph model to cancel out the normalizing constant
    Z_temp_prop = rIsing(X,distM,theta = theta_prop,method = "CFTP",nIter = 100)
    # MH ratio
    Moller_ratio = Moller.ratio(theta_curr ,theta_prop
                                ,Z_curr ,Z_curr
                                ,Z_temp_curr, Z_temp_prop						
                                ,vars_prior
                                ,theta_tuta # this is important, which control the importance sampling 
                                ,envX=X, distM=distM)
    if(Moller_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
    r = runif(1)
    if(r<=Moller_ratio){
      theta_curr=theta_prop
      Z_temp_curr = Z_temp_prop
      accept_theta_occu = accept_theta_occu + 1
    }
    
    # propose detection theta 
    
    
    if(i%%100 == 0) {
      
      cat("Burn in iteration: #",i,"\n")
      cat("# of occupancy theta acceptance: " , accept_theta_occu,"\n")
      cat("# of occupancy acceptance ratio <exp(-10): ",low_acc_theta_occu,"\n\n")
      if(accept_theta_occu==0) cat(theta_curr[c(1:ncov,1:5+(ncov+ncov_det))],"\n\n")
      timing = proc.time()- timing
      cat("Time used in this 100:",timing[1],"s\n")
      cat("\n\n")
      timing = proc.time()
      accept_Z = 0
      low_acc_Z = 0
      accept_theta_occu = 0
      low_acc_theta_occu = 0
      accept_theta_det = 0
      low_acc_theta_det = 0
    }
  }
  cat("Start sampling...\n")
  accept_Z = 0
  low_acc_Z = 0
  accept_theta_occu = 0
  low_acc_theta_occu = 0
  accept_theta_det = 0
  low_acc_theta_det = 0
  timing = proc.time()
  for(i in 1:(mcmc.save)){ # to save
    #propose theta
    # theta_occu first
    walk = rnorm(length(theta_curr),0,sd = sqrt(vars_prop))
    theta_prop = theta_curr+walk
    #theta_prop[c(1:ncov,1:5+(ncov+ncov_det))] = walk [c(1:ncov,1:5+(ncov+ncov_det))] + theta_prop[c(1:ncov,1:5+(ncov+ncov_det))]
    
    # Aux Z proposed from the proposal theta and underlaying graph model
    Z_temp_prop = rIsing(X,distM,theta = theta_prop,method = "CFTP",nIter = 100)
    # MH ratio
    Moller_ratio = Moller.ratio(theta_curr ,theta_prop
                                ,Z_curr ,Z_curr
                                ,Z_temp_curr, Z_temp_prop						
                                
                                ,vars_prior
                                ,theta_tuta
                                ,envX=X, distM=distM)
    r = runif(1)
    if(r<=Moller_ratio){
      theta_curr=theta_prop
      Z_temp_curr = Z_temp_prop
      accept_theta_occu = accept_theta_occu + 1
    }
    if(Moller_ratio<exp(-10)) low_acc_theta_occu = low_acc_theta_occu + 1
    
    # propose Z by single flip, may try Gibbs, this deals with missing obervation of Z
    theta.mcmc[i,]=theta_curr
    if(i%%100 == 0) { # reporting
      cat("Sampling in iteration: #",i,"\n")
      cat("# of occupancy theta acceptance: " , accept_theta_occu,"\n")
      cat("# of occupancy acceptance ratio <exp(-10): ",low_acc_theta_occu,"\n\n")
     # if(accept_theta_occu==0) cat(theta_curr[c(1:ncov,1:5+(ncov+ncov_det))],"\n\n")
      
      timing = proc.time()-timing
      cat("Time used in this 100:",timing[1],"s\n")
      cat("\n\n")
      timing = proc.time()
      accept_Z = 0
      low_acc_Z = 0
      accept_theta_occu = 0
      low_acc_theta_occu = 0
      accept_theta_det = 0
      low_acc_theta_det = 0
    }
  }
  res = list(theta.mcmc = theta.mcmc
             ,theta.mean = apply(theta.mcmc,2,mean)
             ,vars_prior=vars_prior
             )
  class(res)="Moller.test"
  return(res)
}




