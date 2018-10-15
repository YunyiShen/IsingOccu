# this is a sampler using the model, only for the underlaying markov random field, which was used in bootstrap BUT NOT MCEM!!.
rautoccu = function(X, A, A1, A2,theta)
{	
	require(ngspatial)
	p = length(theta)
	A = theta[p - 2] * A + theta[p - 1] * A1 + theta[p - 0] * A2
	theta = as.numeric( matrix(c(theta[1:ncol(X)],1)))
	as.numeric(rautologistic(X, A, theta)) # use ngspatial
}

# sample detection given Z
autooccu_sample.detection = function(theta, X, A, A1, A2, Z, detmat, detX){
	p = length(theta)
  RN = matrix( runif(nrow(detmat)*ncol(detmat)),nrow = nrow(detmat),ncol = ncol(detmat) )
	nperiod = ncol(detmat) # detmat is the data of 0 and 1 for detections
	beta_det = theta[(ncol(X) + 1):(p-3)] # length(beta_det) = ncol(detX[[1]]) + ncol(X) + 3 # beta for detections
	detDesign = lapply(detX,function(x,y){cbind(y,x)},y = X) # This is the full design matrix list of detection probability p at time
	Xbeta_det = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = beta_det) # Xbeta for detections
	P_det = lapply(Xbeta_det,function(W){exp(W) / (1 + exp(W))}) # just a logistic regression for detection
	P_det = (matrix(unlist(P_det),nrow = nrow(X),ncol = nperiod)) # detection probability, row is site i col is period j
	
	occustatus = matrix(rep(Z,nperiod),nrow = length(Z),ncol=nperiod)
	P_det_occu = P_det * occustatus
	
	sample_detHistory = 1 * (RN<P_det_occu)
	sample_detHistory
}

#Sampler for MCEM, sample Z from the total likelihood given all detections and all theta
autooccuMCEM_sampleZ = function(theta, X, A, A1, A2, Z0, detmat, detX, num_sample){
	sample_Z = matrix(,nrow = length(Z0),ncol = num_sample)
	Z_naive = as.matrix(rowSums(detmat)>0)
	Zold = Z0
	#detold = detmat
	for(i in 1:num_sample){ #sample start
		#propose a new Z
		flip_site = sample(which(Z_naive==0),1)
		Znew = Zold
		Znew[flip_site] = 1 - Znew[flip_site] # one spin flip MC
		
		#calculate the Metroplis ratio
		Mratio = -autooccu.logPL(theta, X, A, A1, A2, Znew, detmat, detX) + autooccu.logPL(theta, X, A, A1, A2, Zold, detmat, detX)
		Mratio = min( c(exp(Mratio),1) )
		RN = runif(1)
		if(RN<Mratio){
			sample_Z[,i] = Znew
		}
		else{
			sample_Z[,i] = Zold
		}
		Zold = Znew
	}
	sample_Z
}



# PSEUDO LIKELIHOOD of AUTOOCCU MODEL
# Pseudo Likelihood here, first version with non perfect detection is done Sept/22/2018 and checked NOT CHECKED YET
autooccu.logPL = function(theta, envX, A, A1, A2, Z, detmat, detX, centered) # assume known Z, data should form as nrow(data)==nrow(X),ncol(data)=# of periods. detX is design matrix of detections, WITHOUT 1s, should be a list. length(detX) = ncol(data) = # of periods # A1 is dist matrix for species 2, all about first spc1 should be 0

# add another two parameters, s.t. A = exp(-exp(k)*distMat)
{
	# deal with the underlaying Ising model
	# REMEMBWE TO ADD 1s to envX, detX will automaticlly include because envX in cbinded
	# nsite = nrow(X)/2
  
	p = length(theta) #p should = 2 * ncol(X) + ncol(detX[[1]]) + 3
	beta_full = theta[-c(p-2,p-1,p)]
	eta = theta[p - 2] # spatial relation for spc1
	eta1 = theta[p - 1] # spatial relation for spc2
	eta2 = theta[p - 0] # interspecies relationship constant 
	beta = beta_full[1:ncol(envX)]  
	Xbeta = envX %*% beta # X should be doubled for 2 species cases, since environment is the same to both of them, maybe better to write X = rbind(X,X)
	mu = exp(Xbeta) / (1 + exp(Xbeta))
	#logPL = Xbeta + eta * A %*% (Z - mu) 
	logPL = Xbeta + (eta * A + eta1 * A1 + eta2 * A2) %*% (Z - centered * mu) #  ETA = eta * A + eta1 * A1 + eta2 * A2, in which eta1 is the spatial auto correlation constant for spc.2, eta is the spatial constant of spc.1 eta2 is the interspecies interaction constant
	logPL = t(Z) %*% logPL - sum(log(1 + exp(logPL))) # what need to do is just add detection history likelihoods here, while, Z should be EM.

	# now deal with non-perfect detection
	nperiod = ncol(detmat) # detmat is the data of 0 and 1 for detections
	beta_det = theta[(ncol(envX) + 1):(p-3)] # length(beta_det) = ncol(detX[[1]]) + ncol(X) + 3 # beta for detections
	detDesign = lapply(detX,function(x,y){cbind(y,x)},y = envX) # This is the full design matrix list of detection probability p at time
	Xbeta_det = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = beta_det) # Xbeta for detections
	rm(detDesign)
	P_det = lapply(Xbeta_det,function(W){exp(W) / (1 + exp(W))}) # just a logistic regression for detection
	rm(Xbeta_det)
	P_det = (matrix(unlist(P_det),nrow = nrow(envX),ncol = nperiod)) # log detection probability, row is site i col is period j
	LP_Z1 = as.matrix(rowSums(detmat * log(P_det) + (1-detmat) * log(1-P_det)))
	LP_Z0 = as.matrix(log(1*(rowSums(detmat)==0) + 1e-13 * (1-(rowSums(detmat)==0)))) # I(data = 0), do not want err for those have detections
	logLdata = sum(Z * LP_Z1 + (1-Z) * LP_Z0) # likelihood of data 
	
	# total neg log likelihood
	-logPL-logLdata 
}

# E-STEP of EM Algorithm
autooccu.ElogPL = function(theta, X, A, A1 , A2, Z, detmat, detX) # E step for EM
{
  require('snow')
  cl = makeSOCKcluster(detectCores())
  n_sample = ncol(Z)
  clusterExport(cl,list("Z","theta","X","A","A1","A2","detmat","detX","autooccu.logPL"),envir = environment())
  PLs = parCapply(cl=cl,x=Z,FUN='autooccu.logPL',theta = theta,
                            envX = X, A = A, A1 = A1, A2 = A2,  detmat = detmat, detX = detX)
  #w = as.array(Z)
	sumPL = sum(PLs)#calculate the E of likelihood
	EminuslogPL = sumPL/n_sample
	#print(EminuslogPL)
	stopCluster(cl)
	EminuslogPL
}

autooccu.ElogPLnotpr = function(theta, X, A, A1 , A2, Z, detmat, detX, centered) # E step for EM
{
  n_sample = ncol(Z)
  PLs = apply(X=Z,MARGIN = 2,FUN='autooccu.logPL',theta = theta,
                  envX = X, A = A, A1 = A1, A2 = A2,  detmat = detmat, detX = detX, centered)
  sumPL = sum(PLs)#calculate the E of likelihood
  EminuslogPL = sumPL/n_sample
  EminuslogPL
}


# AUTOOCCU FIT HERE, MCEM is still under developing 2018/9/22
autooccu.fit = function(X, A, A1, A2, detmat, detX, centered = T ,MCEMsample = 10000 ,hessian = T, method = 'BFGS')
{ 
  #require('snow')
  start = glm(as.numeric(rowSums(detmat)>0) ~ X - 1, family = binomial)$coef
	start_det = glm(detmat[,1] ~ cbind(X,detX[[1]]) - 1,family = binomial)$coef
	theta_previous = c(start,start_det,1,1,1)
	theta_previous = as.numeric(theta_previous)
	rownames(theta_previous) = NULL
	cond = T
	toll = 1e-5 # tolerance for EM algorithm
	Z_naive = as.matrix(rowSums(detmat)>0)
	LL0 = autooccu.logPL(theta_previous, X, A, A1, A2, Z_naive, detmat, detX, centered)
	print(paste("model started with -log(L) = ",as.character(LL0)))
	print("Initial MCEM")
	ww = 0
	while(cond){
	  ww = ww + 1
	  print("sampling Z...")
		Z = autooccuMCEM_sampleZ(theta, X, A, A1, A2, Z0 = Z_naive, detmat, detX, num_sample = MCEMsample) # USED SIMPLEST IMPLIMENT here, need to change to perfect sampler
		print("all clear, initial M step...")
		#opt = try(optim(theta_previous, autooccu.ElogPL, autologistic.grad, X = X, A = A, A1=A1, A2 = A2 ,Z = Z, detmat = detmat, detX = detX # gradient needed here, may be we cannot use it before we can find the gradient
		opt = try(optim(theta_previous, autooccu.ElogPLnotpr, NULL, X = X, A = A, A1=A1, A2 = A2 ,Z = Z, detmat = detmat, detX = detX, centered = centered,
              method = method, control = list(maxit = 300,trace=6,REPORT = 1), hessian = hessian),
              silent = F) # M step for EM
		if (class(opt) == "try-error")
		{
			coefficients = NULL
			fitted.values = NULL
			linear.predictors = NULL
			residuals = NULL
			convergence = NULL
			message = opt[1]
			value = NULL
			I.hat = NULL
			AIC = NULL
			err = T
			break
		}
		else
		{
			theta_new = opt$par
			relative_err = sum(abs(theta_new-theta_previous))/sum(abs(theta_previous))
			cond = relative_err>=toll # L1 norm convergence
			theta_previous = theta_new  # iterations here
		}
		print(paste("EM step No.",as.character(ww),"done with expected -log(L)=",as.character(opt$value)))
		#print("Expected -loglikelihood:")
		#print(opt$value)
	}
	
	if(!err)
	{
		convergence = opt$convergence & (!cond)
		if (convergence != 0)
			message = opt$message
		else
		message = NULL
		p = 2 * ncol(X) + ncol(detX[[1]]) + 3# first ncol(X) is environment beta, (ncol(X) + ncol(detX[[1]])) are detections, then 3 etas
		coefficients = opt$par
		names(coefficients)[p - 2] = "eta_spatial_spc1" # eta spatial
		names(coefficients)[p - 1] = "eta_spatial_spc2"
		names(coefficients)[p - 0] = "eta_interspecies"
		Xbeta = X %*% coefficients[1:ncol(X)]
		mu = exp(Xbeta)
		mu = mu / (1 + mu)
		autocovariate_1 = A %*% (Z - centered * mu)
		autocovariate_2 = A1 %*% (Z - centered * mu)
		covariate_interspecies = A2 %*% (Z - centered * mu)
		linear.predictors = Xbeta + coefficients[p - 2] * autocovariate_1 + coefficients[p - 1] * autocovariate_2 + coefficients[p - 0] * covariate_interspecies
		fitted.values = exp(linear.predictors)
		fitted.values = fitted.values / (1 + fitted.values)
		residuals = Z - fitted.values
		value = opt$value
		AIC = 2*value + 2 * length(coefficients)
	}
	if (hessian){
		I.hat = opt$hessian
	}
	
    object = list(coefficients = coefficients, fitted.values = fitted.values,
                  linear.predictors = linear.predictors, residuals = residuals, 
                  convergence = convergence, message = message, value = value, AIC=AIC)
    if (hessian && ! is.null(I.hat))
        object$I.hat = I.hat
    class(object) = "autooccu"
    object
}

## bootstrap to see the CI
autooccu.bootstrap.helper = function(dummy, X, A, A1, A2, theta, detmat,...)
{
    Z = rautoccu(X, A, A1, A2,theta) # sample a Z form the underlaying Ising model
	detmat = autooccu_sample.detection(theta, X, A, A1, A2, Z, detmat, detX) #sample detection history
    fit = autooccu.fit(X, A, A1, A2, detmat, detX, MCEMsample = 10000, centered = T ,hessian = F, method = 'BFGS')
    temp = rep(NA, length(theta))
    if (is.null(fit$convergence) || fit$convergence != 0)
        warning(fit$message)
    else
        temp = fit$coef
    temp
}

autologistic.bootstrap = function(X, A, A1, A2, theta, detmat, centered = T ,bootit)
{
    boot.sample = data.frame(matrix(NA, bootit, length(theta)))
    require(pbapply)
	  cat("\n")
	  flush.console()
	  gathered = pbapply::pblapply(1:bootit, autologistic.bootstrap.helper, X, A, A1, A2, theta, detmat, centered = T)
		
	for (j in 1:bootit){
		boot.sample[j, ] = gathered[[j]]
	}
    boot.sample	
 }
