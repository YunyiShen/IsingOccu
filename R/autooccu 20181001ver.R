
getGraph = function(distM,theta,int_range = "exp")
{
	p = length(theta)
	sites = nrow(distM)
	eta01 = theta[p-4]
	d01 = theta[p-3]
	eta02 = theta[p-2]
	d02 = theta[p-1]
	eta1 = theta[p]
	if(int_range=="arth"){
		D1 = eta01*as.matrix(1/((distM)^(2+d01))) # long range interactions dimension constant d01 and d02 (sigma in some physics papers regarding Ising model)
		D2 = eta02*as.matrix(1/((distM)^(2+d02)))
	}
	else{
	if(int_range=="exp"){
		D1 = eta01*as.matrix(exp(-exp(d01)*distM))
		D2 = eta02*as.matrix(exp(-exp(d02)*distM))
	}
	else{
	print("int_range must be exp or arth, will assume exp")
		D1 = eta01*as.matrix(exp(-exp(d01)*distM))
		D2 = eta02*as.matrix(exp(-exp(d02)*distM))
		}
	}
	I = matrix(0,nrow=sites,ncol=sites)
	diag(I) = eta1
	A = cbind(D1,I)
	B = cbind(I,D2)
	A = rbind(A,B)
	rm(B)
	diag(A)=0
	row.names(A)=colnames(A)
	return(A)
}

rautoccu = function(X, distM,theta,method = "MH",nIter=nIter,n=1,int_range = "exp") #distM is a distance matrix 
{	
	require(IsingSampler)
	p = length(theta)
	sites = nrow(distM)
	ncov = ncol(X)
	zeros = matrix(0,nrow=sites,ncol=ncov)
	beta1 = as.numeric( matrix(c(theta[1:(2*ncol(X))])))
	Xfull = cbind(rbind(X,zeros),rbind(zeros,X))
	thr = Xfull%*%beta1
	rm(Xfull)
	A = getGraph(distM,theta,int_range = int_range)
	as.matrix(IsingSampler(n=n,graph = A,thresholds=thr, responses = c(-1L, 1L),nIter=nIter,method=method)) # use IsingSampler
}

# sample detection given Z
autooccu_sample.detection = function(theta, X, distM, Z ,detmat, detX){
	p = length(theta)
	RN = matrix( runif(nrow(detmat)*ncol(detmat)),nrow = nrow(detmat),ncol = ncol(detmat) )
	nperiod = ncol(detmat) # detmat is the data of 0 and 1 for detections
	beta_det = theta[(2*ncol(X) + 1):(p-5)] # length(beta_det) = 2 * ncol(detX[[1]]) + 2 * ncol(X)  # beta for detections
	detDesign = lapply(detX,function(x,y){cbind(y,x)},y = X) # This is the full design matrix list of detection probability p at time
	npardet = ncol(detDesign[[1]])
	Xbeta_det1 = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = beta_det[1:npardet]) # Xbeta for detections
	Xbeta_det2 = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = beta_det[(npardet+1):(2*npardet)])
	P_det1 = lapply(Xbeta_det1,function(W){exp(W) / (1 + exp(W))}) # just a logistic regression for detection
	rm(Xbeta_det1)
	
	P_det2 = lapply(Xbeta_det2,function(W){exp(W) / (1 + exp(W))})
	rm(Xbeta_det2)
	
	P_det1 = (matrix(unlist(P_det1),nrow = nrow(X),ncol = nperiod)) # detection probability, row is site i col is period j
	P_det2 = (matrix(unlist(P_det2),nrow = nrow(X),ncol = nperiod))
	P_det = rbind(P_det1,P_det2)
	rm(P_det1)
	rm(P_det2)
	occustatus = matrix(rep((Z+1)/2,nperiod),nrow = length(Z),ncol=nperiod)
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

Pdet = function(envX, detmat, detX, theta)
{
	nperiod = ncol(detmat) # detmat is the data of 0 and 1 for detections
	beta_det = theta[(2*ncol(envX) + 1):(p-5)] # length(beta_det) = 2 * ncol(detX[[1]]) + 2 * ncol(X)  # beta for detections
	detDesign = lapply(detX,function(x,y){cbind(y,x)},y = envX) # This is the full design matrix list of detection probability p at time
	npardet = ncol(detDesign[[1]])
	Xbeta_det1 = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = beta_det[1:npardet]) # Xbeta for detections
	Xbeta_det2 = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = beta_det[(npardet+1):(2*npardet)])
	P_det1 = lapply(Xbeta_det1,function(W){exp(W) / (1 + exp(W))}) # just a logistic regression for detection
	rm(Xbeta_det1)
	P_det2 = lapply(Xbeta_det2,function(W){exp(W) / (1 + exp(W))})
	rm(Xbeta_det2)
	P_det1 = (matrix(unlist(P_det1),nrow = nrow(X),ncol = nperiod)) # detection probability, row is site i col is period j
	P_det2 = (matrix(unlist(P_det2),nrow = nrow(X),ncol = nperiod))
	P_det = rbind(P_det1,P_det2)
	return(Pdet)
}


# PSEUDO LIKELIHOOD of AUTOOCCU MODEL
# Pseudo Likelihood here, first version with non perfect detection is done Sept/22/2018 and checked NOT CHECKED YET
autooccu.logPL = function(theta, envX, distM, Z ,detmat, detX, int_range = "exp") # assume known Z, data should form as nrow(data)==nrow(X),ncol(data)=# of periods. detX is design matrix of detections, WITHOUT 1s, should be a list. length(detX) = ncol(data) = # of periods # A1 is dist matrix for species 2, all about first spc1 should be 0
{
	# deal with the underlaying Ising model
	# REMEMBER TO ADD 1s to envX, detX will automaticlly include 1s because envX in cbinded
	# nsite = nrow(X)/2	
	require(IsingSampler)
	p = length(theta)
	sites = nrow(distM)
	ncov = ncol(envX)
	zeros = matrix(0,nrow=sites,ncol=ncov)
	beta1 = as.numeric( matrix(c(theta[1:(2*ncol(envX))])))
	Xfull = cbind(rbind(envX,zeros),rbind(zeros,envX))
	thr = Xfull%*%beta1
	rm(Xfull)
	A = getGraph(distM,theta,int_range = int_range)
	logPL = ( IsingPL(x=Z, graph=A, thresholds=thr, beta=1, responses = c(-1L, 1L)))
	rm(A)
	
	# what need to do is just add detection history likelihoods here, while, Z should be EM.

	# now deal with non-perfect detection
	#nperiod = ncol(detmat) # detmat is the data of 0 and 1 for detections
	#beta_det = theta[(2*ncol(envX) + 1):(p-5)] # length(beta_det) = 2 * ncol(detX[[1]]) + 2 * ncol(X)  # beta for detections
	#detDesign = lapply(detX,function(x,y){cbind(y,x)},y = envX) # This is the full design matrix list of detection probability p at time
	#npardet = ncol(detDesign[[1]])
	#Xbeta_det1 = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = beta_det[1:npardet]) # Xbeta for detections
	#Xbeta_det2 = lapply(detDesign,function(w,beta1){w%*%beta1},beta1 = beta_det[(npardet+1):(2*npardet)])
	#P_det1 = lapply(Xbeta_det1,function(W){exp(W) / (1 + exp(W))}) # just a logistic regression for detection
	#rm(Xbeta_det1)
	#P_det2 = lapply(Xbeta_det2,function(W){exp(W) / (1 + exp(W))})
	#rm(Xbeta_det2)
	#P_det1 = (matrix(unlist(P_det1),nrow = nrow(X),ncol = nperiod)) # detection probability, row is site i col is period j
	#P_det2 = (matrix(unlist(P_det2),nrow = nrow(X),ncol = nperiod))
	#P_det = rbind(P_det1,P_det2)
	P_det = Pdet(envX, detmat, detX, theta)
	LP_Z1 = as.matrix(rowSums(detmat * log(P_det) + (1-detmat) * log(1-P_det)))
	LP_Z0 = as.matrix(log(1*(rowSums(detmat)==0) + 1e-13 * (1-(rowSums(detmat)==0)))) # I(data = 0), do not want err for those have detections
	logLdata = sum(as.numeric((Z+1)/2) * LP_Z1 + as.numeric(1-((Z+1)/2)) * LP_Z0) # likelihood of data 
	
	# total neg log likelihood
	-logPL-logLdata 
}

autooccu.logL.innorm = function(theta, envX, distM, Z ,detmat, detX, int_range = "exp"){ # the in-normalized log likelihood of Ising Model
    require(IsingSampler)
	p = length(theta)
	sites = nrow(distM)
	ncov = ncol(envX)
	zeros = matrix(0,nrow=sites,ncol=ncov)
	beta1 = as.numeric( matrix(c(theta[1:(2*ncol(envX))])))
	Xfull = cbind(rbind(envX,zeros),rbind(zeros,envX))
	thr = Xfull%*%beta1
	rm(Xfull)
	A = getGraph(distM,theta,int_range = int_range)
	negPot = thr%*%Z + t(Z)%*%A%*%Z
	
	P_det = Pdet(envX, detmat, detX, theta)
	LP_Z1 = as.matrix(rowSums(detmat * log(P_det) + (1-detmat) * log(1-P_det)))
	LP_Z0 = as.matrix(log(1*(rowSums(detmat)==0) + 1e-13 * (1-(rowSums(detmat)==0)))) # I(data = 0), do not want err for those have detections
	logLdata = sum(as.numeric((Z+1)/2) * LP_Z1 + as.numeric(1-((Z+1)/2)) * LP_Z0)

	return(negPot+logLdata)
}

# E-STEP of EM Algorithm
# E-STEP of EM Algorithm
autooccu.ElogPL = function(theta, X, distM, Z, detmat, detX, int_range = "exp") # E step for EM a par version using snow
{
  require('snow')
  cl = makeSOCKcluster(detectCores())
  n_sample = ncol(Z)
  clusterExport(cl,list("Z","theta","X","distM","detmat","detX","autooccu.logPL","int_range"),envir = environment())
  PLs = parCapply(cl=cl,x=Z,FUN='autooccu.logPL',theta = theta,
                            envX = X, distM=distM,  detmat = detmat, detX = detX,int_range=int_range)
  #w = as.array(Z)
	sumPL = sum(PLs)#calculate the E of likelihood
	EminuslogPL = sumPL/n_sample
	#print(EminuslogPL)
	stopCluster(cl)
	EminuslogPL
}

autooccu.ElogPLnotpr = function(theta, X, distM, Z, detmat, detX, int_range = "exp") # E step for EM
{
  n_sample = ncol(Z)
  PLs = apply(X=Z,MARGIN = 2,FUN='autooccu.logPL',theta = theta,
                  envX = X, distM=distM,detmat = detmat, detX = detX, int_range=int_range)
  sumPL = sum(PLs)#calculate the E of likelihood
  EminuslogPL = sumPL/n_sample
  EminuslogPL
}


# AUTOOCCU FIT HERE, MCEM is still under developing 2018/9/22
autooccu.fit = function(X, distM, detmat, detX, MCEMsample = 10000 ,hessian = T, method = 'BFGS',int_range = "exp")
{ 
  #require('snow')
  start = glm(as.numeric(rowSums(detmat)>0) ~ X - 1, family = binomial)$coef
	start_det = glm(detmat[,1] ~ cbind(X,detX[[1]]) - 1,family = binomial)$coef
	theta_previous = c(start,start_det,1,0,1,0,1)
	theta_previous = as.numeric(theta_previous)
	rownames(theta_previous) = NULL
	cond = T
	toll = 1e-5 # tolerance for EM algorithm
	Z_naive = as.matrix(rowSums(detmat)>0)
	LL0 = autooccu.logPL(theta_previous, X, distM, Z ,detmat, detX)
	print(paste("model started with -log(L) = ",as.character(LL0)))
	print("Initial MCEM")
	ww = 0
	while(cond){
	  ww = ww + 1
	  print("sampling Z...")
		Z = autooccuMCEM_sampleZ(theta_previous, X, distM, Z0 = Z_naive, detmat, detX, num_sample = MCEMsample) # USED SIMPLEST IMPLIMENT here, need to change to perfect sampler
		print("all clear, initial M step...")
		#opt = try(optim(theta_previous, autooccu.ElogPL, autologistic.grad, X = X, A = A, A1=A1, A2 = A2 ,Z = Z, detmat = detmat, detX = detX # gradient needed here, may be we cannot use it before we can find the gradient [it seems not so hard to get the gradient by hand, so try it]
		opt = try(optim(theta_previous, autooccu.ElogPLnotpr, NULL, X = X, distM = distM ,Z = Z, detmat = detmat, detX = detX, int_range = int_range,
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
			Graph = NULL
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
		Z = autooccuMCEM_sampleZ(theta_previous, X, distM, Z0 = Z, detmat, detX, num_sample = MCEMsample) # sample some Z to calculate the residue
		message = NULL
		p = 2 * ncol(X) + ncol(detX[[1]]) + 5# first ncol(X) is environment beta, (ncol(X) + ncol(detX[[1]])) are detections, then 5 etas
		coefficients = opt$par
		names(coefficients)[p - 4] = "sigma_spatial_spc1" # eta spatial
		names(coefficients)[p - 3] = "eta_spatial_spc1"
		names(coefficients)[p - 2] = "sigma_spatial_spc2" # eta spatial
		names(coefficients)[p - 1] = "eta_spatial_spc2"
		names(coefficients)[p - 0] = "eta_interspecies"
		Xbeta = X %*% coefficients[1:ncol(X)]
		mu = exp(Xbeta)
		mu = mu / (1/mu + mu)
		Graph = getGraph(theta_previous,distM)
		GraphZ = Graph %*% Z
		linear.predictors = Xbeta + apply(GraphZ,1,mean)
		fitted.values = exp(linear.predictors)
		fitted.values = fitted.values / (1/fitted.value + fitted.values)
		residuals =mean(apply(Z,2,function(Z,fitted.values){ Z - fitted.values},fitted.values = fitted.values))
		value = opt$value
		AIC = 2*value + 2 * length(coefficients)
		err = F
	}
	if (hessian){
		I.hat = opt$hessian
	}
	
    object = list(coefficients = coefficients, fitted.values = fitted.values,Graph = Graph,
                  linear.predictors = linear.predictors, residuals = residuals, 
                  convergence = convergence, message = message, value = value, AIC=AIC,err=err)
    if (hessian && ! is.null(I.hat))
        object$I.hat = I.hat
    class(object) = "autooccu"
    object
}

# THIS is the autooccu fitting function using Moller et al. 2006 sampler (if we can only use MCEM to do MPLE, then Baysian is much faster)
autooccu.fit.Moller = function(X,distM, detmat, detX, maxiter = 10000 , vars = rep(1,4*ncol(X)+2*ncol(detX[[1]])+9),int_range = "exp"){
	
	
	
	
	
	
}
## bootstrap to see the CI
autooccu.bootstrap.helper = function(dummy, X, A, A1, A2, theta, detmat,...)
{
    Z = rautoccu(X, A, A1, A2,theta) # sample a Z form the underlaying Ising model
	detmat = autooccu_sample.detection(theta, X, A, A1, A2, Z, detmat, detX) #sample detection history
    fit = autooccu.fit(X, A, A1, A2, detmat, detX, MCEMsample = 10000 ,hessian = F, method = 'BFGS')
    temp = rep(NA, length(theta))
    if (is.null(fit$convergence) || fit$convergence != 0)
        warning(fit$message)
    else
        temp = fit$coef
    temp
}

autologistic.bootstrap = function(X, A, A1, A2, theta, detmat, bootit)
{
    boot.sample = data.frame(matrix(NA, bootit, length(theta)))
    require(pbapply)
	  cat("\n")
	  flush.console()
	  gathered = pbapply::pblapply(1:bootit, autologistic.bootstrap.helper, X, A, A1, A2, theta, detmat)
		
	for (j in 1:bootit){
		boot.sample[j, ] = gathered[[j]]
	}
    boot.sample	
 }
