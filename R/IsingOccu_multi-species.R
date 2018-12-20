getintralayerGraph = function(distM,eta,d,int_range = "exp",nspp)
{
	#eta = eta[1:nspp]
	A = list()
	if(int_range=="arth"){
		for(i in 1:nspp){
			A[i] = eta[i]*as.matrix(1/((distM)^(2+d[i])))
		}
	}
	else{
		if(int_range=="exp"){
			for(i in 1:nspp){
				A[i] = eta[i]*as.matrix(exp(-abs(d[i])*distM))
			}
		}
		else{
			print("int_range must be exp or arth, will assume exp")
			for(i in 1:nspp){
				A[i] = eta[i]*as.matrix(exp(-abs(d[i])*distM))
			}
		}
	}
	return(A)
}

IsingOccu_multispp.logL.innorm = function(beta,eta_intra,d,eta_inter, envX, distM, Z ,detmat, detX, detbeta,int_range = "exp",nspp){ # the in-normalized log likelihood of IsingOccu Model beta is matrix here
# eta_inter should be a symmetric matrix, eta_intra should be vector
# detmat here should be a matrix, with 1:nspp row for first spp
    #require(IsingSampler)
	#nspp = ncol(Z)
    Z_vec=matrix(Z,nrow=length(Z),ncol=1) # Z can be a vector with 1:nsite be first species.
	
	sites = nrow(distM)
	ncov = ncol(envX)
	zeros = matrix(0,nrow=sites,ncol=ncov)
	#beta1 = as.numeric( matrix(c(theta[1:(2*ncol(envX))])))
	#Xfull = cbind(rbind(envX,zeros),rbind(zeros,envX))
	thr = X%*%beta # a matrix
	#rm(Xfull)
	A = getGraph(distM,eta_intra,d,int_range = int_range)
	negPot = 0
	for(i in 1:nspp){ # intralayer terms:
		negPot = negPot + t(thr[,i])%*%Z_vec[1:nsite + (i-1) * nsite] + 0.5*t(Z_vec[1:nsite + (i-1) * nsite])%*%A[[i]]%*%Z_vec[1:nsite + (i-1) * nsite]
	}
	for(i in 2:nspp-1){
		for (j in (i+1):nspp){
			negPot = negPot + eta_inter[i,j] * t(Z_vec[1:nsite + (i-1) * nsite])%*%(Z_vec[1:nsite + (j-1) * nsite])
		}
	}
	
	beta_det = matrix(detbeta,nrow = length(detbeta),ncol = 1)
	P_det = Pdet_multi(envX, detmat, detX, beta_det)
	LP_Z1 = as.matrix(rowSums(detmat * log(P_det) + (1-detmat) * log(1-P_det)))
	LP_Z0 = as.matrix(log(1*(rowSums(detmat)==0) + 1e-13 * (1-(rowSums(detmat)==0)))) # I(data = 0), do not want err for those have detections
	logLdata = sum(as.numeric((Z_vec+1)/2) * LP_Z1 + as.numeric(1-((Z_vec+1)/2)) * LP_Z0)

	return(negPot+logLdata)
}

Pdet_multi = function(envX, detmat, detX, beta_det) # likelihood given Z and detections
{
	# this is still 2 spp case, need to change to multi case
	nperiod = ncol(detmat) # detmat is the data of 0 and 1 for detections
	 # length(beta_det) = 2 * ncol(detX[[1]]) + 2 * ncol(X)  # beta for detections
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
	P_det = rbind(P_det1,P_det2) # col is period, row from 1:nsite for spp1 and nsite+1:nsite for spp2
	return(P_det)
}


