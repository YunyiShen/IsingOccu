negHamiltonian_posterior = function(theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp",Z_vec){
	H = matrix(0,1,2)
	namestr = H
	beta_occu = theta$beta_occu # this will be a matrix for cols are species
	beta_det = theta$beta_det
	eta_intra = theta$eta_intra # intra spp, intra island if apply
	d_intra = theta$d_intra
	eta_inter = theta$eta_inter # assume there is one 
	d_inter = theta$d_inter
	spp_mat = theta$spp_mat
	eta_inter = theta$eta_inter
	nspp = sqrt(length(spp_mat))
	spp_mat = matrix(spp_mat,nspp,nspp)
	ncov = ncol(envX)
	nsites = nrow(envX)
	
	A = getintralayerGraph(distM,link_map$intra,eta_intra,d_intra,int_range = int_range_intra,spp_mat)
	B = getintralayerGraph(distM,link_map$inter,eta_inter,d_inter,int_range = int_range_inter,spp_mat)
	k = 1
	for(i in 1:nspp){ # intralayer terms:
		#H[[i]]=list()
		for(j in 1:ncov){ # species i, envj
			H[k] = sum(beta_occu[j+(i-1)*ncov]*envX[,j]*Z_vec[1:nsites + (i-1) * nsites,])
			
			namestr[k]=paste0("spp_",i,"_beta_",j)
			k = k + 1
		}
		H[k] = .5 * t(as.matrix( Z_vec[1:nsites + (i-1) * nsites,]))%*%A[[i]]%*%as.matrix( Z_vec[1:nsites + (i-1) * nsites,]) # intra-island 
		
		namestr[k]=paste0("spp_",i,"_Intra_Island")
		k = k + 1
		H[k] = .5 * t(as.matrix( Z_vec[1:nsites + (i-1) * nsites,]))%*%B[[i]]%*%as.matrix( Z_vec[1:nsites + (i-1) * nsites,]) # inter-island
		
		namestr[k]=paste0("spp_",i,"_Inter_Island")
		k = k + 1
		thr_mainland = mainland_thr(dist_mainland,link_mainland,eta_inter[i],d_inter[i],int_range_inter)
		H[k] = sum(thr_mainland*Z_vec[1:nsites + (i-1) * nsites,])

		namestr[k]=paste0("spp_",i,"_Mainland")
		k = k + 1
	}
	for(i in 2:nspp-1){
		for(j in (i+1):nspp){
			H[k] = spp_mat[i,j] * (t(Z_vec[1:nsites + (i-1) * nsites,])%*%(Z_vec[1:nsites + (j-1) * nsites,]))
			namestr[k]=paste0("Cor_spp_",i,"_spp_",j)
			k = k + 1
			
		}
	}
	names(H)=namestr
	return(H)
}
  # passed 2019/3/19
