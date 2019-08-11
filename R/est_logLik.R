est_logLik_occu = function(theta_0, Y, theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra="nn",int_range_inter="exp",Z_vec){
	# Descombes 1999 paper, calculate Z_theta/Z_theta0, using samples Y from theta_0, I will use posterior mean.
	n_theta = length(theta_0)
	theta_diff = lapply(1:n_theta,function(i,theta_0,theta){theta_0[[i]]-theta[[i]]},theta_0,theta)
	partition_ratio = mean( 
		exp(-Hamiltonian(theta_diff,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Y))
	)# ratio of partitioning function: Z_theta/Z_theta0
	
	# this is the estimation of a single theta, we need to loop through every sample
	-Hamiltonian(theta,envX,distM,link_map,dist_mainland,link_mainland,int_range_intra,int_range_inter,Z_vec)-log(partition_ratio)
}



