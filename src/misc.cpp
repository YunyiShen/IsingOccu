// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <climits>
#include <string>
using namespace Rcpp;



// Three kind of graph
List getintralayerGraphExpCpp(const arma::sp_mat & distM,
						   const arma::sp_mat & link_map,
						   const arma::mat & eta,
						   const arma::mat & d,
						   const int & nspp){
	List A(nspp);
	arma::sp_mat At;
	
	for(int i=0;i<nspp;++i){ // main loop for intra layer graphs
		At = eta.row(i).col(0) * exp(-exp(d[i]*distM)) % link_map;
		At.diag() = 0;
		A[i] = At;
	}
	return(A);
}

List getintralayerGraphNNCpp(
						   const arma::sp_mat & link_map,
						   const arma::mat & eta,
						   const int & nspp){
	List A(nspp);
	arma::sp_mat At;
	
	for(int i=0;i<nspp;++i){ // main loop for intra layer graphs
		At = eta.row(i).col(0) * link_map;
		At.diag() = 0;
		A[i] = At;
	}
	return(A);
}


List getintralayerGraphArthCpp(const arma::sp_mat & distM,
						   const arma::sp_mat & link_map,
						   const arma::mat & eta,
						   const arma::mat & d,
						   const int & nspp){
	List A(nspp);
	arma::sp_mat At;
	
	for(int i=0;i<nspp;++i){ // main loop for intra layer graphs
		At = eta.row(i).col(0) * (1/distM^(2+d.row(i).col(0))) % link_map;
		At.diag() = 0;
		A[i] = At;
	}
	return(A);
}

// Mainland island process
arma::mat mainland_thrNNCpp(
						   const arma::mat & link_mainland,
						   const double & eta){
	return(eta * link_mainland);

}

arma::mat mainland_thrExpCpp(const arma::mat & dist_mainland,
						   const arma::mat & link_mainland,
						   const double & eta,
						   const double & d){
	return(eta * exp(-exp(d*dist_mainland)) % link_mainland);

}

arma::mat mainland_thrArthCpp(const arma::mat & dist_mainland,
						   const arma::mat & link_mainland,
						   const double & eta,
						   const double & d){
	return(At = eta * (1/dist_mainland^(2+d)) % link_map);

}



// Full graph
arma::sp_mat getfullGraphCpp(const List & A_ex,
							const List & A_in,
							const arma::mat & spp_mat){
	int nspp = spp_mat;
	arma::sp_mat At = A_ex[0];
	int nsite = At.n_rows;
	arma::sp_mat A(nspp*nsite,nspp*nsite);
	
	for(int i = 0; i<(nspp-1) ; ++i){//main loop for full graph
		arma::sp_mat A_ex_temp = A_ex[i]; // inter island
		arma::sp_mat A_in_temp = A_in[i]; // intra island
		A.submat(i*nsite,i*nsite,(i+1)*nsite-1,(i+1)*nsite-1) = A_ex_temp + A_in_temp;
		for(int j = i+1 ; j<nspp ; ++j){
			A.submat((i-1)*nsite,(j-1)*nsite,i*nsite-1,j*nsite-1).diag()=spp_mat(i,j);// inter-specific interactions
			A.submat((j-1)*nsite,(i-1)*nsite,j*nsite-1,i*nsite-1).diag()=spp_mat(j,i);
		}
		//last species
		arma::sp_mat A_ex_temp = A_ex[nspp - 1]; 
		arma::sp_mat A_in_temp = A_in[nspp - 1]; 
		A.submat((nspp - 1)*nsite,(nspp - 1)*nsite,(nspp)*nsite-1,(nspp)*nsite-1) = A_ex_temp + A_in_temp;
		
	}

	return(A);
}


List getMRF(const List & theta, 
			const arma::mat & envX,
			const arma::mat & distM,
		    const List & link_map,
		    const arma::mat & dist_mainland,
		    const arma::mat & link_mainland,
		    const string int_range_intra,
		    const string int_range_inter){
	int nsite = envX.n_rows;
	arma::mat beta_occu = theta["beta_occu"];
	arma::mat eta_intra = theta["eta_intra"]; // intra spp, intra island if apply
	arma::mat d_intra = theta["d_intra"];
	arma::sp_mat link_intra = link_map["intra"];
	arma::mat eta_inter = theta["eta_inter"]; // assume there is a 
	arma::mat d_inter = theta["d_inter"];
	arma::sp_mat link_inter = link_map["inter"];
	arma::mat spp_mat = theta["spp_mat"];
	int nspp = spp_mat.n_rows;
	int npar = envX.n_cols;
	if(int_range_intra.compare("exp")==0){
		List A_in = getintralayerGraphExpCpp(distM,
						    link_intra,
						    eta_intra,
						    d_intra,
						    nspp);
	}
	if(int_range_intra.compare("nn")==0){
		List A_in = getintralayerGraphNNCpp(
						    link_intra,
						    eta_intra,
						    
						    nspp);
	}
	if(int_range_intra.compare("arth")==0){
		List A_in = getintralayerGraphArthCpp(distM,
						    link_intra,
						    eta_intra,
						    d_intra,
						    nspp);
	}
	
	if(int_range_inter.compare("exp")==0){
		List A_ex = getintralayerGraphExpCpp(distM,
						    link_inter,
						    eta_inter,
						    d_inter,
						    nspp);
	}
	if(int_range_inter.compare("nn")==0){
		List A_ex = getintralayerGraphNNCpp(
						    link_inter,
						    eta_inter,
						    
						    nspp);
	}
	if(int_range_inter.compare("arth")==0){
		List A_ex = getintralayerGraphArthCpp(distM,
						    link_inter,
						    eta_inter,
						    d_inter,
						    nspp);
	}

	arma::sp_mat A = getfullGraphCpp(A_ex,A_in,spp_mat);// get full graph
	
	delete A_ex;
	delete A_in;
	
	arma::mat thr(nspp*nsite,1);
	arma::mat mainland_thr;
	
	for(int i=0;i<nspp;++i){
		if(int_range_inter.compare("exp")==0){
			mainland_thr = mainland_thrExpCpp(dist_mainland,
						    link_mainland,
						    eta_inter,
						    d_inter,
						    nspp);
		}
		if(int_range_inter.compare("nn")==0){
			mainland_thr = mainland_thrNNCpp(
						    link_mainland,
						    eta_inter,
						    
						    nspp);
		}
		if(int_range_inter.compare("arth")==0){
			mainland_thr = mainland_thrArthCpp(dist_mainland,
						    link_mainland,
						    eta_inter,
						    d_inter,
						    nspp);
		}
		thr.rows(i*nsite,(i+1)*nsite-1) = envX * beta_occu.rows(i*npar,(i+1)*npar-1) + 
			mainland_thr;
	}
	
	return(Rcpp::List::create(
	        Rcpp::Named("A") = A ,//Graph
	        Rcpp::Named("thr") = thr,// external field 
	));		
}





