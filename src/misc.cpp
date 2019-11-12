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


arma::mat HamiltonianCpp(const List & MRF , const arma::mat & Z_vec){
	nrep = Z_vec.n_cols;
	arma::mat Ham(nrep,1)
	arma::sp_mat A = MRF["A"];
	arma::mat thr = MRF["thr"];
	for(int i=0 ; i<nrep ; ++i){
		Ham.row(i).col(0) = H(A,Z.col(i),thr);
	}
	return(Ham);
}


double IsingStateProbCpp(const arma::mat &s,
						 const arma::mat &graph,
						 const arma::mat &thr,
						 const arma::mat &response){

	int N = graph.n_rows;
	int n_possible = exp(log(2)*N);
	int t=0;
	double Z=0;
	
	arma::mat temp(N,1)
	
	for(int i = 0; i<n_possible ; ++i){
		t = i;
		for(j = 0; j<N ; ++j){
			temp.row(j) = response.row(t mod 2);//use binary number coding
			t = t>>1;
		}
		Z += exp(-H(graph,temp,thr));
	}

	return(exp(-H(graph,s,thr))/Z);
	

}



double Isingdet_single_siteCpp(const arma::mat & thr,
							   const arma::mat & Z,
							   const arma::mat & det,
							   const arma::mat & spp_mat_det,
							   const arma::mat & response){
	uvec ext = find(Z-min(Z));//find which species exist here;
	arma::mat graph = spp_mat_det.submat(ext,ext);//the graph of coexisted species
	arma::mat thr_ext = thr.cols(ext).t();//col matrix of thr
	
	double res;
	res = IsingStateProbCpp(det,graph,thr,response);
	return(res);

}


double Pdet_Ising_repCpp(const int & nrep, // number of repeats
					  const int & nperiod, // number of detection periods per rep
					  const arma::mat & envX,
					  const List & detX, // make sure this design matrix for det is well prepared in R
					  const arma::mat & beta_det,
					  const arma::mat & sppmat_det,
					  const arma::mat & Z, // columes as repeats
					  const List & detmat, // elements are detections at certain rep
					  const arma::mat & response
					 ){
	double Pdets=0;
	double P_rep = 0;
	double P_site = 0;
	int n_site = Z.n_rows; // number of sites
	int n_spp = sppmat_det.n_rows;
	int n_par = beta_det.n_rows/n_spp;
	
	arma::mat beta_det_m = beta_det.reshape(n_par,n_spp)
	
	for(int i = 0; i<nrep ; ++i){
		List detX_i = detX[i]; // detection design matrix at repeat i;
		arma::mat detmat_i = detmat[i];// detection history matrix
		arma::mat Z_i = Z.col(i);
		Z_i = Z_i.reshape(n_site,n_spp);
		for(int j = 0 ; j<n_period ; ++j){
			
			arma::mat detX_j = detX_i[j]; // detection design matrix at repeat i and period j 			
			arma::mat thr_j = detX_j * beta_det_m; // this will give the external field in detection at each site, for each species			
			for(int w = 0 ; w<n_site ; ++w){
				if( sum(detmat_i.col(j))!=INT_MIN ){
				P_site += Isingdet_single_siteCpp(thr_j.row(w), Z_i.row(w)
											   , detmat_i.col(j), sppmat_det,response);
				}// at least have observations
							
			}
			P_rep += P_site;
			P_site = 0;
		}
		Pdets += P_rep;	
		P_rep = 0;
	}
	return(Pdets);
}















