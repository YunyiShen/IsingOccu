// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <climits>
#include <string>
#include <math.h> // to use pow
using namespace Rcpp;

double Murray_ratio_Ising_occuDetCpp(const List & MRF_curr,
							   const List & MRF_prop,
							   const List & log_pi, // this is the prior density
							   const arma::mat & Z,
							   const arma::mat & Z_temp,
							   const arma::mat & Z_prop,
							   const arma::mat & beta_det_curr,
							   const arma::mat & beta_det_curr,
							   const List & detX,
							   const arma::sp_mat sppmat_det_curr,
							   const arma::sp_mat sppmat_det_prop
								){







}





