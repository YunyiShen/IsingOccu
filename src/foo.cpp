// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <climits>
#include <string>
using namespace Rcpp;


arma::mat allNAs(const arma::mat & m1){
	arma::mat res;
	res = (m1!=m1);
	return(res);

}