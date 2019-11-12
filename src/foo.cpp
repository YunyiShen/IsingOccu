// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <climits>
#include <complex>
#include <math.h>
using namespace Rcpp;



// [[Rcpp::export]]
double H(const arma::sp_mat& J, IntegerVector s, NumericVector h)
{
  double Res = 0;
  int N = J.n_rows;
  for (int i=0;i<N;i++)
  {
    Res -= h[i] * s[i];
    for (arma::sp_mat::const_col_iterator it = J.begin_col(i); it != J.end_col(i); ++it)
    {
      if (it.row()!=i) Res -= *it * s[i] * s[it.row()] * .5;
    }
  }
  return(Res);
}


// Likelihood without Z
// [[Rcpp::export]]
double f(IntegerMatrix Y, const arma::sp_mat& J, NumericVector h)
{
  double Res = 1;
  int Np = Y.nrow();
  int Ni = J.n_cols;
  IntegerVector s(Ni);
  for (int p=0;p<Np;p++)
  {
    for (int i=0;i<Ni;i++)
    {
      s[i] = Y(p,i);
    }
    Res *= exp(-1.0 * H(J, s, h));
  }
  return(Res);
}



// [[Rcpp::export]]
arma::mat AllStateCpp(
						 int N,
						 const IntegerVector &response){

	int n_possible = pow(2,N);
	int t=0;
	double Z=0;
	arma::mat res(n_possible,N);
	
	arma::mat temp(1,N);
	
	for(int i = 0; i<n_possible ; ++i){
		t = i;
		for(int j = 0; j<N ; ++j){
			temp[j] = response[t % 2];//use binary number coding
			t = t>>1;
		}
		res.row(i)=temp;
	}

	return(res);
	

}


// [[Rcpp::export]]
double PartitionCpp(
    arma::sp_mat graph,
    NumericVector thr,
    const IntegerVector &response){
  
  int N = graph.n_rows;
  int n_possible = pow(2,N);
  int t=0;
  double Z=0;
  arma::mat res(n_possible,N);
  
  IntegerVector temp(N);
  
  for(int i = 0; i<n_possible ; ++i){
    t = i;
    for(int j = 0; j<N ; ++j){
      temp[j] = response[t % 2];//use binary number coding
      t = t>>1;
    }
    Z+=exp(-H(graph,temp,thr));
  }
  
  return(Z);
  
  
}



