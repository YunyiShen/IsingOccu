// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;

#include <RcppArmadillo.h>
#include <iostream>

using namespace std;
using namespace arma;

// [[Rcpp::export]]
double bmse(const vec& vals)
{
  int N = vals.size();
  double result;
  if (N < 10)
    result = -1;
  else
  {
    int b = floor(sqrt(N));
    int a = floor(N / b);
    colvec Ys(a);
    colvec temp(b);
    for (int i = 0; i < a; i++)
    {
      temp = vals.subvec(i * b, (i + 1) * b - 1);
      Ys[i] = mean(temp);
    }
    double muhat = mean(Ys);
    double sigmahatsq = b * as_scalar(sum(pow((Ys - muhat), 2))) / (a - 1);
    result = sqrt(sigmahatsq / N);
  }
  return result;
}// this comes from nagspatial package, for bootstrap

double autologis_logL(const mat& X, const mat& A, const mat& A1, const mat& A2, const colvec& Z ,const colvec& theta)
{
	colvec beta;
	colvec Xbeta;
	colvec mu;
	int n_env;
	n_env = X.n_cols;
	beta = theta.subvec(0,n_env - 1);
	Xbeta = X * beta;
	mu = exp(Xbeta)
	


}// this is simple -log(L) for the underlaying Ising model

double autodet_logL(const colvec& Z ,const colvec& theta, List detX, const mat& det_history)
{


}// this is for the detection -log(L) given Z



RCPP_MODULE(Pseudolikelihood)
{
  Rcpp::function("bmse", &bmse);
}

