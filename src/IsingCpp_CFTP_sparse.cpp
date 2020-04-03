// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <climits>
using namespace Rcpp;
using namespace arma;
// FUNCTIONS FOR EXACT SAMPLING //

// Inner function to resize list:
List resize( const List& x, int n ){
    int oldsize = x.size() ;
    List y(n) ;
    for( int i=0; i<oldsize; i++) y[i] = x[i] ;
    return y ;
}

// Inner function to simulate random uniforms in a matrix:
arma::mat RandMat(int nrow, int ncol)
 {
  arma::mat Res = arma::randu(nrow,ncol);
  return(Res);
 }

// Computes maximal and minimal probability of node flipping:
NumericVector PplusMinMax(int i, const arma::sp_mat& J, IntegerVector s, NumericVector h, double beta, IntegerVector responses)
{
  // The function computes the probability that node i is in Response 1 instead of 0, given all other nodes, which might be missing.
  // Output: minimal and maximal probablity
  
  NumericVector H0(2, h[i] * responses[0]); // relevant part of the Hamiltonian for state = 0
  NumericVector H1(2, h[i] * responses[1]); // relevant part of the Hamiltonian for state = 1
  
  NumericVector Res(2);
  
  //int N = J.n_rows;
  NumericVector TwoOpts(2);
  
  for (arma::sp_mat::const_col_iterator it = J.begin_col(i); it != J.end_col(i); ++it)
  {
    if (i != it.row())
    {
      if (s[it.row()] != INT_MIN)
      {
       H0[0] += *it * responses[0] * s[it.row()];
       H0[1] += *it * responses[0] * s[it.row()];
       H1[0] += *it * responses[1] * s[it.row()];
       H1[1] += *it * responses[1] * s[it.row()]; 
      } else 
      {
               
        TwoOpts[0] = *it * responses[1] * responses[0];
        TwoOpts[1] = *it * responses[1] * responses[1];

        if (TwoOpts[1] > TwoOpts[0])
        {
          H1[0] += TwoOpts[0];
          H1[1] += TwoOpts[1];
          
          H0[0] += *it * responses[0] * responses[0];
          H0[1] += *it * responses[0] * responses[1];
        } else 
        {
          H1[0] += TwoOpts[1];
          H1[1] += TwoOpts[0];          
          
          H0[0] += *it * responses[0] * responses[1];
          H0[1] += *it * responses[0] * responses[0];
        }
      }
    }
  }

  Res[0] = exp(beta * H1[0]) / ( exp(beta * H0[0]) + exp(beta * H1[0]) );
  Res[1] = exp(beta * H1[1]) / ( exp(beta * H0[1]) + exp(beta * H1[1]) );
  
  
  return(Res);
}
       
// Inner function:
IntegerVector IsingEx(const arma::sp_mat& graph, NumericVector thresholds, double beta, int nIter, IntegerVector responses, bool exact,
IntegerVector constrain)
{
  // Parameters and results vector:
  int N = graph.n_rows;
  IntegerVector state(N, INT_MIN);
  double u;
  NumericVector P(2);
  int maxChain = 100;
  List U(1);
  int minT = 0;
  bool anyNA = true;
    
  do
  { 
    // Resize U if needed:
    if (minT > 0)
    {
      U = resize(U, minT+1);
    }
    
    // Generate new random numbers:
    U[minT] = RandMat(nIter, N);
    
    // Initialize states:
    for (int i=0; i<N; i++)
    {
      if (exact)
      {
        state[i] = INT_MIN;
      } else 
      {
        state[i] = ifelse(runif(1) < 0.5, responses[1], responses[0])[0];
      }
    }    

    // START ALGORITHM
    for (int t=minT; t > -1;  t--)
    {
      for (int it=0;it<nIter;it++)
      {
        arma::mat Ucur = U[t];
        for (int node=0;node<N;node++)
        {
          u = Ucur(it, node);
          P = PplusMinMax(node, graph, state, thresholds, beta, responses);
          if (u < P[0])
          {
            state[node] = responses[1];
          } else if (u >= P[1])
          {
            state[node] = responses[0];
          } else 
          {
            state[node] = INT_MIN;
          }
        }
      }
    }
    
    anyNA = false;
    if (exact)
    {
      if (minT < maxChain)
      {
       for (int i=0; i<N; i++)
       {
        if (state[i] == INT_MIN)
        {
          anyNA = true;
        }
        } 
      } 
    }    
    minT++;
    
  } while (anyNA);

  // Rf_PrintValue(wrap(minT));
  return(state);
}



// FUNCTIONS FOR METROPOLIS SAMPLER //
double Pplus(int i, const arma::sp_mat& J, IntegerVector s, NumericVector h, double beta, IntegerVector responses)
{
  // The function computes the probability that node i is in Response 1 instead of 0, given all other nodes, which might be missing.
  // Output: minimal and maximal probablity
  
  double H0 = h[i] * responses[0]; // relevant part of the Hamiltonian for state = 0
  double H1 = h[i] * responses[1]; // relevant part of the Hamiltonian for state = 1

  //double Res;

  
  //int N = J.n_rows;
  
  
  for (arma::sp_mat::const_col_iterator it = J.begin_col(i); it != J.end_col(i); ++it)
  {
    if (i != it.row())
    {
       H0 += *it * responses[0] * s[it.row()];
       H1 += *it * responses[1] * s[it.row()];
    }
  }
  
  return(exp(beta * H1) / ( exp(beta * H0) + exp(beta * H1) ));// MH ratio here, need changing
}


IntegerVector IsingMet(const arma::sp_mat& graph, NumericVector thresholds, double beta, int nIter, IntegerVector responses,
IntegerVector constrain)
{
  // Parameters and results vector:
  int N = graph.n_rows;
  IntegerVector state =  ifelse(runif(N) < 0.5, responses[1], responses[0]);
  for (int i=0; i<N; i++)
  {
    if (constrain[i] != INT_MIN)
    {
      state[i] = constrain[i];
    }
  }
  double u;
  double P;
    
    // START ALGORITHM
    for (int it=0;it<nIter;it++)
    {
      for (int node=0;node<N;node++)
      {
        if (constrain[node] == INT_MIN)
        {
         u = runif(1)[0];
         P = Pplus(node, graph, state, thresholds, beta, responses);
          if (u < P)
         {
           state[node] = responses[1];
         } else 
         {
           state[node] = responses[0];
         } 
        }
      }
    }
   
  return(state);
}


///ISING PROCESS SAMPLER:
// [[Rcpp::export]]
IntegerMatrix IsingProcess(int nSample, const arma::sp_mat& graph, NumericVector thresholds, double beta, IntegerVector responses)
{
  // Parameters and results vector:
  int N = graph.n_rows;
  IntegerVector state =  ifelse(runif(N) < 0.5, responses[1], responses[0]);
  double u;
  double P;
  IntegerMatrix Res(nSample,N);
  int node;
    
    // START ALGORITHM
    for (int it=0;it<nSample;it++)
    {
      node = floor(R::runif(0,N));
        u = runif(1)[0];
        P = Pplus(node, graph, state, thresholds, beta, responses);
        if (u < P)
        {
          state[node] = responses[1];
        } else 
        {
          state[node] = responses[0];
        }
        for (int k=0; k<N; k++) Res(it,k) = state[k];
    }
   
  return(Res);
}





// OVERAL FUNCTION //
// [[Rcpp::export]]
IntegerMatrix IsingSamplerCpp(int n, const arma::sp_mat& graph, NumericVector thresholds, double beta, int nIter, IntegerVector responses, bool exact,
IntegerVector constrain)
{
  int Ni = graph.n_rows;
  IntegerMatrix Res(n,Ni);
  IntegerVector state(Ni);
  //IntegerVector constrainVec(Ni);
  if (exact)
  {
    for (int s=0;s<n;s++)
    {
      //for (int i=0;i<Ni;i++) constrainVec[i] = constrain(s,i);
      state = IsingEx(graph, thresholds, beta, nIter, responses, exact, constrain);
      for (int i=0;i<Ni;i++) Res(s,i) = state[i];
    }
  } else 
  {
    for (int s=0;s<n;s++)
    {
      //for (int i=0;i<Ni;i++) constrainVec[i] = constrain(s,i);
      state = IsingMet(graph, thresholds, beta, nIter, responses, constrain);
      for (int i=0;i<Ni;i++) Res(s,i) = state[i];
    }
  }
  
  return(Res);
}


// HELPER FUNCTIONS //
// Hamiltonian:
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


// [[Rcpp::export]]
double PartitionCpp(
    const arma::sp_mat & graph,
    const NumericVector & thr,
    const double & beta,
    const IntegerVector &responses){
  
  int N = graph.n_rows;
  int n_possible = pow(2,N);
  int t=0;
  double Z=0;
  
  IntegerVector temp(N);
  
  for(int i = 0; i<n_possible ; ++i){
    t = i;
    for(int j = 0; j<N ; ++j){
      temp[j] = responses[t % 2];//use binary number coding
      t = t>>1;
    }
    Z+=exp(-beta * H(graph,temp,thr));
  }
  
  return(Z);
  
  
}

// calculating state probability directly 
// [[Rcpp::export]]
double IsingStateProbCpp(const arma::vec & Z, 
                         const arma::sp_mat & graph,
                         const arma::vec & thr,
                         const IntegerVector &responses){
  double partion = PartitionCpp(graph,
                          Rcpp::NumericVector(thr.begin(),
                                              thr.end()),
                          1,responses);
  return(exp(- H(graph, IntegerVector(Z.begin(),
                                     Z.end()),                       
                       Rcpp::NumericVector(thr.begin(),
                                              thr.end())))/partion);
} // passed 2/4/2020



// [[Rcpp::export]]
double ColSumsarma(arma::mat A){
  return(sum( is_na(as<NumericMatrix>(wrap(A)))));
}


// [[Rcpp::export]]
double Pdet_Ising_single_siteCpp(
    const arma::mat & det_thr_temp, 
    const arma::vec & Z_temp, 
    const arma::mat & dethis, 
    const arma::mat & sppmat_det,
    const IntegerVector & responses){
  
  int nspp = sppmat_det.n_rows;
  if(sum(sum(Z_temp))==nspp*responses[0] || (sum( is_na(as<NumericMatrix>(wrap(dethis)))))==dethis.n_elem) return(0); // no species there, probability 1 to be all 0s.
  arma::uvec spp_exist = find(Z_temp == 1);
  arma::mat thr_exist = det_thr_temp.cols(spp_exist);
  arma::mat graph = sppmat_det.submat(spp_exist,spp_exist);
  arma::vec row_sums_dethis = sum(dethis.t()).t();
  
  LogicalVector has_obs_ornot = !is_na(NumericVector(row_sums_dethis.begin(),
                                                    row_sums_dethis.end())) ;
  arma::uvec has_obs = find(as<arma::vec>(has_obs_ornot));
  //arma::uvec has_obs = find(sum(dethis.t() != INT_MIN)==nspp);
  arma::mat dethis_temp = dethis.submat(has_obs,spp_exist); // delete those rows without observation and columns without species exist
  arma::mat thr_exist_temp = thr_exist.rows(has_obs);
  
  double res = 0;
  
  for(int i = 0; i < dethis_temp.n_rows ; ++i ){
    arma::vec dethis_this = dethis_temp.row(i).t();
    arma::vec thr_this = thr_exist_temp.row(i).t();
    res += log(IsingStateProbCpp( dethis_this,sp_mat(graph),
                              thr_this,
                              responses)+1e-15);
  
  
  }
  return(res);  
} // passed 20200403

// [[Rcpp::export]]
arma::mat extract_thrCpp(const int & which_site,// different from R version, start from 0
                         const List & det_thr,
                         const int & nspp,
                         const int & nperiod,
                         const int & nsite
                        ){
  arma::mat res(nperiod,nspp);
  
  for(int i = 0; i<nspp; ++i){
    arma::mat det_temp = det_thr[i];
    res.col(i) = det_temp.row(which_site).t();
  }
  return(res);
} // passed 20200403

// [[Rcpp::export]]
arma::mat Gibbs_Z_helperCpp(const arma::mat & Z_curr, // make sure Z_curr was column vector
                             const arma::vec & scans, // vector of non detection site/species numbers
                             const arma::mat & detmat,
                             const arma::sp_mat & A,
                             const arma::vec & thr,
                             const arma::mat & sppmat_det,
                             const List & det_thr,
                             const int & nsite,
                             const IntegerVector & responses
                            ) {
  int n_scans = scans.n_elem;
  int nspp = sppmat_det.n_rows;
  arma::mat Z_new = Z_curr;
  int nperiod = detmat.n_cols;
  
  for(int i = 0; i < n_scans ; ++i){
    int scan = scans(i); // the one working on this scan, should be one that was -1 in absolute
    double Ham_plus1 =  (A.row(scan-1) * Z_curr)(0,0) + thr(scan-1); // prior part, Ising model
    
    int which_site = scan % nsite;
    
    if(which_site==0) which_site = nsite;
    int which_spp = (scan-which_site)/nsite;
    uvec adder = regspace<uvec>(0, nspp-1);
    uvec which_row = adder*nsite+which_site-1;// start from 0
    //Rcout << which_row<<std::endl;
    
    arma::vec Z_temp = Z_curr.rows(which_row);//get relavent Z
    
    arma::mat dethist = detmat.rows(which_row);
    dethist = dethist.t();// row as period, columns as species, note that all species were included here.
    arma::mat det_thr_temp = extract_thrCpp(which_site-1, //extract_thr ask index to start with 0,
                                det_thr,nspp,nperiod,nsite);
    
    Z_temp(which_spp) = responses[1];
    double Pplus1 = Pdet_Ising_single_siteCpp(det_thr_temp, 
                                Z_temp, dethist, sppmat_det,responses);
    Z_temp(which_spp) = responses[0];
    double Pminus1 = Pdet_Ising_single_siteCpp(det_thr_temp, 
                                Z_temp, dethist, sppmat_det,responses);

    
    double PZiplus1=(exp(Ham_plus1+Pplus1))/
    (exp(-Ham_plus1+Pminus1)+exp(Ham_plus1+Pplus1));
    
    Z_new.row(scan-1) = ifelse(runif(1) < PZiplus1, 
                         responses[1], responses[0])[0];
  }
  return(Z_new);
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


// VECTOR VERSIONS //


// Hamiltonian:
// [[Rcpp::export]]
double Hvec(IntegerVector s, NumericVector Theta, int N)
{
  double Res = 0;
  int c=0;
  for (int i=0;i<N;i++)
  {
    Res -= Theta[c] * s[i];
    c++;
  }
  for (int i=0;i<N;i++)
  {
    for (int j=i; j<N;j++)
    {
     if (j!=i) 
     {
       Res -= Theta[c] * s[i] * s[j]; 
       c++;
     }
    }
  }
  return(Res);
}

// Likelihood without Z
double fvec(IntegerMatrix Y, NumericVector Theta)
{
  double Res = 1;
  int Np = Y.nrow();
  int Ni = Y.ncol();
  IntegerVector s(Ni);
  for (int p=0;p<Np;p++)
  {
    for (int i=0;i<Ni;i++)
    {
      s[i] = Y(p,i);
    }
    Res *= exp(-1.0 * Hvec(s, Theta, Ni));
  }
  return(Res);
}

// Log-Likelihood without Z
double fveclog(IntegerMatrix Y, NumericVector Theta)
{
  double Res = 1;
  int Np = Y.nrow();
  int Ni = Y.ncol();
  IntegerVector s(Ni);
  for (int p=0;p<Np;p++)
  {
    for (int i=0;i<Ni;i++)
    {
      s[i] = Y(p,i);
    }
    Res -= Hvec(s, Theta, Ni);
  }
  return(Res);
}
//

// Function to compute expected values:
// [[Rcpp::export]]
NumericVector expvalues(IntegerMatrix x){
  // Sample size:
  int N = x.nrow();
  // Number of nodes:
  int P = x.ncol();
  int nPar = P + P*(P-1)/2;
  // Results vector:
  NumericVector Res(nPar, 0.0);
  
  int par = 0;
  
  // Fill
  while (par < nPar){
    
    // Means:
    for (int j=0; j<P;j++){
      double mean = 0;
      for (int k=0; k<N; k++){
        mean += x(k,j) ;
      }
      Res[par] = mean / N ;
      par++;
    }
    
    // Covariances:
    for (int i=0; i<P; i++){
      for (int j=i; j<P;j++){
        if (i != j){
         double squared = 0;
          for (int k=0; k<N; k++){
            squared += x(k,i) * x(k,j) ;
          }
          Res[par] = squared / N;
          par++; 
        }
      }
    }
  }
  
  return(Res);
}

// Function to obtain thresholds from vector:
// [[Rcpp::export]]
NumericVector vec2Thresh(NumericVector vec, int P){
  NumericVector Res(P);
  
  for (int i=0; i<P; i++){
    Res[i] = vec[i];
  }
  
  return(Res);
}

// [[Rcpp::export]]
arma::sp_mat vec2Graph(NumericVector vec, int P){
  arma::sp_mat Res(P, P);
  
  int par = P;
  
  for (int i=0; i<P; i++){
      for (int j=i; j<P; j++){
        if (i != j){
           Res(i,j) = Res(j,i) = vec[par];   
           par++;
        }
      }
  }
  
  return(Res);
}

// Main optimisation function:
// [[Rcpp::export]]
NumericVector Broderick2013(
  IntegerMatrix x, // Data matrix
  int M, // Number of samples  to draw
  int T, // Number of iterations
  int nIter, // Temporary: number of sequences, replace with convergence test
  IntegerVector responses
  )
{
  // Sample size:
  int N = x.nrow();
  // Number of nodes:
  int P = x.ncol();
  // Number of parameters:
  int nPar = P + P*(P-1)/2;
  // Current estimtes and new estimates:
  NumericVector curEsts(nPar, 0.0);
  NumericVector newEsts(nPar, 0.0);
  
    // Dummy constraints mat (ugly, should be removed):
  IntegerMatrix cons(M, P);
  std::fill(cons.begin(), cons.end(), INT_MIN);

  // Observed statistics:
  NumericVector obsStats = expvalues(x);
  
  // Thresholds to mimic margins:
  for (int i=0; i<P; i++){
    curEsts[i] = (responses[1] - responses[0]) *log(obsStats[i]);
  }

  double step = 1;

  // Start iterating:
  for (int s=0; s<nIter; s++){
    // Set new ests:
    for (int i=0;i<nPar;i++){
      curEsts[i] = newEsts[i];
    }
    
    // Generate monte carlo samples:
    IntegerMatrix Samples =  IsingSamplerCpp(M, vec2Graph(curEsts, P), vec2Thresh(curEsts, P), 1.0, 1000, responses, false,cons);
    
    // Statistics:
    NumericVector sampStats = expvalues(Samples);
    
    for (int t=0; t<T; t++){
      // For each statistic, move up if expected value too low, down if expected value too high:
      for (int par=0;par<nPar;par++){
        // Estimate sampStat:
        double stat = (sampStats[par] * exp(-(newEsts[par] - curEsts[par]) * sampStats[par]) ) / (exp(-(newEsts[par] - curEsts[par]) * sampStats[par]) );
        
        
        // too high:
        if (stat > obsStats[par]){
          newEsts[par] -= step;
        } else {
          // Too low:
          newEsts[par] += step;
        }
      }
    }

    step *= 0.5;
  }
  return(newEsts);
}





