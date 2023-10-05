#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::colvec cVecUdSp(const arma::sp_mat & XX){
  
  int n = XX.n_rows;
  int count =0;
  arma::colvec Y = arma::zeros(n*(n-1)/2,1);
  
  for(int cc=1;cc<n;cc++){
    for(int rr=0;rr<cc;rr++){
      Y(count) = XX(rr,cc);
      count += 1;
    }
  }
  
  return Y;
}