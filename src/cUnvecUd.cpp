#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat cUnvecUd(const NumericVector & YY){
  
  arma::colvec Y = Rcpp::as<arma::colvec>(YY);
  int n = (1+sqrt(1+8*Y.n_rows))/2;
  int count =0;
  arma::mat X = arma::zeros(n,n);
  
  for(int cc=1;cc<n;cc++){
    for(int rr=0;rr<cc;rr++){
      X(rr,cc) = Y(count);
      X(cc,rr) = X(rr,cc);
      count += 1;
    }
  }
                           
  return X;
}