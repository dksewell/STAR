#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double cUpperMRt(const NumericVector & XX){       
  arma::colvec X = Rcpp::as<arma::colvec>(XX);
  int n = (1+sqrt(1+4*X.n_rows))/2;
  double sum = 0;
  for(int rr=0;rr<n-1;rr++){
    for(int cc=rr+1;cc<n;cc++){
        sum = sum + pow(X((n-1)*cc+rr),2);
    }
  }
  return sum;
}