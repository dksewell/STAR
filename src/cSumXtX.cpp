#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::colvec cSumXtX(const NumericVector & XX){
  arma::colvec X = Rcpp::as<arma::colvec>(XX);
  int n = (1+sqrt(1+4*X.n_rows))/2;
  arma::colvec Y = arma::zeros(n*(n-1),1);

 for(int cc=0;cc<n;cc++){
   for(int rr=0;rr<n;rr++){
      if(cc != rr){
        int ind = (rr>cc);
        int ind2 = (cc>rr);
        Y((n-1)*cc+rr-ind) = X((n-1)*cc+rr-ind) + X((n-1)*rr+cc-ind2);
      }
    }
  }
  
  return Y;
}
