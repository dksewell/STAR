#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
  
// [[Rcpp::export]]
arma::sp_mat cUnvecSp(const NumericVector & YY){
  
  arma::colvec Y = Rcpp::as<arma::colvec>(YY);
  int n = (1+sqrt(1+4*Y.n_rows))/2;
  int count =0;
  arma::sp_mat X = arma::sp_mat(n,n);
  for(int cc=0;cc<n;cc++){
    for(int rr=0;rr<n;rr++){
      if(cc != rr){
        X(rr,cc)=Y(count);
        count += 1;
      }
    }
  }
  return X;
}