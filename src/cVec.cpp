#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::colvec cVec(const NumericMatrix & XX){       
  arma::mat X = Rcpp::as<arma::mat>(XX);
  int n = X.n_rows;
  int count =0;
  arma::colvec Y = arma::zeros(n*(n-1),1);
  for(int cc=0;cc<n;cc++){
    for(int rr=0;rr<n;rr++){
      if(cc != rr){
        Y(count) = X(rr,cc);
        count += 1;
      }
    }
  }
                       
return Y;
}
