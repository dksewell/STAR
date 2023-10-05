#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::colvec cSumXtXAndDraw(const NumericVector & XX,
                            const double & S2MRt){
  arma::colvec X = Rcpp::as<arma::colvec>(XX);
  int n = (1+sqrt(1+4*X.n_rows))/2;
  arma::colvec Y = arma::zeros(n*(n-1),1);
  double sqrtS2MRt = sqrt(S2MRt);

 for(int rr=0;rr<n-1;rr++){
   for(int cc=rr+1;cc<n;cc++){
      Y((n-1)*cc+rr) = X((n-1)*cc+rr) + X((n-1)*rr+cc-1);
      Y((n-1)*cc+rr) = S2MRt*Y((n-1)*cc+rr) + sqrtS2MRt*arma::as_scalar(arma::randn(1));
      Y((n-1)*rr+cc-1) = Y((n-1)*cc+rr);
    }
  }
  
  return Y;
}
