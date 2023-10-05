#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat cBOmegaGibbs(const NumericMatrix & BBOm0,
                  const NumericMatrix & SS2,
                  const NumericMatrix & RR2){
  
  arma::mat BOm0 = Rcpp::as<arma::mat>(BBOm0);
  arma::mat s2 = Rcpp::as<arma::mat>(SS2);
  arma::mat r2 = Rcpp::as<arma::mat>(RR2);
                         
  int TT = s2.n_cols;
  int n = s2.n_rows;
  arma::mat BOm = BOm0;
                         
  for(int i=0;i<n;i++){
    for(int tt=0;tt<TT;tt++){
      BOm(0,0) += pow(s2(i,tt),2);
      BOm(1,1) += pow(r2(i,tt),2);
      BOm(0,1) += s2(i,tt)*r2(i,tt);
    }
  }
  BOm(1,0) = BOm(0,1);
                         
  return BOm;
}

