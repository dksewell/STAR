#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat cRoundMatrix(const NumericMatrix & XX){
  arma::mat X = Rcpp::as<arma::mat>(XX);
  double temp=0;
  int n= X.n_rows;
                              
  for(int i=0; i<n-1;i++){
    for(int j=i+1;j<n;j++){
      temp = 0.5*(X(i,j)+X(j,i));
      X(i,j)=temp;
      X(j,i)=temp;
    }
  }
  return X;
}