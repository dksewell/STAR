#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat cBOmega(const NumericMatrix & BBOm0,
                  const NumericMatrix & SSigsr,
                  const NumericMatrix & mmusr){
  
  arma::mat BOm0 = Rcpp::as<arma::mat>(BBOm0);
  arma::mat Sigsr = Rcpp::as<arma::mat>(SSigsr);
  arma::mat musr = Rcpp::as<arma::mat>(mmusr);
                         
  int TT = musr.n_cols;
  int n = musr.n_rows/2;
  arma::mat BOm = BOm0;
                         
  for(int i=0;i<n;i++){
    BOm(0,0) += TT*Sigsr(i,i);
    BOm(1,1) += TT*Sigsr(n+i,n+i);
    BOm(0,1) += TT*Sigsr(i,n+i);
    for(int tt=0;tt<TT;tt++){
      BOm(0,0) += musr(i,tt)*musr(i,tt);
      BOm(1,1) += musr(n+i,tt)*musr(n+i,tt);
      BOm(0,1) += musr(i,tt)*musr(n+i,tt);
    }
  }
  BOm(1,0) = BOm(0,1);
                         
  return BOm;
}

// [[Rcpp::export]]
arma::mat cBOmega1(const NumericMatrix & BBOm0,
                  const NumericVector & SSigsr,
                  const NumericMatrix & mmusr){
  
  arma::mat BOm0 = Rcpp::as<arma::mat>(BBOm0);
  arma::mat musr = Rcpp::as<arma::mat>(mmusr);
                         
  int TT = musr.n_cols;
  int n = musr.n_rows/2;
  arma::mat BOm = BOm0;
  
  arma::cube Sigsr(SSigsr.begin(),2*n,2*n,TT);
  
                         
  for(int i=0;i<n;i++){
    for(int tt=0;tt<TT;tt++){
      BOm(0,0) += Sigsr.slice(tt)(i,i);
      BOm(1,1) += Sigsr.slice(tt)(n+i,n+i);
      BOm(0,1) += Sigsr.slice(tt)(i,n+i);
      BOm(0,0) += musr(i,tt)*musr(i,tt);
      BOm(1,1) += musr(n+i,tt)*musr(n+i,tt);
      BOm(0,1) += musr(i,tt)*musr(n+i,tt);
    }
  }
  BOm(1,0) = BOm(0,1);
                         
  return BOm;
}


// [[Rcpp::export]]
arma::mat cBOmegaSimple(const NumericMatrix & BBOm0,
                        const NumericVector & SSigsr,
                        const NumericMatrix & mmusr){
  
  arma::mat BOm0 = Rcpp::as<arma::mat>(BBOm0);
  arma::mat musr = Rcpp::as<arma::mat>(mmusr);
                         
  int TT = musr.n_cols;
  int n = musr.n_rows/2;
  arma::mat BOm = BOm0;
  
  arma::mat Sigsr = Rcpp::as<arma::mat>(SSigsr);
  
                         
  for(int i=0;i<n;i++){
    for(int tt=0;tt<TT;tt++){
      BOm(0,0) += Sigsr(i,i);
      BOm(1,1) += Sigsr(n+i,n+i);
      BOm(0,1) += Sigsr(i,n+i);
      BOm(0,0) += musr(i,tt)*musr(i,tt);
      BOm(1,1) += musr(n+i,tt)*musr(n+i,tt);
      BOm(0,1) += musr(i,tt)*musr(n+i,tt);
    }
  }
  BOm(1,0) = BOm(0,1);
                         
  return BOm;
}