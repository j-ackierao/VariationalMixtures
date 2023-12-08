#include <Rcpp.h>
using namespace Rcpp;

// Functions in Rcpp for variational mixture model with variable selection

// [[Rcpp::export]]
NumericMatrix rnkCalc(NumericMatrix rhonk, double N, double K) {
  NumericMatrix rnk(N, K );
  for (int n = 0; n < N; n++){
    double sum = 0; // Sum value
    for(int k = 0; k < K; k++){
      sum += rhonk(n, k);
    }
    for(int k = 0; k < K; k++){
      rnk(n, k) = rhonk(n, k) / sum;
    }
  }
  return rnk;
}

// [[Rcpp::export]]
NumericMatrix cmatrixCalc(NumericMatrix nullphi, NumericMatrix X, NumericVector c, double N, double D) {
  NumericMatrix cmatrix(N, D );
  for (int n = 0; n < N; n++){
    for(int d = 0; d < D; d++){
      cmatrix(n, d) = (1 - c(d)) * log(nullphi(d, X(n, d) - 1));
    }
  }
  return cmatrix;
}