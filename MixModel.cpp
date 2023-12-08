#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>
using namespace Rcpp;

// Functions in Rcpp for original variational mixture model

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


