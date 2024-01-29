#include <RcppArmadillo.h> // new 'lighter' header
#include <boost/math/special_functions/digamma.hpp>

// Rcpp Armadillo functions for profile regression with categorical response

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat CpriorepsCalc(arma::cube prioreps, double K, double D, arma::vec nCat){
  arma::mat v(K, D);
  for (int d = 0; d < D; d++){
    double varCat = nCat(d);
    for (int k = 0; k < K; k++){
      double sum1 = 0; // Sum value
      for(int j = 0; j < varCat; j++){
        sum1 += prioreps(k, j, d);
      }
      double sum2 = 0;
      for (int j = 0; j < varCat; j++){
        sum2 += lgamma(prioreps(k, j, d));
      }
      v(k, d) = lgamma(sum1) - sum2;
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat CpostepsCalc(arma::cube eps, double K, double D, arma::vec nCat){
  arma::mat v(K, D);
  for (int d = 0; d < D; d++){
    double varCat = nCat(d);
    for (int k = 0; k < K; k++){
      double sum1 = 0; // Sum value
      for(int j = 0; j < varCat; j++){
        sum1 += eps(k, j, d);
      }
      double sum2 = 0;
      for (int j = 0; j < varCat; j++){
        sum2 += lgamma(eps(k, j, d));
      }
      v(k, d) = lgamma(sum1) - sum2;
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec CpostdeltaCalc(arma::vec c, double a, double D){
  arma::vec v(D);
  for (int d = 0; d < D; d++){
    v(d) = lgamma(1 + 2*a) - lgamma(c(d) + a) - lgamma(1 - c(d) + a);
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat logrhonkCalcProfCat(arma::vec Elogpi, arma::mat Elogtheta, arma::vec y, arma::cube carray, arma::mat cmatrix, double K, double D, double N){
  arma::mat v(N, K);
  for (int n = 0; n < N; n++){
    for (int k = 0; k < K; k++){
      double sum1 = 0; // Sum value
      for(int d = 0; d < D; d++){
        sum1 += carray(k, d, n);
      }
      double sum2 = 0; // Sum value
      for(int d = 0; d < D; d++){
        sum2 += cmatrix(n, d);
      }
      v(n, k) = Elogpi(k) + Elogtheta(k, y(n) - 1) + sum1 + sum2;
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube ElogphiCalc(arma::cube eps, double K, double D, double N, double maxNCat, arma::mat X){
  arma::cube v(K, D, N);
  for (int k = 0; k < K; k++){
    for (int d = 0; d < D; d++){
      double sum = 0; // Sum value
      for(int i = 0; i < maxNCat; i++){
        sum += eps(k, i, d);
      }
      for (int n = 0; n < N; n++){
        double j = X(n, d);
        v(k, d, n) = boost::math::digamma(eps(k, j-1, d)) - boost::math::digamma(sum);
      }
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube ElogphiLCalc(arma::cube eps, double K, double D, double maxNCat, arma::vec nCat){
  arma::cube v(K, maxNCat, D);
  for (int d = 0; d < D; d++){
    double varCat = nCat(d);
    for (int k = 0; k < K; k++){
      double sum = 0; // Sum value
      for(int i = 0; i < maxNCat; i++){
        sum += eps(k, i, d);
      }
      for (int l = 0; l < varCat; l++){
        v(k, l, d) = boost::math::digamma(eps(k, l, d)) - boost::math::digamma(sum);
      }
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube lognullphiCalc(arma::mat nullphi, arma::mat X, double K, double D, double N){
  arma::cube v(N, D, K);
  for (int n = 0; n < N; n++){
    for (int k = 0; k < K; k++){
      for(int d = 0; d < D; d++){
        v(n, d, k) = log(nullphi(d, X(n, d) - 1));
      }
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec eta1Calc(arma::cube Elogphi, arma::mat rnk, arma::vec Elogdelta, double K, double D, double N){
  arma::vec v(D);
  for (int d = 0; d < D; d++){
    double sum = 0; // Sum value
    for (int n = 0; n < N; n++){
      for(int k = 0; k < K; k++){
        sum += Elogphi(k, d, n) * rnk(n, k);
      }
    }
    v(d) = exp(sum + Elogdelta(d));
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec eta2Calc(arma::cube lognullphi, arma::mat rnk, arma::vec Elogminusdelta, double K, double D, double N){
  arma::vec v(D);
  for (int d = 0; d < D; d++){
    double sum = 0; // Sum value
    for (int n = 0; n < N; n++){
      for(int k = 0; k < K; k++){
        sum += lognullphi(n, d, k) * rnk(n, k);
      }
    }
    v(d) = exp(sum + Elogminusdelta(d));
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec cCalc(arma::vec eta1, arma::vec eta2, double D){
  arma::vec v(D);
  for (int d = 0; d < D; d++){
    v(d) = eta1(d) / (eta1(d) + eta2(d));
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube epsCalcVarSel(double K, double maxNCat, double D, double N, arma::mat prioreps, arma::mat X, arma::mat rnk, arma::vec c){
  arma::cube v(K, maxNCat, D);
  for (int k = 0; k < K; k++){
    for (int d = 0; d < D; d++){
      for (int l = 0; l < maxNCat; l++){
        double sum = 0; // Sum value
        for(int n = 0; n < N; n++){
          if(X(n, d) == l+1){
            sum += rnk(n, k) * c(d);
          } 
        }
        v(k, l, d) = prioreps(d, l) + sum;
      }
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sumDElogphiCalcVarSel(arma::cube carray, arma::mat cmatrix, double K, double D, double N){
  arma::mat v(N, K);
  for (int n = 0; n < N; n++){
    for (int k = 0; k < K; k++){
      double sum = 0; // Sum value
      for(int d = 0; d < D; d++){
        sum += carray(k, d, n) + cmatrix(n, d);
      }
      v(n, k) = sum;
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube priorepsminusoneCalc(arma::cube prioreps, double K, double D, double maxNCat){
  arma::cube v(K, maxNCat, D);
  for (int l = 0; l < maxNCat; l++){
    for (int k = 0; k < K; k++){
      for(int d = 0; d < D; d++){
        if (prioreps(k, l, d) != 0){
          v(k, l, d) = prioreps(k, l, d) - 1;
        } else{
          v(k, l, d) = 0;
        }
      }
    }
  }
  return v;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube epsminusoneCalc(arma::cube eps, double K, double D, double maxNCat){
  arma::cube v(K, maxNCat, D);
  for (int l = 0; l < maxNCat; l++){
    for (int k = 0; k < K; k++){
      for(int d = 0; d < D; d++){
        if (eps(k, l, d) != 0){
          v(k, l, d) = eps(k, l, d) - 1;
        } else{
          v(k, l, d) = 0;
        }
      }
    }
  }
  return v;
}
