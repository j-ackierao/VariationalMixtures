#Variational E step - profile regression w/ categorical response
library(matrixStats)
expectStepProfCat <- function(X, y, model){
  #Model should contain all current parameters: parameters alpha, epsilon, c, labels
  #Add parameter rnk (responsibilities) in this step; this is the first step before maximisation
  
  alpha <- model$alpha
  eps <- model$eps
  c <- model$c
  nullphi <- model$nullphi
  beta <- model$beta
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  J = dim(beta)[2]
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #k dim vector, expectation of log pi k
  Elogphi <- ElogphiCalc(eps, K, D, N, maxNCat, X)
  Elogtheta <- ElogthetaCalcCat(beta, K, J)
  
  carray <- replicate(N, matrix(rep(c, K), ncol = D, byrow = TRUE), simplify="array") * Elogphi
  #array of c_i * Elogphi
  cmatrix <- cmatrixCalc(nullphi, X, c, N, D)
  
  logrhonk <- logrhonkCalcProfCat(Elogpi, Elogtheta, y, carray, cmatrix, K, D, N)
  lse <- rowLogSumExps(logrhonk)
  rnk <- rnkCalc(logrhonk, lse, N, K)
  
  labels <- apply(rnk, 1, which.max) #k with the highest responsibility is the cluster that zn is assigned to
  
  model$rnk <- rnk #update responsibilities in model
  model$labels <- labels #update labels in model
  
  return(model)
  
  
}
