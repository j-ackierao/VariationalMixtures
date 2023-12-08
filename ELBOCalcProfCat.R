#ELBO calculation function - profile regression
library(Rcpp)

ELBOCalcProfCat <- function(X, y, model, prior){
  N = dim(X)[1]
  D = dim(X)[2]
  J = dim(model$beta)[2]
  K = length(model$alpha)
  
  EPSreshape = t(prior$eps) 
  dim(EPSreshape) = c(1,maxNCat,D)
  prior2 = list(alpha = rep(prior$alpha, K),
                eps = EPSreshape[rep(1,K),,])
  prioralpha <- prior2$alpha
  prioreps <- prior2$eps
  a <- prior$a
  priorbeta <- prior$beta
  
  alpha <- model$alpha
  eps <- model$eps
  rnk <- model$rnk
  c <- model$c
  nullphi <- model$nullphi
  beta <- model$beta
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #Taken from E step
  Elogphi <- ElogphiCalc(eps, K, D, N, maxNCat, X)
  
  ElogphiL <- ElogphiLCalc(eps, K, D, maxNCat, nCat)
  
  Elogdelta <- digamma(c + a) - digamma(2*a + 1)
  Elogminusdelta <- digamma(1 - c + a) - digamma(2*a + 1)
  Elogtheta <- ElogthetaCalcCat(beta, K, J)
  
  #(log) normalising constants of Dirichlet - STAYS THE SAME COMPARED TO NO VARIABLE SELECTION
  Cprioralpha <- lgamma(sum(prioralpha)) - sum(lgamma(prioralpha))
  Cpostalpha <- lgamma(sum(alpha)) - sum(lgamma(alpha))
  #K X D matrix for epsilon constants as thats how many phi priors (and posteriors) we have, all dim L (no ofcategories)
  #NEED TO MAKE SURE ZEROES ARE NOW IGNORED
  Cprioreps <- CpriorepsCalc(prioreps, K, D, nCat)
  Cposteps <- CpostepsCalc(eps, K, D, nCat)
  Cpriordelta <- lgamma(a + a) - 2 * lgamma(a)
  Cpostdelta <- CpostdeltaCalc(c, a, D)
  Cpriortheta <- CpriorbetaCalc(priorbeta, K, J)
  Cposttheta <- CpostbetaCalc(beta, K, J)
  
  #Calculations
  
  carray <- replicate(N, matrix(rep(c, K), ncol = D, byrow = TRUE), simplify="array") * Elogphi
  #array of c_i * Elogphi
  cmatrix <- cmatrixCalc(nullphi, X, c, N, D) #c_i * phi_0ixni
  
  sumDElogphi <- sumDElogphiCalcVarSel(carray, cmatrix, K, D, N)
  
  #Matrix of epsilon parameters -1, where all 0's remain 0 (as these are unused parameters)
  priorepsminusone <- priorepsminusoneCalc(prioreps, K, D, maxNCat)
  epsminusone <- epsminusoneCalc(eps, K, D, maxNCat)
  
  resptheta <- respthetaCalc(Elogtheta, rnk, y, N, K)
  matExp1 <- rnk * sumDElogphi
  
  Exp1 <- sum(resptheta) + sum(matExp1) #E(logp(X, Y|Z,phi, theta, gamma)) #DONE
  
  matExp2 <- rnk * matrix(rep(Elogpi,N), ncol = K, byrow = TRUE)
  
  Exp2 <- sum(matExp2) #E(logp(Z|pi)) #STAYS THE SAME
  
  Exp3 <- sum((prioralpha - 1)*Elogpi) + Cprioralpha #E(logp(pi)) STAYS THE SAME
  
  Exp4 <- sum((priorepsminusone)*ElogphiL) + sum(Cprioreps) #E(logp(phi)) STAYS THE SAME
  
  Exp5 <- sum((c * Elogdelta) + (1-c)*Elogminusdelta) #E(logp(gamma|delta)) STAYS THE SAME
  
  Exp6 <- sum((a - 1) * Elogdelta + (a - 1) * Elogminusdelta + Cpriordelta) #E(logp(delta)) STAYS THE SAME
  
  Exp7 <- sum((priorbeta - 1) * Elogtheta) + sum(Cpriortheta) #Elogptheta DONE
  
  Exp8 <- sum(rnk * log(rnk)) #E(q(Z)) #STAYS THE SAME
  
  Exp9 <- sum((alpha - 1)*Elogpi) + Cpostalpha #E(logq(pi)) #STAYS THE SAME
  
  Exp10 <- sum((epsminusone)*ElogphiL) + sum(Cposteps) #Elogq(phi) #STAYS THE SAME
  
  matExp11 <- (c * log(c)) + ((1-c) * log(1-c))
  matExp11[is.na(matExp11)] <- 0
  
  Exp11 <- sum(matExp11) #Elogq(gamma) #unsure on this...
  
  Exp12 <- sum((c + a - 1) * Elogdelta + (a - c) * Elogminusdelta + Cpostdelta) #Elogq(delta) STAYS THE SAME
  
  Exp13 <- sum((beta - 1) * Elogtheta) + sum(Cposttheta) #E(logq(theta)) DONE
  
  ELBO <- Exp1 + Exp2 + Exp3 + Exp4 + Exp5 + Exp6 + Exp7 - Exp8 - Exp9 - Exp10 - Exp11 - Exp12 - Exp13
  
}
