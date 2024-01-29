#ELBO calculation function

ELBOCalc <- function(X, model, prior){
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  
  EPSreshape = t(prior$eps) 
  dim(EPSreshape) = c(1,maxNCat,D)
  prior2 = list(alpha = rep(prior$alpha, K),
               eps = EPSreshape[rep(1,K),,])
  prioralpha <- prior2$alpha
  prioreps <- prior2$eps
  
  alpha <- model$alpha
  eps <- model$eps
  rnk <- model$rnk
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #Taken from E step
  Elogphi <- ElogphiCalc(eps, K, D, N, maxNCat, X)
  ElogphiL <- ElogphiLCalc(eps, K, D, maxNCat, nCat)
  
  #(log) normalising constants of Dirichlet
  Cprioralpha <- lgamma(sum(prioralpha)) - sum(lgamma(prioralpha))
  Cpostalpha <- lgamma(sum(alpha)) - sum(lgamma(alpha))
  #K X D matrix for epsilon constants as thats how many phi priors (and posteriors) we have, all dim L (no ofcategories)
  #NEED TO MAKE SURE ZEROES ARE NOW IGNORED
  Cprioreps <- CpriorepsCalc(prioreps, K, D, nCat)
  Cposteps <- CpostepsCalc(eps, K, D, nCat)
    
  #Calculations
  
  sumDElogphi <- sumDElogphiCalc(Elogphi, K, D, N) #nxK matrix of sum of E(log phi) over i = 1, ..., D, used in calculation of first expression
  
  #Matrix of epsilon parameters -1, where all 0's remain 0 (as these are unused parameters)
  priorepsminusone <- priorepsminusoneCalc(prioreps, K, D, maxNCat)
  epsminusone <- epsminusoneCalc(eps, K, D, maxNCat)
  
  matExp1 <- rnk * sumDElogphi
  
  Exp1 <- sum(matExp1) #E(logp(X|Z,phi))
  
  matExp2 <- rnk * matrix(rep(Elogpi,N), ncol = K, byrow = TRUE)
  
  Exp2 <- sum(matExp2) #E(logp(Z|pi))
  
  Exp3 <- sum((prioralpha - 1)*Elogpi) + Cprioralpha #E(logp(pi)), sum will sum over all k. CHECKED
  
  Exp4 <- sum((priorepsminusone)*ElogphiL) + sum(Cprioreps) #E(logp(phi)) I think this is correct but double check ElogphiL
  
  logrnk <- log(rnk)
  logrnk[logrnk == -Inf] <- 0
  Exp5 <- sum(rnk * logrnk) #E(q(Z)) this should be correct - element-wise multiplication
  
  Exp6 <- sum((alpha - 1)*Elogpi) + Cpostalpha #E(logq(pi))
  
  Exp7 <- sum((epsminusone)*ElogphiL) + sum(Cposteps) #Elogq(phi)
  
  ELBO <- Exp1 + Exp2 + Exp3 + Exp4 - Exp5 - Exp6 - Exp7
  
}
