#Variational M Step - profile regression w/categorical response
library(matrixStats)
maxStepProfCat <- function(X, y, model, prior){
  #Model should contain all current parameters
  #Need prior in this step as parameters from priors appear in the optimisation equations
  prioralpha <- prior$alpha
  prioreps <- prior$eps
  priorbeta <- prior$beta
  a <- prior$a
  
  rnk <- model$rnk
  c <- model$c
  nullphi <- model$nullphi
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  J = dim(model$beta)[2]
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  #Parameters for pi update - Dirichlet
  alpha <- prioralpha + colSums(rnk)
  
  #Parameters for phi and theta update - Dirichlet
  eps <- epsCalcVarSel(K, maxNCat, D, N, prioreps, X, rnk, c)
  beta <- betaCalc(priorbeta, y, K, J, N, rnk)
  
  #First calculate c_i
  Elogphi <- ElogphiCalc(eps, K, D, N, maxNCat, X)
  lognullphi <- lognullphiCalc(nullphi, X, K, D, N) #phi_0ixni
  
  Elogdelta <- digamma(c + a) - digamma(2*a + 1)
  Elogminusdelta <- digamma(1 - c + a) - digamma(2*a + 1)
  
  logeta1 <- as.vector(logeta1Calc(Elogphi, rnk, Elogdelta, K, D, N))
  logeta2 <- as.vector(logeta2Calc(lognullphi, rnk, Elogminusdelta, K, D, N))
  logetas <- matrix(c(logeta1, logeta2), nrow = D, ncol = 2)
  
  clse <- rowLogSumExps(logetas)
  c <- as.vector(cCalc(logeta1, clse, D))
  
  model$alpha <- alpha #update alpha* in model
  model$eps <- eps #update epsilon* in model
  model$c <- c #update c in model
  model$beta <- beta #update beta* in model
  return(model)
}
