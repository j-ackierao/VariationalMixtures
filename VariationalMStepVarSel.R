#Variational M Step

maxStepVarSel <- function(X, model, prior){
  #Model should contain all current parameters I assume? parameters alpha, epsilon, rnk, labels
  #Need prior in this step as parameters from priors appear in the optimisation equations
  prioralpha <- prior$alpha
  prioreps <- prior$eps
  a <- prior$a
  
  rnk <- model$rnk
  c <- model$c
  nullphi <- model$nullphi
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  #Parameters for pi update - Dirichlet (just need to update parameters)
  alpha <- prioralpha + colSums(rnk)
  
  #Parameters for phi update - Dirichlet
  eps <- array(0, dim = c(K, maxNCat, D))
  for(k in 1:K){
    for(i in 1:D){
      for (l in 1:maxNCat){
        eps[k, l, i] <- prioreps[i, l] + sum((X[,i] == l)*rnk[,k]*c[i])
      }
    }
  }
  
  #First calculate c_i
  Elogphi <- array(0, c(K, D, N))
  #Make an array 
  for (n in 1:N){
    for (i in 1:D){
      for (k in 1:K){
        Elogphi[k,i, n] <- digamma(eps[k, X[n,i], i]) - digamma(sum(eps[k, 1:maxNCat, i])) #last term sum over Li variables
      }
    }
  }
  
  lognullphi <- array(0, c(N, D, K)) #phi_0ixni
  for (n in 1:N){
    for (i in 1:D){
      for (k in 1:K){
        lognullphi[n, i, k] <- log(nullphi[i, X[n, i]])
      }
    }
  }
  
  Elogdelta <- digamma(c + a) - digamma(2*a + 1)
  Elogminusdelta <- digamma(1 - c + a) - digamma(2*a + 1)
  
  logeta1 <- c()
  for (i in 1:D){
    logeta1[i] <- sum(t(Elogphi[1:K, i, 1:N]) * rnk) + Elogdelta[i]
  }
  
  logeta2 <- c()
  for (i in 1:D){
    logeta2[i] <- sum(lognullphi[1:N, i, 1:K] * rnk) + Elogminusdelta[i]
  }
  
  eta1 <- exp(logeta1)
  eta2 <- exp(logeta2)
  
  c <- c()
  for (i in 1:D){
    c[i] <- eta1[i] / (eta1[i] + eta2[i])
  }
  c[is.na(c)] <- 0
  
  model$alpha <- alpha #update alpha* in model
  model$eps <- eps #update epsilon* in model
  model$c <- c #update c in model
  return(model)
}
