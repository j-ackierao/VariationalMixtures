#Variational E step

expectStep <- function(X, model){
  #Model should contain all current parameters I assume? parameters alpha, epsilon, labels
  #Add parameter rnk (responsibilities) in this step; this is the first step before maximisation
  
  alpha <- model$alpha
  eps <- model$eps
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #k dim vector, expectation of log pi k
  
  #This is probably super inefficient?
  Elogphi <- array(0, c(K, D, N))
  #Make an array 
  for (n in 1:N){
    for (i in 1:D){
      for (k in 1:K){
        Elogphi[k,i, n] <- digamma(eps[k, X[n,i], i]) - digamma(sum(eps[k, 1:maxNCat, i])) #last term sum over Li variables
      }
    }
  }

  logrhonk <- array(0, c(N, K)) #calculate log rho_nk
  for (n in 1:N){
    for (k in 1:K){
      logrhonk[n, k] <- Elogpi[k] + sum(Elogphi[k, 1:D, n])
    }
  }
  
  rhonk <- exp(logrhonk) #calculate rho_nk
  
  rnk <- array(0, c(N, K)) #calculate responsibilities
  for (n in 1:N){
    for (k in 1:K){
      rnk[n, k] <- rhonk[n, k] / sum(rhonk[n, 1:K])
    }
  }
  
  #Check this bit! is the k with the highest responsibility the cluster that zn is assigned to?
  
  labels <- apply(rnk, 1, which.max)
  
  model$rnk <- rnk #update responsibilities in model
  model$labels <- labels #update labels in model

  return(model)

}