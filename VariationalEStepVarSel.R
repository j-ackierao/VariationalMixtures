#Variational E step - model selection

expectStepVarSel <- function(X, model){
  #Model should contain all current parameters: parameters alpha, epsilon, c, labels
  #Add parameter rnk (responsibilities) in this step; this is the first step before maximisation
  
  alpha <- model$alpha
  eps <- model$eps
  c <- model$c
  nullphi <- model$nullphi
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
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
  
  carray <- replicate(N, matrix(rep(c, K), ncol = D, byrow = TRUE), simplify="array") * Elogphi
  #array of c_i * Elogphi
  cmatrix <- array(0, c(N, D)) #1 - c_i * phi_0ixni
  for (n in 1:N){
    for (i in 1:D){
      cmatrix[n, i] <- (1 - c[i]) * log(nullphi[i, X[n, i]])
    }
  }
  
  logrhonk <- array(0, c(N, K)) #calculate log rho_nk
  for (n in 1:N){
    for (k in 1:K){
      logrhonk[n, k] <- Elogpi[k] + sum(carray[k, 1:D, n]) + sum(cmatrix[n, 1:D])
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
