#Variational M Step

maxStep <- function(X, model, prior){
  #Model should contain all current parameters I assume? parameters alpha, epsilon, rnk, labels
  #Need prior in this step as parameters from priors appear in the optimisation equations
  prioralpha <- prior$alpha
  prioreps <- prior$eps
  rnk <- model$rnk
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  #Parameters for pi update - Dirichlet (just need to update parameters, don't need full dist at this point)
  alpha <- prioralpha + colSums(rnk)
  
  #Parameters for phi update - Dirichlet
  eps <- array(0, dim = c(K, maxNCat, D))
  for(k in 1:K){
    for(i in 1:D){
      for (l in 1:maxNCat){
        eps[k, l, i] <- prioreps[i, l] + sum((X[,i] == l)*rnk[,k])
      }
    }
  }
  
  model$alpha <- alpha
  model$eps <- eps
  return(model)
}
