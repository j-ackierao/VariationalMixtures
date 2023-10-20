#ELBO calculation function

ELBOCalc <- function(X, model, prior){
  
  prior2 = list(alpha = rep(prior$alpha, K),
               eps = EPSreshape[rep(1,K),,])
  prioralpha <- prior2$alpha
  prioreps <- prior2$eps
  
  alpha <- model$alpha
  eps <- model$eps
  rnk <- model$rnk
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #Taken from E step
  Elogphi <- array(0, c(K, D, N))
  #Make an array 
  for (n in 1:N){
    for (i in 1:D){
      for (k in 1:K){
        Elogphi[k,i, n] <- digamma(eps[k, X[n,i], i]) - digamma(sum(eps[k, 1:maxNCat, i])) #last term sum over Li variables
      }
    }
  }
  
  ElogphiL <- array(0, c(K, maxNCat, D))
  #Make an array 
  for (i in 1:D){
    varCat <- nCat[i]
    for (l in 1:varCat){
      for (k in 1:K){
        ElogphiL[k,l, i] <- digamma(eps[k, l, i]) - digamma(sum(eps[k, 1:maxNCat, i])) #last term sum over Li variables
      }
    }
  }
  
  #(log) normalising constants of Dirichlet
  Cprioralpha <- lgamma(sum(prioralpha)) - sum(lgamma(prioralpha))
  Cpostalpha <- lgamma(sum(alpha)) - sum(lgamma(alpha))
  #K X D matrix for epsilon constants as thats how many phi priors (and posteriors) we have, all dim L (no ofcategories)
  #NEED TO MAKE SURE ZEROES ARE NOW IGNORED
  Cprioreps <- array(0, c(K, D))
  for (k in 1:K){
    for (i in 1:D){
      varCat <- nCat[i]
      Cprioreps[k, i] <- lgamma(sum(prioreps[k, 1:varCat, i])) - sum(lgamma(prioreps[k, 1:varCat, i]))
    }
  }
  Cposteps <- array(0, c(K, D))
  for (k in 1:K){
    for (i in 1:D){
      varCat <- nCat[i]
      Cposteps[k, i] <- lgamma(sum(eps[k, 1:varCat, i])) - sum(lgamma(eps[k, 1:varCat, i]))
    }
  }
    
  #Calculations
  
  sumDElogphi <- array(0, c(N, K)) #nxK matrix of sum of E(log phi) over i = 1, ..., D, used in calculation of first expression
  for (n in 1:N){
    for (k in 1:K){
      sumDElogphi[n, k] <- sum(Elogphi[k, 1:D, n])
    }
  }
  
  #Matrix of epsilon parameters -1, where all 0's remain 0 (as these are unused parameters)
  priorepsminusone <- array(0, c(K, maxNCat, D))
  for (k in 1:K){
    for (l in 1:maxNCat){
      for (i in 1:D){
        if (prioreps[k, l, i] != 0){
          priorepsminusone[k,l,i] = prioreps[k,l,i] - 1
        }
        else{
          priorepsminusone[k,l,i] == 0
        }
      }
    }
  }
  
  epsminusone <- array(0, c(K, maxNCat, D))
  for (k in 1:K){
    for (l in 1:maxNCat){
      for (i in 1:D){
        if (eps[k, l, i] != 0){
          epsminusone[k,l,i] = eps[k,l,i] - 1
        }
        else{
          epsminusone[k,l,i] == 0
        }
      }
    }
  }
  
  matExp1 <- rnk * sumDElogphi
  
  Exp1 <- sum(matExp1) #E(logp(X|Z,phi))
  
  matExp2 <- rnk * matrix(rep(Elogpi,N), ncol = K, byrow = TRUE)
  
  Exp2 <- sum(matExp2) #E(logp(Z|pi))
  
  Exp3 <- sum((prioralpha - 1)*Elogpi) + Cprioralpha #E(logp(pi)), sum will sum over all k
  
  Exp4 <- sum((priorepsminusone)*ElogphiL) + sum(Cprioreps) #E(logp(phi)) I think this is correct but double check ElogphiL
  
  Exp5 <- sum(rnk * log(rnk)) #E(q(Z)) this should be correct - element-wise multiplication
  
  Exp6 <- sum((alpha - 1)*Elogpi) + Cpostalpha #E(logq(pi))
  
  Exp7 <- sum((epsminusone)*ElogphiL) + sum(Cposteps) #Elogq(phi)
  
  ELBO <- Exp1 + Exp2 + Exp3 + Exp4 - Exp5 - Exp6 - Exp7
  
}
