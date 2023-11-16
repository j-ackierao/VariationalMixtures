#ELBO calculation function - variable selection

ELBOCalcVarSel <- function(X, model, prior){
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  
  EPSreshape = t(prior$eps) 
  dim(EPSreshape) = c(1,maxNCat,D)
  prior2 = list(alpha = rep(prior$alpha, K),
                eps = EPSreshape[rep(1,K),,])
  prioralpha <- prior2$alpha
  prioreps <- prior2$eps
  a <- prior$a
  
  alpha <- model$alpha
  eps <- model$eps
  rnk <- model$rnk
  c <- model$c
  nullphi <- model$nullphi
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #Taken from E step
  Elogphi <- array(0, c(K, D, N)) #sum over all xni
  #Make an array 
  for (n in 1:N){
    for (i in 1:D){
      for (k in 1:K){
        Elogphi[k,i, n] <- digamma(eps[k, X[n,i], i]) - digamma(sum(eps[k, 1:maxNCat, i])) #last term sum over Li variables
      }
    }
  }
  
  ElogphiL <- array(0, c(K, maxNCat, D)) #sum over all l
  #Make an array 
  for (i in 1:D){
    varCat <- nCat[i]
    for (l in 1:varCat){
      for (k in 1:K){
        ElogphiL[k,l, i] <- digamma(eps[k, l, i]) - digamma(sum(eps[k, 1:maxNCat, i])) #last term sum over Li variables
      }
    }
  }
  
  Elogdelta <- digamma(c + a) - digamma(2*a + 1)
  Elogminusdelta <- digamma(1 - c + a) - digamma(2*a + 1)
  
  #(log) normalising constants of Dirichlet - STAYS THE SAME COMPARED TO NO VARIABLE SELECTION
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
  Cpriordelta <- lgamma(a + a) - 2 * lgamma(a)
  Cpostdelta <- c()
  for (i in 1:D){
    Cpostdelta[i] <- lgamma(1 + 2*a) - lgamma(c[i] + a) - lgamma(1 - c[i] + a)
  }
  
  #Calculations
  
  carray <- replicate(N, matrix(rep(c, K), ncol = D, byrow = TRUE), simplify="array") * Elogphi
  #array of c_i * Elogphi
  cmatrix <- array(0, c(N, D)) #c_i * phi_0ixni
  for (n in 1:N){
    for (i in 1:D){
      cmatrix[n, i] <- (1 - c[i]) * log(nullphi[i, X[n, i]])
    }
  }
  
  sumDElogphi <- array(0, c(N, K)) #nxK matrix used in calculation of first expression: the sum from 1:D
  for (n in 1:N){
    for (k in 1:K){
      sumDElogphi[n, k] <- sum(carray[k, 1:D, n]) + sum(cmatrix[n, 1:D])
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
  
  Exp1 <- sum(matExp1) #E(logp(X|Z,phi)) #DONE
  
  matExp2 <- rnk * matrix(rep(Elogpi,N), ncol = K, byrow = TRUE)
  
  Exp2 <- sum(matExp2) #E(logp(Z|pi)) #STAYS THE SAME
  
  Exp3 <- sum((prioralpha - 1)*Elogpi) + Cprioralpha #E(logp(pi)) STAYS THE SAME
  
  Exp4 <- sum((priorepsminusone)*ElogphiL) + sum(Cprioreps) #E(logp(phi)) STAYS THE SAME
  
  Exp5 <- sum((c * Elogdelta) + (1-c)*Elogminusdelta) #E(logp(gamma|delta)) DONE
  
  Exp6 <- sum((a - 1) * Elogdelta + (a - 1) * Elogminusdelta + Cpriordelta) #E(logp(delta)) DONE
  
  Exp7 <- sum(rnk * log(rnk)) #E(q(Z)) #STAYS THE SAME
  
  Exp8 <- sum((alpha - 1)*Elogpi) + Cpostalpha #E(logq(pi)) #STAYS THE SAME
  
  Exp9 <- sum((epsminusone)*ElogphiL) + sum(Cposteps) #Elogq(phi) #STAYS THE SAME
  
  matExp10 <- (c * log(c)) + ((1-c) * log(1-c))
  matExp10[is.na(matExp10)] <- 0
  
  Exp10 <- sum(matExp10) #Elogq(gamma) #unsure on this...
  
  Exp11 <- sum((c + a - 1) * Elogdelta + (a - c) * Elogminusdelta + Cpostdelta) #Elogq(delta) DONE
  
  ELBO <- Exp1 + Exp2 + Exp3 + Exp4 + Exp5 + Exp6 - Exp7 - Exp8 - Exp9 - Exp10 - Exp11
  
}
