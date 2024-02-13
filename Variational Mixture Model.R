library(klaR)
library(tidyverse)
library(matrixStats)
source('VariationalEStep.R')
source('VariationalMStep.R')
source('ELBOCalc.R')

library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp('MixModel.cpp')
Rcpp::sourceCpp('MixModelArma.cpp')

#Variational mixtures for categorical distributions
check_convergence<- function(ELBO, iter, maxiter, tol){
  if (iter > 1 && abs(ELBO[iter * 2] - ELBO[iter * 2-1]) < tol && abs(ELBO[iter * 2-1] - ELBO[iter * 2-2] < tol) && 
      abs(ELBO[iter * 2-2] - ELBO[iter * 2-3] < tol)){
    print(paste("Stopped after iteration ",iter)) #make sure the last 3 ELBOs close to each other
    return(TRUE)
  }
  if (iter == maxiter){
    print(paste("Not converged after maximum number of iterations"))
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

#data = 'X' observed variables
#K = number of clusters
#alpha = prior parameters for clusters
#maxiter = maximum number of iterations
#tol = convergence criteria

mixturemodel <- function(data, K, alpha, maxiter, tol){
  ################################################
  #Process dataset
  
  #First make sure X (data) is a data frame and make every column into factors
  X <- as.data.frame(data)
  X[colnames(X)] <- lapply(X[colnames(X)], as.factor)
  
  #Data frame mapping factor labels to numbers
  factor_labels = data.frame()
  for (i in 1:length(colnames(X))){
    factorlist <- data.frame(factors = unique(X[,i]), value = as.numeric(unique(X[,i])))
    factordict <- cbind(data.frame(variable = colnames(X)[i]), factorlist)
    factor_labels <- rbind(factor_labels, factordict)
  }
  
  #Check for any completely empty factors as these may cause problems
  categories <- lapply(1:ncol(X), function(j)levels(X[, j])) 
  cat_lengths <- sapply(categories, length)
  if(any(cat_lengths == 0)) {
    stop("Column(s) ", paste0(colnames(X)[cat_lengths == 0], collapse = ","), " is empty")
  }
  #Currently work under assumption there are no missing values
  
  #Create numeric matrix for data
  X <- data.matrix(X)
  #note that now eg. binary factors are now 1 and 2 instead of 0 and 1
  
  ##########################################
  
  N = dim(X)[1]
  D = dim(X)[2]
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  ELBO = Cl = rep(-Inf,maxiter*2) #record ELBO and no of clusters at each iteration
  
  ########################################
  #Define priors - given input alpha, all variables have a symmetric Dirichlet prior with parameter alpha
  
  prior = list(alpha = alpha) #prior for clustering pi
  prior$eps = matrix(0, D, maxNCat) #prior for clustering phi
  for(d in 1:D){
    prior$eps [d,1:nCat[d]] <- 1/nCat[d]
  }
  
  #Initialising cluster labels
  #cluster = sample(1:K, N, replace=T) random assignment!
  clusterInit <- klaR::kmodes(X, modes = K)$cluster #k modes: analogous to k means
  
  #######################################
  
  EPSreshape = t(prior$eps) 
  dim(EPSreshape) = c(1,maxNCat,D) #sets dimension of matrix, D arrays (variable) where each one has 1 row and 
                                  #maxNCat columns (parameter for each category in the variable, all unused are 0)
  model = list(alpha = rep(prior$alpha, K),
               eps = EPSreshape[rep(1,K),,],
               labels = clusterInit) #current parameters for alpha, epsilon: makes D arrays
            #ie. array for each variable, each array has a row for each cluster and a column for each category
            #each array represents phi_ki k*max(Li) number of columns and rows, zeroes to stand for categories unused
            #model alpha and model epsilon is therefore pi parameters and phi parameters

  #######################################
  model$eps <- firstepsCalc(K, maxNCat, D, N, prior$eps, X, clusterInit)
  #Update the epsilons based on the initial cluster assignment, based on how many observations are in each 
  #cluster with a given variable value - breaks symmetries between clusters as all priors are symmetric 
  #Similar to our usual phi (epsilon) update based on initial cluster assignments (rnk = 1 if n is in cluster k)
  
  ########################################

  for (iter in 1:maxiter){
    print(paste("Iteration number ",iter))
    model = expectStep(X, model) #Expectation step
    .GlobalEnv$maxNCat <- maxNCat
    ELBO[iter * 2-1] = ELBOCalc(X, model, prior) #ELBO
    print(ELBO[iter * 2-1])

    Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
    model = maxStep(X, model, prior) #Maximisation step
    ELBO[iter * 2] = ELBOCalc(X, model, prior) #ELBO
    print(ELBO[iter * 2])

    Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
    
    if(check_convergence(ELBO, iter, maxiter, tol)) break
    
  }
  
  output <- list(ELBO = ELBO[1:(iter*2)], Cl = Cl[1:iter], model = model, factor_labels = factor_labels)
  return(output)
  
}


