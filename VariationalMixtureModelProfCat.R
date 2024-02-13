library(klaR)
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
source('VariationalEStepProf.R')
source('VariationalMStepProf.R')
source('ELBOCalcProf.R')

library(Rcpp)
Rcpp::sourceCpp('ProfCat.cpp')
Rcpp::sourceCpp('ProfCatArma.cpp')

#Variational mixtures for categorical distributions with variable selection - outcome included
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
#y = outcome (nx1 vector) 
#K = number of clusters
#alpha = prior parameters for clusters
#a = prior parameter for delta
#beta = prior parameters for theta. Can be 1 dim or 2 dim (beta = (beta_0, beta_1)
#maxiter = maximum number of iterations
#tol = convergence criteria

mixturemodelprofCat <- function(data, outcome, K, alpha, a, maxiter, tol){
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
  
  y <- as.data.frame(outcome)
  y <- as.numeric(y[,1])
  ##########################################
  
  N = dim(X)[1]
  D = dim(X)[2]
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  ELBO = Cl = rep(-Inf,maxiter*2) #record ELBO and no of clusters at each iteration
  J <- length(unique(y))
  
  ########################################
  #Define priors - given input alpha, all variables have a symmetric Dirichlet prior with parameter alpha
  
  prior = list(alpha = alpha) #prior for clustering pi
  prior$eps = matrix(0, D, maxNCat) #prior for clustering phi
  for(d in 1:D){
    prior$eps [d,1:nCat[d]] <- 1/nCat[d]
  }
  prior$a = a
  prior$beta = rep(1/J,each=J)
  #could change this to be unequal across different categories
  
  #Initialising
  #cluster = sample(1:K, N, replace=T) random assignment!
  clusterInit <- klaR::kmodes(X, modes = K)$cluster
  
  #######################################
  
  EPSreshape = t(prior$eps) 
  dim(EPSreshape) = c(1,maxNCat,D) #sets dimension of matrix, D arrays (variable) where each one has 1 row and 
  #maxNCat columns (parameter for each category in the variable, all unused are 0)
  model = list(alpha = rep(prior$alpha, K),
               eps = EPSreshape[rep(1,K),,],
               c = rep(1, D), #initialise c_i - all variables initially included
               beta = matrix(rep(prior$beta,K),ncol=J,byrow=TRUE),
               labels = clusterInit) #current parameters for alpha, epsilon: makes D arrays
  #ie. array for each variable, each array has a row for each cluster and a column for each category
  #each array represents phi_ki k*max(Li) number of columns and rows, zeroes to stand for categories unused
  #model alpha and model epsilon is therefore pi parameters and phi parameters
  
  #######################################
  
  #Set null phi to the rate of the parameter value in the dataset - max likelihood
  nullphi <- array(0, dim = c(D, maxNCat))
  for (i in 1:D){
    for (j in 1:nCat[i]){
      nullphi[i, j] <- sum((X[,i] == j)) / N
    }
  }
  model$nullphi <- nullphi
  
  #######################################
  model$eps <- firstepsCalc(K, maxNCat, D, N, prior$eps, X, clusterInit)
  #Update the epsilons based on the initial cluster assignment, based on how many observations are in each 
  #cluster with a given variable value - breaks symmetries between clusters as all priors are symmetric 
  #Similar to our usual phi (epsilon) update based on initial cluster assignments (rnk = 1 if n is in cluster k)
  
  #Also update betas based on initial cluster assignment
  for(j in 1:J){
    for(k in 1:K){
      model$beta[k, j] <- prior$beta[j] + sum((y==j)*(clusterInit==k))
    }
  }
  
  ########################################
  
  
  for (iter in 1:maxiter){
    print(paste("Iteration number ",iter))
    model = expectStepProfCat(X, y, model) #Expectation step
    .GlobalEnv$maxNCat <- maxNCat
    ELBO[iter * 2-1] = ELBOCalcProfCat(X, y, model, prior) #ELBO
    print(ELBO[iter * 2-1])
    
    Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
    model = maxStepProfCat(X, y, model, prior) #Maximisation step
    ELBO[iter * 2] = ELBOCalcProfCat(X, y, model, prior) #ELBO
    print(ELBO[iter * 2])
    
    Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
    
    if(check_convergence(ELBO, iter, maxiter, tol)) break
    
  }
  
  output <- list(ELBO = ELBO[1:(iter*2)], Cl = Cl[1:iter], model = model, factor_labels = factor_labels)
  return(output)
  
}
