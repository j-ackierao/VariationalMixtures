#Variational mixtures for categorical distributions

mixturemodel <- function(X, K, alpha, maxiter){
  
  ################################################
  #Process dataset - use inspiration from vimix and mixdir codes
  
  #First make sure X is a data frame and make every column into factors
  X <- as.data.frame(X)
  X[colnames(X)] <- lapply(X[colnames(X)], as.factor)
  
  
  #Data frame mapping factor labels to numbers
  factor_labels = data.frame()
  for (i in 1:length(colnames(X))){
    factorlist <- data.frame(factors = unique(X[,i]), value = as.numeric(unique(X[,i])))
    factordict <- cbind(variable = colnames(X)[1], factorlist)
    factor_labels <- rbind(factor_labels, factordict)
  }
  
  #Check for any completely empty factors as these may cause problems
  categories <- lapply(1:ncol(X), function(j)levels(X[, j])) #list out categories
  
  cat_lengths <- sapply(categories, length)
  if(any(cat_lengths == 0)) {
    stop("Column(s) ", paste0(colnames(X)[cat_lengths == 0], collapse = ","), " is empty (i.e. all values are NA). Please fix this.")
  }
  #Currently work under assumption there are no missing values? unsure how to deal with this now lol
  
  #Create numeric matrix for data
  X <- data.matrix(X)
  #note that now eg. binary factors are now 1 and 2 instead of 0 and 1
  
  
  ##########################################
  
  N = dim(X)[1]
  D = dim(X)[2]
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  L = Cl = rep(-Inf,maxiter*2) #prepare for L (ELBO) and number of clusters at each iteration to be recorded
  
  ########################################
  #Define priors
  #Assume there is an input alpha and all variables have a symmetric Dirichlet prior with parameters alpha/Li where
  #there are Li different potential categories for variable i
  
  #Prior for parameters in each cluster
  
  prior = list(alpha = 1/K) #default value of alpha - prior for clustering pi
  prior$eps = matrix(0, D, maxNCat) #prior for clustering phi. Eps stands for epsilon 
  for(d in 1:D){
    prior$eps [d,1:nCat[d]] <- 1/nCat[d]
  }
  
  #Initialising
  #cluster = sample(1:K, N, replace=T) random assignment!
  clusterInit <- klaR::kmodes(X, modes = K)$cluster
  
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
  model$eps <- model$eps
  for(i in 1:D){
    for(j in 1:nCat[i]){
      for(k in 1:K){
        model$eps[k,j,i] <- prior$eps[i,j] + sum((X[,i]==j)*(clusterInit==k))
      }
    }
  }
  #I think this tries to update the epsilons based on the initial cluster assignment, based on how many observations are in each 
  #cluster with a given variable value - breaks symmetries between clusters as all priors are symmetric 
  #Look more into why this is the calculation done
  
  ########################################
  
  #Define ELBO calculation
  
  #Initialise model
  
  #Set function to iterate between E and M and calculate ELBO between each iteration (after E and after M)
  #Record these in a little list
  #Also record number of clusters each time
  
  check_convergence<- function(L){L + 1
    }
  
  for (iter in 1:maxiter){
    print(paste("Iteration number ",iter))
    model = expectStep(X, model) #Expectation step
    ELBO[iter * 2-1] = ELBOCalc(X, model, prior) #ELBO
    print(ELBO[iter * 2-1])
    if(iter > 1 && ELBO[iter * 2-1] < ELBO[iter * 2-2]){
      print("SOMETHING WENT WRONG")
    }
    Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
    model = maxStep(X, model, prior) #Maximisation step
    ELBO[iter * 2] = ELBOCalc(X, model, prior) #ELBO
    print(ELBO[iter * 2])
    if(ELBO[iter * 2] < ELBO[iter * 2-1]){
      print("SOMETHING WENT WRONG")
    }
    Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
    
    #if(check_convergence(L)) break
  }
  
  #Define how to check for convergence
  
}

