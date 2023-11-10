#Generate multiview sample data
source("GenerateSampleData.R")

#Example settings for parameters in function
n    <- 1000 # The overall population size
views <- 2 #The number of clustering strutures present
K    <- c(4, 5)    # The number of clusters/subpopulations
w    <- list(c(0.1, 0.2, 0.3, 0.4), c(0.1, 0.2, 0.2, 0.2, 0.3)) # The mixture weights (proportion of population in each cluster)
#w    <- w/sum(w)   # This ensures that the elements of w sum to 1 
p    <- c(15, 15) # The number of variables. dim(p) = views 
Irrp <- 10 #The number of irrelevant variables
#yout is true for one view by default
#Total data dimension is D = p + Irrp

GenerateMultiSampleData <- function(n, views, K, w, p, Irrp){
  #Create a list of the sample data and merge the data into one, maintain different cluster labellings
  dataMatrices <- list()
  dataMatrices[[1]] <- GenerateSampleData(n, K[1], w[[1]], p = p[1], Irrp = Irrp, yout = TRUE)
  #This has the outcome column and also the irrelevant variables (null view)
  for (v in 2:views){
    dataMatrices[[v]] <- GenerateSampleData(n, K[v], w[[v]], p = p[v], Irrp = 0, yout = FALSE)
    #irrelevant views
  }
  
  #Merge data into one data matrix
  #First entry in list: y is a separate entry (third entry) in the list
  dataMatrix <- list()
  for (v in 1:views){
    dataMatrix[[v]] <- dataMatrices[[v]][[1]]
  }
  dataMatrix <- do.call(cbind, dataMatrix)
  
  final_output <- list()
  final_output$dataMatrix <- dataMatrix
  final_output$outcome <- dataMatrices[[1]][[3]]
  final_output$rel_labels <- dataMatrices[[1]][[2]]
  for (v in 2:views){
    final_output[[v + 2]] <- dataMatrices[[v]][[2]]
  }
  return(final_output)
}