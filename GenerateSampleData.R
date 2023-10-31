#Generate example binary clustered data

library(pheatmap)
library(RColorBrewer)

set.seed(463)

#Example settings for parameters in function
n    <- 1000 # The overall population size
K    <- 4    # The number of clusters/subpopulations
w    <- c(0.1, 0.2, 0.3, 0.4) # The mixture weights (proportion of population in each cluster)
w    <- w/sum(w)   # This ensures that the elements of w sum to 1 
p    <- 20 # The number of variables 
Irrp <- 5 #The number of irrelevant variables
yout <- TRUE #Indicate whether a binary outcome is required
#Total data dimension is D = p + Irrp

#Returns a list of the data matrix as well as the 'true' cluster labels for each observation
GenerateSampleData <- function(n, K, w, p, Irrp, yout = FALSE){
  # Variable called "clusterLabels" to store the cluster labels, 

  clusterLabels <- sample(1:K, n, replace = T, prob = w)
  
  clusterParameters <- matrix(nrow = K, ncol = p)
  for(i in 1:K)
  {
    for(j in 1:p)
    {
      clusterParameters[i,j] <- rbeta(1, 1, 5) # Can tweak the Beta parameters to control, e.g. sparsity
      #ie. generate one parameter for each variable in each cluster via Beta dist
    }
  }
  
  if(Irrp != 0){
    unclusteredParameters <- rbeta(Irrp, 1, 5)
  }
  
  # Use cluster labels and cluster parameters to generate sample data

  dataMatrix    <- matrix(nrow = n, ncol = p + Irrp)
  for(i in 1:n)
  {
    currentClusterLabel <- clusterLabels[i]  # Get the cluster label for the i-th person
    for(j in 1:p)
    {
      # Simulate data according to the parameters for the current cluster
      dataMatrix[i,j] <- rbinom(1,1, clusterParameters[currentClusterLabel,j]) 
      #Generated via Bernoulli trial and uses beta probability
    }
    if (Irrp != 0){
      for (j in 1:Irrp){
        #Simulate data according to random probability, no clustering structure
        dataMatrix[i,p + j] <- rbinom(1, 1, unclusteredParameters[j])
      }
    }
  }
  
  if (yout == TRUE){
    y <- vector(mode = "numeric", length = n)
    yParams <- vector(mode = "numeric", length = K) 
    
    prob <- runif(1)
    for (k in 1:K){
      yParams[k] <- rbeta(1, prob, prob) #generate random Beta probability
    }
    
    for (i in 1:n){
      currentClusterLabel <- clusterLabels[i]
      y[i] <- rbinom(1, 1, yParams[currentClusterLabel])
    }
    
  }
  
  #Create an "annotation row" to show the cluster label for each person
  annotationRow <- data.frame(
    Cluster = factor(clusterLabels)
  )
  #Create data frame for outcome for each person if needed
  if (yout == TRUE){
    y <- data.frame(
      outcome = factor(y)
    )
  }
  
  #We require the data, outcome and annotation row to have the same rownames:
  rownames(annotationRow) <- rownames(dataMatrix) <- paste0("Person", seq(1,n))
  if(yout == TRUE){
    rownames(y) <- rownames(dataMatrix)
  }
  
  if(yout == TRUE){
    return(list(dataMatrix, annotationRow, y))
  } else{
    return(list(dataMatrix, annotationRow))
  }
  }

######Visualising our simulated data

#If binary outcome is present:
annotationRow <- cbind(outcome = dataMatrix[[3]], Clusters = dataMatrix[[2]])
annotationColor <- list(
  outcome  = c("0" = "white", "1" = "black"),
  Clusters = c("1" = "#1B9E77", "2" = "#D95F02", "3" = "#7570B3", "4" = "#E7298A")
)


pheatmap(dataMatrix[[1]], show_rownames = F, annotation_row = annotationRow, 
         color = colorRampPalette(colors = c("white", "black"))(2), 
         annotation_colors = annotationColor,
         cluster_rows = F, cluster_cols = F)
pheatmap(dataMatrix[[1]][sort(as.numeric(dataMatrix[[2]]$Cluster), index.return = T)$ix,], 
         show_rownames = F, annotation_row = annotationRow, 
         annotation_colors = annotationColor,
         color = colorRampPalette(colors = c("white", "black"))(2),
         cluster_rows = F)

#If binary outcome is not present
pheatmap(dataMatrix[[1]], show_rownames = F, annotation_row = dataMatrix[[2]], 
         color = colorRampPalette(colors = c("white", "black"))(2), 
         cluster_rows = F, cluster_cols = F)
pheatmap(dataMatrix[[1]][sort(as.numeric(dataMatrix[[2]]$Cluster), index.return = T)$ix,], 
         show_rownames = F, annotation_row = dataMatrix[[2]], 
         color = colorRampPalette(colors = c("white", "black"))(2),
         cluster_rows = F)

#Visualising with outcome

#Run this after fitting model 

currentAnnotationRow <- data.frame(
  Cluster = factor(vidataMatrix$model$labels),
  trueClusters = factor(dataMatrix[[2]]$Cluster)
)
rownames(currentAnnotationRow) <- rownames(dataMatrix[[1]]) 
pheatmap(dataMatrix[[1]][sort(vidataMatrix$model$labels, index.return = T)$ix,], 
         cluster_rows = F, show_rownames = F, show_colnames = F, 
         color = colorRampPalette(colors = c("white", "black"))(2),
         annotation_row = currentAnnotationRow)

