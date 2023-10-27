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
#Total data dimension is D = p + Irrp

#Returns a list of the data matrix as well as the 'true' cluster labels for each observation
GenerateSampleData <- function(n, K, w, p, Irrp){
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
  
  unclusteredParameters <- rbeta(Irrp, 1, 5)
  
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
    for (j in 1:Irrp){
      #Simulate data according to random probability, no clustering structure
      dataMatrix[i,p + j] <- rbinom(1, 1, unclusteredParameters[j])
    }
  }
  
  #Create an "annotation row" to show the cluster label for each person
  annotationRow <- data.frame(
    Cluster = factor(clusterLabels)
  )
  #We require the data and annotation row to have the same rownames:
  rownames(annotationRow) <- rownames(dataMatrix) <- paste0("Person", seq(1,n))
  return(list(dataMatrix, annotationRow))
  }

#Visualising our simulated data
pheatmap(dataMatrix[[1]], show_rownames = F, annotation_row = annotationRow, 
         color = colorRampPalette(colors = c("white", "black"))(2), 
         cluster_rows = F, cluster_cols = F)
pheatmap(dataMatrix[[1]][sort(as.numeric(dataMatrix[[2]]$Cluster), index.return = T)$ix,], 
         show_rownames = F, annotation_row = dataMatrix[[2]], 
         color = colorRampPalette(colors = c("white", "black"))(2),
         cluster_rows = F)

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

