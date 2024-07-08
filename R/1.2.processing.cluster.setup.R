#### Define Cluster for Parallel Processing ####

## CLUSTER SETUP: creates a cluster for parallelization of script execution using 50% of available cores or length of listed data files


## Specify number of nodes

clus1 <- makeCluster(num.core)

## Pass libraries and functions to cluster

clusterEvalQ(clus1,{
  
  ### Libraries
  
  library(Seurat)
  library(ggplot2)
  library(SoupX)
  library(SingleCellExperiment)
  library(scDblFinder)

  ### Functions
  
  
  
  })


## Export global environment variables to cluster

clusterExport(clus1,
              varlist = c("list.params",
                          "thm.mult","thm.univ",
                          "thm.gen","col1a",
                          "col1b","col2a",
                          "col3a","num.core"))








