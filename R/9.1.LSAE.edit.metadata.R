#### Assign Metadata Columns for Merged Original and NKX2.1 KO Dataset ####

# Read RDS (if not analyzing sequentially)

# d.analysis <- readRDS("Processed/1.all.lsae.rds")
  
  # Extract metadata for modifying columns
  
  d.analysis.meta <- list.analysis@meta.data
  
  ## Return columns that contain NAs
  
  d.na <- colnames(
    d.analysis.meta
  )[apply(
    d.analysis.meta,
    2,
    anyNA
  )]
  
  setNames(
    lapply(
      d.na,
      function(x) 
        unique(
          d.analysis.meta[,x]
        )
    ),
    d.na
  )
  
  
  ## Edit columns
  
  d.analysis.meta[["Time"]] <- factor(
    ifelse(
      is.na(d.analysis.meta[["Time"]]),
      "D28",
      d.analysis.meta[["Time"]]
    ),
    levels = c(
      "D28")
  )
  
  
  d.analysis.meta[["Knockout"]] <- factor(
    ifelse(
      is.na(
        d.analysis.meta[["Knockout"]]
      ) |
        d.analysis.meta[["Knockout"]] == "NG",
      "Control",
      "KO"
    ),
    levels = c(
      "Control",
      "KO")
  )
  
  
  d.analysis.meta[["Batch"]] <- factor(
    ifelse(
      is.na(d.analysis.meta[["Batch"]]),
      "nkx",
      d.analysis.meta[["Batch"]]
    ),
    levels = c(
      "b1",
      "b2",
      "nkx")
  )
  
  
  ## Check revised columns
  
  d.na <- colnames(
    d.analysis.meta
  )[apply(
    d.analysis.meta,
    2,
    anyNA
  )]
  
  
  # Assign CellTypes based on original LSAE dataset
  
  ## View unique cell type clusters
  
  setNames(
    lapply(
      names(
        d.analysis.meta
      ),
      function(x) 
        unique(
          d.analysis.meta[,x]
        )
    ),
    names(
      d.analysis.meta
    )
  )
  
  
  ## Assign preliminary cell types based on Seurat clusters from original LSAE dataset
  
  ### Clusters
  
  d.clus <- setNames(
    d.analysis.meta[d.analysis.meta$Batch == "b1" |
                      d.analysis.meta$Batch == "b2",
                    c(
                      "seurat_clusters",
                      "CellType"
                    )],
    c(
      "seurat_clusters",
      "Type"
    )
  )
  
  d.clus <- unique(d.clus)
  
  
  # Merge updated meta data columns with main analysis object
  
  list.analysis <- AddMetaData(
    list.analysis,
    d.analysis.meta[["Batch"]],
    col.name = "Batch"
  )
  
  list.analysis <- AddMetaData(
    list.analysis,
    d.analysis.meta[["Knockout"]],
    col.name = "Knockout"
  )
  
  list.analysis <- AddMetaData(
    list.analysis,
    d.analysis.meta[["Time"]],
    col.name = "Time"
  )

  list.analysis <- AddMetaData(
    list.analysis,
    as.factor(d.analysis.meta[["Airway"]]),
    col.name = "Airway"
  )














