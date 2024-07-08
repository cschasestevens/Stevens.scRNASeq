#### Cell type predictions and marker genes for cluster annotation ####

fun.predict.clusters <- function(
    d.ref,
    d.qry,
    f.size
    ) {
  
  #---- Input Parameters ----
  
  ## cluster number and cores to use
  
  num.clus <- length(
    as.numeric(
      unique(
        levels(
          d.qry@meta.data[["seurat.clusters"]]
        )
      )
    )
  )
  
  num.core <- detectCores()
  
  ## Change default assay to `RNA`
  
  DefaultAssay(d.ref) <- "RNA"
  DefaultAssay(d.qry) <- "RNA"
  
  ## Change to multisession (runs in parallel)
  
  plan(
    "multisession",
    workers = ifelse(
      num.clus < num.core,
      ceiling(
        num.clus*0.8
        ),
      ceiling(
        num.core*0.8
        )
      )
    )
  
  options(
    future.globals.maxSize = f.size * 1024^2
    )
  
  
  ## Find marker genes for all clusters (for manual annotation of clusters)
  
  if(
    file.exists(
      paste(
        list.p$a.path,
        "marker.genes.txt",
        sep = ""
      ) 
    )
  ) {
    
    cl.mark <- read.table(
      paste(
        list.p$a.path,
        "marker.genes.txt",
        sep = ""
      ),
      sep = "\t",
      header = T
    )
    
  }
  
  if(
    !file.exists(
      paste(
        list.p$a.path,
        "marker.genes.txt",
        sep = ""
      ) 
    )
  ) {
    
    cl.mark <- FindAllMarkers(
      d.qry,
      verbose = T
    )
    
    write.table(
      cl.mark,
      paste(
        list.p$a.path,
        "marker.genes.txt",
        sep = ""
      ),
      col.names = T,
      row.names = F,
      sep = "\t"
    )
    
  }
  
  
  
  
  
  #---- Cell Type Prediction ----
  
  ## Reference features
  
  data.ref <- FindVariableFeatures(
    d.ref,
    selection.method = "vst",
    nfeatures = 4000,
    assay = "RNA"
    )
  
  ## Query features
  
  data.qry <- FindVariableFeatures(
    d.qry,
    selection.method = "vst",
    nfeatures = 4000,
    assay = "RNA"
  )
  
  
  ## Transfer anchors
  
  data.transfer.anchors <- FindTransferAnchors(
    reference = data.ref,
    query = data.qry,
    features = VariableFeatures(object = data.ref),
    reference.assay = "RNA",
    query.assay = "RNA",
    reduction = "pcaproject")
  
  
  ## Transfer cell type predictions to query set
  
  data.predicted <- TransferData(
    anchorset = data.transfer.anchors,
    refdata = data.ref@meta.data[["CellType"]],
    weight.reduction = "pcaproject"
    )
  
  
  ## change back to sequential
  
  plan(
    "sequential"
    )
  
  options(
    future.globals.maxSize = 500 * 1024^2
    )
  
  
  ## add metadata and return the query dataset with cell type predictions
  
  data.qry <- AddMetaData(
    data.qry,
    metadata = data.predicted)
  
  return(
    list(
      "Predicted Clusters" = data.qry,
      "Marker Genes" = cl.mark
      )
    )
  
}
























