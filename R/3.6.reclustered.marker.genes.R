#### Cell type predictions and marker genes for cluster annotation ####

fun.marker.gene.recluster <- function(
    d.qry,
    title1
    ) {
  
  #---- Input Parameters ----

  ## Change default assay to `RNA`
  
  DefaultAssay(d.qry) <- "RNA"
  
  
  ## Find marker genes for all clusters
  
  if(
    file.exists(
      paste(
        title1,
        sep = ""
      ) 
    )
  ) {
    
    cl.mark <- read.table(
      paste(
        title1,
        sep = ""
      ),
      sep = "\t",
      header = T
    )
    
  }
  
  if(
    !file.exists(
      paste(
        title1,
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
        title1,
        sep = ""
      ),
      col.names = T,
      row.names = F,
      sep = "\t"
    )
    
  }
  
  return(
    cl.mark
  )
  
}
























