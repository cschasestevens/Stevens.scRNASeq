#' Cell Type Prediction
#'
#' Predicts cell type identities of clusters in a Seurat object given a previously annotated data set. The reference data
#' must have a cell type column named 'CellType' for predictions to be successful.
#'
#' @param so An object of class Seurat.
#' @param d.ref An annotated Seurat object containing cell type identities for cell type prediction.
#' @param f.size Memory allocated for objects created during cell type prediction (in MB), provided as a numeric value. For larger datasets, set to 5000 or higher.
#' @param cl.var Character string containing the name of the cluster variable for cell type predictions.
#' @param md.list A vector of character strings indicating metadata columns to include with the prediction score output table.
#' @param parl Logical indicating whether processing should be run in parallel (Linux and WSL2 only). Set to FALSE if running sequentially.
#' @param core.perc Percentage of available cores to use if running in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @return A list containing a table of marker genes per cluster and an annotated Seurat object with putative cell identities.
#' @examples
#'
#' sc.predict.clusters(d.seurat,d.ref,5000,"seurat_clusters",c("col1","col2","col3"),TRUE,0.5)
#'
#' @export
sc.predict.clusters <- function(
    so,
    d.ref,
    f.size,
    cl.var,
    md.list,
    parl,
    core.perc
    ) {
  dr <- d.ref
  d <- so
  # Change assay and active identity
  SeuratObject::DefaultAssay(dr) <- "RNA"
  SeuratObject::DefaultAssay(d) <- "RNA"
  d <- SeuratObject::SetIdent(d,value = d@meta.data[[cl.var]])






  if(Sys.info()[["sysname"]] != "Windows" & parl == TRUE) {
    ## Change to multisession (runs in parallel)
    future::plan("multisession",workers = parallel::detectCores()*core.perc)
    options(future.globals.maxSize = f.size * 1024^2)
    ## Find marker genes for all clusters (for manual annotation of clusters)
    if(file.exists("analysis/table.marker.genes.txt")) {
      cl.mark <- read.table(
        "analysis/table.marker.genes.txt",
        sep = "\t",
        header = T
        )}
    if(!file.exists("analysis/table.marker.genes.txt")) {
      cl.mark <- Seurat::FindAllMarkers(
        d,
        verbose = T)
      write.table(
        cl.mark,
        "analysis/table.marker.genes.txt",
        col.names = T,
        row.names = F,
        sep = "\t"
        )
      }
    ## Cell Type Prediction
    if(
      file.exists("analysis/data.with.predictions.rds") &
      file.exists("analysis/table.predicted.types.txt")) {
      list.d <- readRDS("analysis/data.with.predictions.rds")
      list.d[["Prediction Scores"]] <- read.table(
        "analysis/table.predicted.types.txt",
        sep = "\t",
        header = T
      )
    }
    if(
      !file.exists("analysis/data.with.predictions.rds") &
      !file.exists("analysis/table.predicted.types.txt")) {
    ### Reference features
    data.ref <- Seurat::FindVariableFeatures(
      dr,
      selection.method = "vst",
      nfeatures = 4000,
      assay = "RNA"
      )
    ### Query features
    data.qry <- Seurat::FindVariableFeatures(
      d,
      selection.method = "vst",
      nfeatures = 4000,
      assay = "RNA"
      )
    ## Transfer anchors
    data.transfer.anchors <- Seurat::FindTransferAnchors(
      reference = data.ref,
      query = data.qry,
      features = Seurat::VariableFeatures(object = data.ref),
      reference.assay = "RNA",
      query.assay = "RNA",
      reduction = "pcaproject")
    ## Transfer cell type predictions to query set
    data.predicted <- Seurat::TransferData(
      anchorset = data.transfer.anchors,
      refdata = data.ref@meta.data[["CellType"]],
      weight.reduction = "pcaproject"
      )
    # change to sequential
    future::plan("sequential")
    options(future.globals.maxSize = 500 * 1024^2)

    ## add metadata and return the query dataset with cell type predictions
    data.qry <- Seurat::AddMetaData(
      data.qry,
      metadata = data.predicted
    )

    # Combine marker genes and Seurat object as list
    list.d <- list(
      "Predicted Clusters" = data.qry,
      "Marker Genes" = cl.mark
    )
    ## Save predictions as table and Seurat object with predictions
    saveRDS(list.d,"analysis/data.with.predictions.rds")
    list.d[["Prediction Scores"]] <- dplyr::select(
      list.d$`Predicted Clusters`@meta.data,
      c(md.list,cl.var,names(
        list.d$`Predicted Clusters`@meta.data[grepl(
          "predicted|prediction",
          names(list.d$`Predicted Clusters`@meta.data))])))
    write.table(
      list.d$`Prediction Scores`,
      "analysis/table.predicted.types.txt",
      col.names = T,
      row.names = F,
      sep = "\t"
      )
    }
  }

  if(Sys.info()[["sysname"]] == "Windows") {
    ## Change to multisession (runs in parallel)
    ## Find marker genes for all clusters (for manual annotation of clusters)
    if(file.exists("analysis/table.marker.genes.txt")) {
      cl.mark <- read.table(
        "analysis/table.marker.genes.txt",
        sep = "\t",
        header = T
        )}
    if(!file.exists("analysis/table.marker.genes.txt")) {
      cl.mark <- Seurat::FindAllMarkers(
        d,
        verbose = T)
      write.table(
        cl.mark,
        "analysis/table.marker.genes.txt",
        col.names = T,
        row.names = F,
        sep = "\t"
        )
      }
    ## Cell Type Prediction
    if(
      file.exists("analysis/data.with.predictions.rds") &
      file.exists("analysis/table.predicted.types.txt")) {
      list.d <- readRDS("analysis/data.with.predictions.rds")
      list.d[["Prediction Scores"]] <- read.table(
        "analysis/table.predicted.types.txt",
        sep = "\t",
        header = T
      )
    }
    if(
      !file.exists("analysis/data.with.predictions.rds") &
      !file.exists("analysis/table.predicted.types.txt")) {
    ### Reference features
    data.ref <- Seurat::FindVariableFeatures(
      dr,
      selection.method = "vst",
      nfeatures = 4000,
      assay = "RNA"
      )
    ### Query features
    data.qry <- Seurat::FindVariableFeatures(
      d,
      selection.method = "vst",
      nfeatures = 4000,
      assay = "RNA"
      )
    ## Transfer anchors
    data.transfer.anchors <- Seurat::FindTransferAnchors(
      reference = data.ref,
      query = data.qry,
      features = Seurat::VariableFeatures(object = data.ref),
      reference.assay = "RNA",
      query.assay = "RNA",
      reduction = "pcaproject")
    ## Transfer cell type predictions to query set
    data.predicted <- Seurat::TransferData(
      anchorset = data.transfer.anchors,
      refdata = data.ref@meta.data[["CellType"]],
      weight.reduction = "pcaproject"
      )

    ## add metadata and return the query dataset with cell type predictions
    data.qry <- Seurat::AddMetaData(
      data.qry,
      metadata = data.predicted
      )

    # Combine marker genes and Seurat object as list
    list.d <- list(
      "Predicted Clusters" = data.qry,
      "Marker Genes" = cl.mark
      )
    ## Save predictions as table and Seurat object with predictions
    saveRDS(list.d,"analysis/data.with.predictions.rds")
    list.d[["Prediction Scores"]] <- dplyr::select(
      list.d$`Predicted Clusters`@meta.data,
      c(md.list,cl.var,names(
        list.d$`Predicted Clusters`@meta.data[grepl(
          "predicted|prediction",
          names(list.d$`Predicted Clusters`@meta.data))])))
    write.table(
      list.d$`Prediction Scores`,
      "analysis/table.predicted.types.txt",
      col.names = T,
      row.names = F,
      sep = "\t"
      )

    }
  }

  if(Sys.info()[["sysname"]] != "Windows" & parl == FALSE) {
    ## Change to multisession (runs in parallel)
    ## Find marker genes for all clusters (for manual annotation of clusters)
    if(file.exists("analysis/table.marker.genes.txt")) {
      cl.mark <- read.table(
        "analysis/table.marker.genes.txt",
        sep = "\t",
        header = T
      )}
    if(!file.exists("analysis/table.marker.genes.txt")) {
      cl.mark <- Seurat::FindAllMarkers(
        d,
        verbose = T)
      write.table(
        cl.mark,
        "analysis/table.marker.genes.txt",
        col.names = T,
        row.names = F,
        sep = "\t"
      )
    }

    ## Cell Type Prediction
    if(
      file.exists("analysis/data.with.predictions.rds") &
      file.exists("analysis/table.predicted.types.txt")) {
      list.d <- readRDS("analysis/data.with.predictions.rds")
      list.d[["Prediction Scores"]] <- read.table(
        "analysis/table.predicted.types.txt",
        sep = "\t",
        header = T
        )
      }
    if(
      !file.exists("analysis/data.with.predictions.rds") &
      !file.exists("analysis/table.predicted.types.txt")) {
    ### Reference features
    data.ref <- Seurat::FindVariableFeatures(
      dr,
      selection.method = "vst",
      nfeatures = 4000,
      assay = "RNA"
      )
    ### Query features
    data.qry <- Seurat::FindVariableFeatures(
      d,
      selection.method = "vst",
      nfeatures = 4000,
      assay = "RNA"
      )
    ## Transfer anchors
    data.transfer.anchors <- Seurat::FindTransferAnchors(
      reference = data.ref,
      query = data.qry,
      features = Seurat::VariableFeatures(object = data.ref),
      reference.assay = "RNA",
      query.assay = "RNA",
      reduction = "pcaproject")
    ## Transfer cell type predictions to query set
    data.predicted <- Seurat::TransferData(
      anchorset = data.transfer.anchors,
      refdata = data.ref@meta.data[["CellType"]],
      weight.reduction = "pcaproject"
      )

  ## add metadata and return the query dataset with cell type predictions
  data.qry <- Seurat::AddMetaData(
    data.qry,
    metadata = data.predicted
    )

  # Combine marker genes and Seurat object as list
  list.d <- list(
    "Predicted Clusters" = data.qry,
    "Marker Genes" = cl.mark
    )
  ## Save predictions as table and Seurat object with predictions
  saveRDS(list.d,"analysis/data.with.predictions.rds")
  list.d[["Prediction Scores"]] <- dplyr::select(
    list.d$`Predicted Clusters`@meta.data,
    c(md.list,cl.var,names(
      list.d$`Predicted Clusters`@meta.data[grepl(
        "predicted|prediction",
        names(list.d$`Predicted Clusters`@meta.data))])))
  write.table(
    list.d$`Prediction Scores`,
    "analysis/table.predicted.types.txt",
    col.names = T,
    row.names = F,
    sep = "\t"
    )
  }
  }

  print(
    table(
      list.d$`Predicted Clusters`$prediction.score.max >
        0.5
      )
    )
  ### Cell Type Prediction Score Distribution
  fun.dist.score <- function(x)
    {
    p.score.dist <- ggplot2::ggplot(
      x,
      ggplot2::aes(
        x = .data[["prediction.score.max"]]
        )
      ) +
      ggplot2::geom_density(
        color = "black",
        fill = col.univ()[[2]]
        ) +
      ggplot2::labs(
        x = "Prediction Score",
        y = "Density",
        title = "Prediction Score Distribution"
      ) +
      sc.theme1()
    }
  ggplot2::ggsave(
    "analysis/plot.predicted.scores.dist.png",
    fun.dist.score(
      list.d$`Prediction Scores`
    ),
    width = 8,
    height = 8,
    dpi = 700
    )

  ### Predicted Cell Type Proportions for Each Cluster
  # Counting function
  fun.predict.prop.alt <- function(
    x,
    c,
    md
    ) {
    data.pred.prop2 <- setNames(
      dplyr::count(
        x,
        .data[[c]]
        ),
      c(c,"Total Cells")
      )
    data.pred.prop <- dplyr::count(
      x,
      x[,c(md,c)]
      )
    data.pred.prop <- dplyr::left_join(
      data.pred.prop,
      data.pred.prop2,
      by = c
      )
    data.pred.prop[["Proportion"]] <- round(
      data.pred.prop$n/
        data.pred.prop$`Total Cells`,
      digits = 3
      )
    return(data.pred.prop)
    }

  ## Cell Type Proportion Summary and Consensus Type
  list.d[["Predicted Proportions"]] <- dplyr::filter(
    fun.predict.prop.alt(
      list.d[["Predicted Clusters"]]@meta.data,
      "predicted.id",
      c(md.list,cl.var)
      ),
      .data[["Proportion"]] >
        0.001
      )
  write.table(
    list.d[["Predicted Proportions"]],
    "analysis/table.predicted.prop.txt",
    col.names = T,
    row.names = F,
    sep = "\t"
    )
  list.d[["pred.prop.summary"]] <- setNames(
    aggregate(
      list.d[["Predicted Proportions"]][["Proportion"]],
      list(
        list.d[["Predicted Proportions"]][[cl.var]],
        list.d[["Predicted Proportions"]][["predicted.id"]]
        ),
      FUN = sum
      ),
    c(cl.var,"predicted.id","Proportion")
    )
  write.table(
    list.d[["pred.prop.summary"]],
    "analysis/table.predicted.prop.summary.txt",
    col.names = T,
    row.names = F,
    sep = "\t"
    )

  ## Return cell type predictions and assign highest ranked type to each cell
  list.d[["cluster.proportions"]] <- dplyr::select(
    unique(
      dplyr::left_join(
        setNames(
          aggregate(
            list.d$pred.prop.summary[["Proportion"]],
            list(
              list.d$pred.prop.summary[[cl.var]]),
            FUN = max
            ),
          c(cl.var,
            "Proportion"
            )
          ),
        list.d$pred.prop.summary[,c("predicted.id","Proportion")],
        by = c("Proportion")
        )
      ),
    c(cl.var,"predicted.id","Proportion")
    )
  list.d[["cluster.proportions"]] <- list.d[["cluster.proportions"]][!duplicated(list.d[["cluster.proportions"]][[cl.var]]),]
  list.d[["cluster.assignments"]] <- dplyr::mutate(
    list.d$cluster.proportions,
    "predicted.id" = paste(
      seq(1:nrow(list.d$cluster.proportions)),
      gsub("FOXN4","",
           gsub("\\.","",
                gsub("^*..","",
                     list.d$cluster.proportions$predicted.id
                     )
                )
           ),
      sep = "."),
    "CellGroup" = as.factor(
      gsub("FOXN4","",
           gsub("\\.","",
                gsub("^*..","",
                     list.d$cluster.proportions$predicted.id
                     )
                )
           )
      )
    )
  ### Change grouping columns to factors
  list.d[["cluster.assignments"]][["predicted.id"]] <- factor(
    list.d[["cluster.assignments"]][["predicted.id"]],
    levels = c(gtools::mixedsort(list.d[["cluster.assignments"]][["predicted.id"]]))
    )
  list.d[["cluster.assignments"]][["CellGroup"]] <- factor(
    list.d[["cluster.assignments"]][["CellGroup"]],levels = c(list.ct))
  if("CellGroup" %in% names(list.d[["Predicted Clusters"]]@meta.data) == TRUE){
  list.d[["cluster.assignments"]] <- dplyr::select(
    dplyr::left_join(
      list.d$`Predicted Clusters`@meta.data,
      list.d$cluster.assignments,
      by = cl.var
      ),
    c("predicted.id.y","CellGroup.y"))

    list.d$`Predicted Clusters` <- Seurat::AddMetaData(
      list.d$`Predicted Clusters`,
      list.d$cluster.assignments[["CellGroup.y"]],
      col.name = "CellGroup"
    )

  }
  if("CellGroup" %in% names(list.d[["Predicted Clusters"]]@meta.data) == FALSE){
    list.d[["cluster.assignments"]] <- dplyr::select(
      dplyr::left_join(
        list.d$`Predicted Clusters`@meta.data,
        list.d$cluster.assignments,
        by = cl.var
      ),
      c("predicted.id.y","CellGroup"))

      list.d$`Predicted Clusters` <- Seurat::AddMetaData(
        list.d$`Predicted Clusters`,
        list.d$cluster.assignments[["CellGroup"]],
        col.name = "CellGroup"
      )
    }
  ### Add Cell Type and Cell Group columns to seurat object
  list.d$`Predicted Clusters` <- Seurat::AddMetaData(
    list.d$`Predicted Clusters`,
    list.d$cluster.assignments[["predicted.id.y"]],
    col.name = "CellType"
    )

  return(list.d)
  }




