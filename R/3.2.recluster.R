#' Reclustering Analysis
#'
#' Performs reclustering and basic visualization of a subsetted Seurat Object. Subsets are determined based on CellGroup.
#'
#' @param so An object of class Seurat. Must contain the columns 'CellType' and 'CellGroup' in the metadata slot.
#' @param ct Character string vector containing the name(s) of up to 2 cell type groups present in the CellGroup column.
#' The provided Seurat object is then subsetted to include only clusters from the provided group names.
#' @param g.list Path to an existing text file containing genes for plotting expression.
#' @param md.list A vector of character strings indicating metadata columns for overlaying on a loadings plot.
#' @param h.w Numeric value for marker gene heatmap width (passed to ComplexHeatmap).
#' @param h.h Numeric value for marker gene heatmap height (passed to ComplexHeatmap).
#' @param fs.c Numeric value for marker gene heatmap column fontsize (passed to ComplexHeatmap).
#' @param fs.r Numeric value for marker gene heatmap row fontsize (passed to ComplexHeatmap).
#' @param parl Logical indicating whether processing should be run in parallel (Linux and WSL2 only). Set to FALSE if running sequentially.
#' @param core.perc Percentage of available cores to use if running in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @param slot1 A character string indicating which assay should be used for dimension reduction.
#' @param slot2 A character string indicating the reduction name for plotting gene expression UMAPs.
#' @param run.qc A logical indicating whether gene expression QC plots should be generated for the reclustered data.
#' @return A reclustered Seurat Object with summary UMAP plots and a marker gene list.
#' @examples
#'
#' # p.umap <- sc.top10.marker.heatmap(d.annotated,"seurat.clusters",18,24,6,8)
#'
#' @export
sc.recluster.data <- function(
    so,
    ct,
    g.list,
    md.list,
    h.w,
    h.h,
    fs.c,
    fs.r,
    parl,
    core.perc,
    slot1,
    slot2,
    run.qc
    ) {
  # Load gene list and Seurat object
  list.genes <- read.table(
    g.list,
    sep = "\t",
    header = T)
  d <- so

  # Subset Seurat
  if(length(ct) == 2) {
    d <- BiocGenerics::subset(
      d,
      subset = CellGroup == ct[[1]]|CellGroup == ct[[2]]
      )
    }
  if(length(ct) == 1) {
    d <- BiocGenerics::subset(
      d,
      subset = CellGroup == ct
      )
    }

  # Run PCA/UMAP
  d <- sc_pca(d,"RNA")
  d.pca <- sc_pca_plot(
    d,
    c(md.list,"seurat_clusters")
    )
  ggplot2::ggsave(
    "analysis/recluster/plot.pca.panel.png",
    d.pca,
    width = 36,
    height = 12,
    dpi = 700)

  # Determine optimal number of dimensions to use, then run UMAP
  d <- SeuratObject::SetIdent(d, value = d@meta.data$CellType)
  d <- Seurat::RunUMAP(
    object = d,
    dims = 1:20,
    n.components = 3
    )
  SeuratObject::DefaultAssay(d) <- slot1
  d <- Seurat::FindNeighbors(
    d,
    reduction = "umap",
    dims = 1:3,
    verbose = T
    )
  d <- Seurat::FindClusters(
    d,
    cluster.name = "recluster",
    resolution = 0.3
    )
  d.umap1 <- sc.umap.panel(
    d,
    c("CellType","recluster",md.list),
    slot2
    )
  ggplot2::ggsave(
    "analysis/recluster/plot.umap.panel.png",
    d.umap1,
    height = 16,
    width = 32,
    dpi = 400
    )

  # Recluster QC
  if(run.qc == T){
    ggplot2::ggsave(
      "analysis/recluster/plot.recluster.qc.png",
      sc_integration_qc(d,"recluster"),
      width = 24,
      height = 8,
      dpi = 700)
    }
  
  # Recluster UMAP
  ### Cell Type
  ggplot2::ggsave(
    "analysis/recluster/plot.umap.CellType.png",
    sc_umap_standard(
      # Seurat object
      d,
      # metadata column
      "CellType",
      slot2
      ),
    height = 8,
    width = 8,
    dpi = 700
    )
  ### Cell Group
  ggplot2::ggsave(
    "analysis/recluster/plot.umap.reclustering.png",
    sc_umap_standard(
      # Seurat object
      d,
      # metadata column
      "recluster",
      slot2
      ),
    height = 8,
    width = 8,
    dpi = 700
    )

  ## Create Heatmap
  p.heatmap.top10 <- sc.top10.marker.heatmap.rc(
    # Seurat object
    d,
    # Cluster column
    "recluster",
    # Heatmap width
    h.w,
    # Heatmap height
    h.h,
    # Column fontsize
    fs.c,
    # Row fontsize
    fs.r
    )

  ## Save
  png(
    "analysis/recluster/plot.heatmap.top10.markers.png",
    width = 40,
    height = 14,
    units = "cm",
    res = 1000
  )
  print(p.heatmap.top10)
  dev.off()

  # Predicted Cell Type Proportions for Each Cluster
  ## Counting function
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
  d.prop <- dplyr::filter(
    fun.predict.prop.alt(
      d@meta.data,
      "recluster",
      c(md.list)
      ),
    .data[["Proportion"]] >
      0.001
    )
  write.table(
    d.prop,
    "analysis/recluster/table.recluster.prop.txt",
    col.names = T,
    row.names = F,
    sep = "\t"
    )
  d.prop.sum <- setNames(
    aggregate(
      d.prop[["Proportion"]],
      list(
        d.prop[["recluster"]]
        ),
      FUN = sum
      ),
    c("recluster","Proportion")
    )
  write.table(
    d.prop.sum,
    "analysis/recluster/table.recluster.prop.summary.txt",
    col.names = T,
    row.names = F,
    sep = "\t"
    )

  # Reclustered gene expression plots
  sc.umap.panel.gene.list <- function(list.g,so,md.var,col.scheme,col.names,leg.x,leg.y,parl,core.perc,slot2) {
    lg <- list.g
    d <- so
    lg <- unique(lg[lg %in% SeuratObject::Features(d)])
    lg.abs <- subset(lg, !(lg %in% SeuratObject::Features(d)))
    # Create plots
    if(Sys.info()[["sysname"]] != "Windows" &
       parl == TRUE){
      parallel::mclapply(
        mc.cores = ceiling(
          parallel::detectCores()*
            core.perc
        ),
        lg,
        function(x) {
          pg <- sc_umap_panel_gene(
            d,
            md.var,
            x,
            col.scheme,
            col.names,
            leg.x,
            leg.y,
            slot2
          )
          # Save each plot
          ggplot2::ggsave(
            paste("analysis/recluster/gene/plot.umap.exp.",x,".png",sep = ""),
            pg,
            height = 12,
            width = 36,
            dpi = 700
          )
        }
      )

    }
    if(Sys.info()[["sysname"]] == "Windows"){
      lapply(
        lg,
        function(x) {
          pg <- sc_umap_panel_gene(
            d,
            md.var,
            x,
            col.scheme,
            col.names,
            leg.x,
            leg.y,
            slot2
          )
          # Save each plot
          ggplot2::ggsave(
            paste("analysis/recluster/gene/plot.umap.exp.",x,".png",sep = ""),
            pg,
            height = 12,
            width = 36,
            dpi = 700
          )
        }
      )

    }
    if(Sys.info()[["sysname"]] != "Windows" & parl == FALSE){
      lapply(
        lg,
        function(x) {
          pg <- sc_umap_panel_gene(
            d,
            md.var,
            x,
            col.scheme,
            col.names,
            leg.x,
            leg.y,
            slot2
          )
          # Save each plot
          ggplot2::ggsave(
            paste("analysis/recluster/gene/plot.umap.exp.",x,".png",sep = ""),
            pg,
            height = 12,
            width = 36,
            dpi = 700
          )
        }
      )

    }
    if(length(lg.abs) > 0) {
      print(
        paste(
          lg.abs,
          "was not found; plots for this gene will be excluded from the final list...",
          sep = " "
        )
      )
    }
  }
  SeuratObject::DefaultAssay(d) <- "RNA"
  sc.umap.panel.gene.list(
    # Gene list
    list.genes[[1]],
    # Seurat Object
    d,
    # Group variable
    "recluster",
    # Color scheme
    col_univ()[1:length(levels(d@meta.data[["recluster"]]))],
    # Color Names
    c(levels(d@meta.data[["recluster"]])),
    # legend x-position
    0.95,
    # legend y-position
    0.95,
    # Run in parallel? (set to FALSE if on Windows)
    parl,
    # Percentage of available cores to use
    core.perc,
    # reduction to plot
    slot2
    )

  return(d)
}
