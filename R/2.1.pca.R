#' Seurat PCA
#'
#' Performs PCA on a Seurat object.
#'
#' @param so An object on class Seurat.
#' @return A Seurat object containing PCA results.
#' @examples
#'
#' # d.integrated <- sc.pca(d.integrated)
#'
#' @export
sc.pca <- function(
    so
    ) {
  d <- so
  # Change assay, scale, and run PCA
  Seurat::DefaultAssay(
    object = d
    ) <- "integrated"
  d <- Seurat::RunPCA(
    object = Seurat::ScaleData(
      object = d,
      verbose = T
      ),
    verbose = T
    )
  # Visualize dim loadings and elbow plot
  Seurat::VizDimLoadings(
    object = d,
    dims = 1:2,
    reduction = "pca")
  Seurat::ElbowPlot(
    d,
    reduction = "pca",
    ndims = 50)
  return(d)
  }





#' scRNA-Seq PCA Loadings Plot
#'
#' Generates a panel of loadings plots given a Seurat object containing PCA results.
#'
#' @param so An object of class Seurat.
#' @param md.list A vector of character strings indicating metadata columns for overlaying on a loadings plot.
#' @return A series of loadings plots with specified metadata overlays.
#' @examples
#'
#' # p.pca <- sc.pca.plot(d.integrated,c("col1","col2","col3"))
#'
#' @export
sc.pca.plot <- function(
    so,
    md.list
    ) {
  # Format input data
  d <- so
  d2 <- data.frame(
    d@meta.data,
    PC1 = d@reductions$pca@cell.embeddings[,1],
    PC2 = d@reductions$pca@cell.embeddings[,2]
  )
  # Generate df for each grouping
  d2.list <- setNames(
    lapply(
      c(
        md.list
      ),
      function(x)
        d2[,c(
          x,
          "PC1",
          "PC2"
          )]
    ),
    c(
      md.list
    )
  )
  # Generate plots
  d2.plot <- lapply(
    c(
      md.list
    ),
    function(x)
    {
      ggplot2::ggplot(
        d2.list[[x]],
        ggplot2::aes(
          x=PC1,
          y=PC2,
          color = .data[[x]]
          )
        ) +
        ggplot2::geom_point(
        shape=16,
        size = 2,
        alpha = 0.5
        ) +
        ggplot2::scale_color_manual(
        paste(""),
        values = col.univ()
        ) +
      sc.theme1()
      }
    )
  # Combine output
  d2.out <- ggpubr::ggarrange(
    plotlist = d2.plot,
    labels = c(
      names(
        d2.plot
        )
      ),
    ncol = ifelse(
      length(
        d2.plot
      ) <=3,
      length(
        d2.plot
      ),
      3
    ),
    nrow = ifelse(
      length(
        d2.plot
      ) > 3,
      ceiling(
        length(
          d2.plot
          )/
          3
        ),
      1
      )
    )

  return(d2.out)
  }




