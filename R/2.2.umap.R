#' scRNA-Seq UMAP Plot
#'
#' Generates a panel of UMAPs given a Seurat object containing PCA and UMAP results.
#'
#' @param so An object of class Seurat.
#' @param md.list A vector of character strings indicating metadata columns for overlaying on a loadings plot.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A series of UMAPs with specified metadata overlays.
#' @examples
#'
#' # p.umap <- sc.umap.panel(d.integrated,c("col1","col2","col3"),"wnn.umap")
#'
#' @export
sc.umap.panel <- function(
  so,
  md.list,
  slot1
  ) {
  d <- so
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 3){
  d2 <- data.frame(
    d@meta.data,
    `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[,1],
    `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[,2],
    `UMAP.3` = d@reductions[[slot1]]@cell.embeddings[,3]
  )
  d2.list <- setNames(
    lapply(
      c(md.list),
      function(x)
        d2[,c(
          x,
          "UMAP.1",
          "UMAP.2",
          "UMAP.3"
        )]
    ),
    c(md.list)
  )
  
  }
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 2){
    d2 <- data.frame(
      d@meta.data,
      `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[,1],
      `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[,2]
    )
    d2.list <- setNames(
      lapply(
        c(md.list),
        function(x)
          d2[,c(
            x,
            "UMAP.1",
            "UMAP.2"
          )]
      ),
      c(md.list)
    )
    
    }
  # Generate plots
  d2.plot <- lapply(
    c(md.list),
    function(x){
      p <-  ggplot2::ggplot(
        d2.list[[x]],
        ggplot2::aes(
          x=`UMAP.1`,
          y=`UMAP.2`,
          color = .data[[x]],
          label = .data[[x]]
          )
        ) +
        ggplot2::scale_color_manual(
          paste(""),
          values = col.univ()
          ) +
        # Add points
        ggplot2::geom_point(
          shape = 16,
          size = 1,
          alpha = 0.6
          ) +
        sc.theme1() +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          plot.margin = ggplot2::unit(
            c(0.1,0.1,0.1,0.1),
            "cm"
            ),
          legend.position = c(
            0.9,
            0.85)
          )
      }
    )
  # Combine output
  if(length(d2.plot) > 1) {
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
    }
  if(
    length(
      d2.plot
      ) == 1
    ) {
  d2.out <- d2.plot[[1]]
    }
  return(d2.out)
  }




#' Visualize Individual Gene Expression
#'
#' Generates a series of plots to visualize individual gene expression per cluster.
#'
#' @param so An object of class Seurat.
#' @param md.var A character string indicating the clustering column for overlaying on a UMAP plot.
#' @param g.name A character string indicating the gene name.
#' @param col.scheme The color scheme to be used for distinguishing between groups, provided as a vector.
#' @param col.names A vector of the same length as the provided color scheme for assigning colors to each group.
#' @param leg.x A numeric value indicating the placement of the figure legend on the x-axis.
#' @param leg.y A numeric value indicating the placement of the figure legend on the y-axis.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A series of plots stored as a ggplot2 object for visualizing cluster gene expression.
#' @examples
#'
#' # p.umap <- sc.umap.panel.gene(
#' # d.integrated,c("col1","col2","col3"),"CFTR",col.univ,c("group1","group2"),0.95,0.95,"wnn.umap")
#'
#' @export
sc.umap.panel.gene <- function(
    so,
    md.var,
    g.name,
    col.scheme,
    col.names,
    leg.x,
    leg.y,
    slot1
) {
  # Format input data
  cols <- setNames(col.scheme,
                   col.names)
  d <- so
  d2 <- data.frame(
    Seurat::FetchData(
      d,
      vars = c(
        gsub(
          "-",
          "-",
          g.name
        ),
        md.var
      )
    ),
    `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[,1],
    `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[,2]
    )

  # Generate plots
  d2.plot <- ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x=`UMAP.1`,
      y=`UMAP.2`,
      color = .data[[gsub(
        "-",
        ".",
        g.name
      )]]
    )
  ) +

    ggplot2::scale_color_gradientn(
      name = "Relative Expression",
      colors = col.grad()
    ) +

    # Add 3D points, axes, and axis-labels

    ggplot2::geom_point(shape = 16,
               size = 1,
               alpha = 0.6) +

    ggplot2::ggtitle(g.name) +

    # Add general multivariate plot theme and adjust axis text
    sc.theme1() +

    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(
        c(
          0.1,0.1,0.1,0.1
        ),
        "cm"
      ),
      legend.position = c(
        leg.x,
        leg.y
      )
    )




  ## Metadata overlay

  p.md <-  ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x=`UMAP.1`,
      y=`UMAP.2`,
      color = .data[[md.var]],
      label = .data[[md.var]]
    )
  ) +

    ggplot2::scale_color_manual(
      paste(""),
      values = cols
    ) +

    ggplot2::geom_point(shape = 16,
               size = 1,
               alpha = 0.6) +

    ggrepel::geom_text_repel(data = setNames(
      aggregate(
        d2[,c(
          "UMAP.1",
          "UMAP.2"
        )],
        list(
          d2[[md.var]]
        ),
        FUN = median
      ),
      c(
        md.var,
        names(
          d2[,c(
            "UMAP.1",
            "UMAP.2"
          )]
        )
      )
    ),
    size = 4,
    bg.color = "white") +

    ggplot2::ggtitle(md.var) +

    # Add general multivariate plot theme and adjust axis text
    sc.theme1() +

    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(
        c(
          0.1,0.1,0.1,0.1
        ),
        "cm"
      ),
      legend.position = "none"
    )


  ## Violin Plot

  plot.v <- ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x = .data[[md.var]],
      y = .data[[gsub(
        "-",
        ".",
        g.name
      )]],
      fill = .data[[md.var]]
    )
  ) +

    ggplot2::scale_fill_manual(
      name = md.var,
      values = col.scheme
    ) +

    # Add violin plot and dotplot
    ggplot2::geom_violin(
      trim = T
    ) +

    ggplot2::geom_jitter(
      ggplot2::aes(
        alpha = 0.2
      ),
      shape = 16,
      size = 0.2,
      position = ggplot2::position_jitter(
        width = 0.4
      ),
      show.legend = F
    ) +

    # Add Theme

    sc.theme1() +
    ggplot2::labs(
      y = "Relative Expression"
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::unit(
        c(
          0.1,
          0.1,
          0.1,
          0.1
        ),
        "cm"
      ),
      legend.position = "none"
    )

  # Combine output
  d2.out <- ggpubr::ggarrange(
    d2.plot,
    p.md,
    plot.v,
    ncol = 3,
    nrow = 1,
    common.legend = F
  )

  return(d2.out)

  }




#' Visualize Gene List Expression
#'
#' Generates a series of plots to visualize the expression of a list of genes per cluster.
#'
#' @param list.g A vector containing a list of genes to be plotted for a given Seurat object.
#' @param so An object of class Seurat.
#' @param md.var A character string indicating the clustering column for overlaying on a UMAP plot.
#' @param col.scheme The color scheme to be used for distinguishing between groups, provided as a vector.
#' @param col.names A vector of the same length as the provided color scheme for assigning colors to each group.
#' @param leg.x A numeric value indicating the placement of the figure legend on the x-axis.
#' @param leg.y A numeric value indicating the placement of the figure legend on the y-axis.
#' @param parl Logical indicating whether processing should be run in parallel (Linux and WSL2 only). Set to FALSE if running sequentially.
#' @param core.perc Percentage of available cores to use if running in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A list of plots saved as ggplot2 objects for visualizing cluster gene expression.
#' @examples
#'
#' # sc.umap.panel.gene.list(
#' # list.genes,
#' # d.seurat,
#' # "seurat_clusters",
#' # col.vec,
#' # col.vec.names,
#' # 0.95,
#' # 0.95,
#' # TRUE,
#' # 0.5,
#' # "umap"
#' # )
#'
#' @export
sc.umap.panel.gene.list <- function(
    list.g,so,md.var,col.scheme,col.names,leg.x,leg.y,parl,core.perc,slot1
    ) {
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
        pg <- sc.umap.panel.gene(
          d,
          md.var,
          x,
          col.scheme,
          col.names,
          leg.x,
          leg.y,
          slot1
        )
        # Save each plot
        ggplot2::ggsave(
          paste("analysis/gene/plot.umap.exp.",x,".png",sep = ""),
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
        pg <- sc.umap.panel.gene(
          d,
          md.var,
          x,
          col.scheme,
          col.names,
          leg.x,
          leg.y,
          slot1
          )
        # Save each plot
        ggplot2::ggsave(
          paste("analysis/gene/plot.umap.exp.",x,".png",sep = ""),
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
        pg <- sc.umap.panel.gene(
          d,
          md.var,
          x,
          col.scheme,
          col.names,
          leg.x,
          leg.y,
          slot1
        )
        # Save each plot
        ggplot2::ggsave(
          paste("analysis/gene/plot.umap.exp.",x,".png",sep = ""),
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






#' Standard UMAP Plot
#'
#' Generates a single UMAP plot with a specific metadata overlay for a Seurat object.
#'
#' @param so An object of class Seurat.
#' @param md.var A character string indicating the clustering column for overlaying on a UMAP plot.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A UMAP plot with points grouped by a specific metadata column.
#' @examples
#'
#' # p.umap <- sc.umap.standard(d.integrated,"col1","wnn.umap")
#'
#' @export
sc.umap.standard <- function(
  so,
  md.var,
  slot1
  ) {
  # Format input data
  d <- so
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 3){
  d2 <- data.frame(
    d@meta.data,
    `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[,1],
    `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[,2],
    `UMAP.3` = d@reductions[[slot1]]@cell.embeddings[,3],
    md.var = d@meta.data[[md.var]]
  )}
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 2){
    d2 <- data.frame(
      d@meta.data,
      `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[,1],
      `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[,2],
      md.var = d@meta.data[[md.var]]
    )}
  
  # Generate plot
  d2.plot <- ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x=`UMAP.1`,
      y=`UMAP.2`,
      color = .data[[md.var]],
      label = .data[[md.var]]
      )
    ) +
    ggplot2::geom_point(
      shape = 16,
      size = 1,
      alpha = 0.6
      ) +
    ggrepel::geom_text_repel(
      data = setNames(
        aggregate(
          d2[,c("UMAP.1","UMAP.2")],
          list(d2[[md.var]]),
          FUN = median
          ),
        c(md.var,names(
          d2[,c("UMAP.1","UMAP.2")])
          )
        ),
      size = 4,
      bg.color = "white"
      ) +
    ggplot2::scale_color_manual(
        paste(""),
        values = col.univ()
        ) +
    sc.theme1() +
    ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          plot.margin = ggplot2::unit(
            c(0.1,0.1,0.1,0.1),"cm"),
          legend.position = "none"
          )
    return(d2.plot)
    }




# 
# sc.umap.panel.gene(
#   # Seurat Object
#   d,
#   # Group variable
#   "CellType",
#   # Gene name
#   "CFTR",
#   # Color scheme
#   col.univ()[1:length(levels(d@meta.data[["CellType"]]))],
#   # Color Names
#   c(levels(d@meta.data[["CellType"]])),
#   # legend x-position
#   0.95,
#   # legend y-position
#   0.95,
#   # reduction to plot
#   "umap_harmony"
#   )
# 
# sc.umap.panel.gene <- function(
#     so,
#     md.var,
#     g.name,
#     col.scheme,
#     col.names,
#     leg.x,
#     leg.y,
#     slot1
# ) {
#   # Format input data
#   cols <- setNames(col.univ()[1:length(levels(d@meta.data[["CellType"]]))],
#                    c(levels(d@meta.data[["CellType"]])))
#   d <- d
#   d2 <- data.frame(
#     Seurat::FetchData(
#       d,
#       vars = c(
#         diff.output.activity[diff.output.activity[["gene"]] == "GRHL1","ID"][[1]],
#         "CellType"
#       )
#     ),
#     `UMAP.1` = d@reductions[["wnn.umap.cor"]]@cell.embeddings[,1],
#     `UMAP.2` = d@reductions[["wnn.umap.cor"]]@cell.embeddings[,2]
#     )
# 
#   # Generate plots
#   d2.plot <- ggplot2::ggplot(
#     d2,
#     ggplot2::aes(
#       x=`UMAP.1`,
#       y=`UMAP.2`,
#       color = .data[[diff.output.activity[diff.output.activity[["gene"]] == "GRHL1","ID"][[1]]]]
#     )
#   ) +
# 
#     ggplot2::scale_color_gradientn(
#       name = "Motif Activity",
#       colors = col.grad()
#     ) +
# 
#     # Add 3D points, axes, and axis-labels
# 
#     ggplot2::geom_point(shape = 16,
#                size = 1,
#                alpha = 0.6) +
# 
#     ggplot2::ggtitle("GRHL1") +
# 
#     # Add general multivariate plot theme and adjust axis text
#     sc.theme1() +
# 
#     ggplot2::theme(
#       panel.grid.major.y = ggplot2::element_blank(),
#       axis.text.x = ggplot2::element_blank(),
#       axis.text.y = ggplot2::element_blank(),
#       axis.title.x = ggplot2::element_blank(),
#       axis.title.y = ggplot2::element_blank(),
#       axis.ticks = ggplot2::element_blank(),
#       plot.margin = ggplot2::unit(
#         c(
#           0.1,0.1,0.1,0.1
#         ),
#         "cm"
#       ),
#       legend.position = c(
#         0.95,
#         0.95
#       )
#     )
# 
# 
#   ## Metadata overlay
# 
#   p.md <-  ggplot2::ggplot(
#     d2,
#     ggplot2::aes(
#       x=`UMAP.1`,
#       y=`UMAP.2`,
#       color = .data[["CellType"]],
#       label = .data[["CellType"]]
#     )
#   ) +
# 
#     ggplot2::scale_color_manual(
#       paste(""),
#       values = cols
#     ) +
# 
#     ggplot2::geom_point(shape = 16,
#                size = 1,
#                alpha = 0.6) +
# 
#     ggrepel::geom_text_repel(data = setNames(
#       aggregate(
#         d2[,c(
#           "UMAP.1",
#           "UMAP.2"
#         )],
#         list(
#           d2[["CellType"]]
#         ),
#         FUN = median
#       ),
#       c(
#         "CellType",
#         names(
#           d2[,c(
#             "UMAP.1",
#             "UMAP.2"
#           )]
#         )
#       )
#     ),
#     size = 4,
#     bg.color = "white") +
# 
#     ggplot2::ggtitle("CellType") +
# 
#     # Add general multivariate plot theme and adjust axis text
#     sc.theme1() +
# 
#     ggplot2::theme(
#       panel.grid.major.y = ggplot2::element_blank(),
#       axis.text.x = ggplot2::element_blank(),
#       axis.text.y = ggplot2::element_blank(),
#       axis.title.x = ggplot2::element_blank(),
#       axis.title.y = ggplot2::element_blank(),
#       axis.ticks = ggplot2::element_blank(),
#       plot.margin = ggplot2::unit(
#         c(
#           0.1,0.1,0.1,0.1
#         ),
#         "cm"
#       ),
#       legend.position = "none"
#     )
# 
# 
#   ## Violin Plot
# 
#   plot.v <- ggplot2::ggplot(
#     d2,
#     ggplot2::aes(
#       x = .data[["CellType"]],
#       y = .data[[diff.output.activity[diff.output.activity[["gene"]] == "GRHL1","ID"][[1]]]],
#       fill = .data[["CellType"]]
#     )
#   ) +
# 
#     ggplot2::scale_fill_manual(
#       name = "CellType",
#       values = col.univ()[1:length(levels(d@meta.data[["CellType"]]))]
#     ) +
# 
#     # Add violin plot and dotplot
#     ggplot2::geom_violin(
#       trim = T
#     ) +
# 
#     ggplot2::geom_jitter(
#       ggplot2::aes(
#         alpha = 0.2
#       ),
#       shape = 16,
#       size = 0.2,
#       position = ggplot2::position_jitter(
#         width = 0.4
#       ),
#       show.legend = F
#     ) +
# 
#     # Add Theme
# 
#     sc.theme1() +
#     ggplot2::labs(
#       y = "Motif Activity"
#     ) +
#     ggplot2::theme(
#       plot.margin = ggplot2::unit(
#         c(
#           0.1,
#           0.1,
#           0.1,
#           0.1
#         ),
#         "cm"
#       ),
#       legend.position = "none"
#     )
# 
#   # Combine output
#   d2.out <- ggpubr::ggarrange(
#     d2.plot,
#     p.md,
#     plot.v,
#     ncol = 3,
#     nrow = 1,
#     common.legend = F
#   )
# d2.out
#   return(d2.out)
# 
#   }
# 
# ggplot2::ggsave(
#           paste("analysis/","test.activity.plot",".png",sep = ""),
#           d2.out,
#           height = 12,
#           width = 36,
#           dpi = 700
#           )