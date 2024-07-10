#' scRNA-Seq UMAP Plot
#'
#' Generates a panel of UMAPs given a Seurat object containing PCA and UMAP results.
#'
#' @param so An object of class Seurat.
#' @param md.list A vector of character strings indicating metadata columns for overlaying on a loadings plot.
#' @return A series of UMAPs with specified metadata overlays.
#' @examples
#'
#' p.umap <- sc.umap.panel(d.integrated,c("col1","col2","col3"))
#'
#' @export
sc.umap.panel <- function(
  so,
  md.list
  ) {
  d <- so
  d2 <- data.frame(
    d@meta.data,
    `UMAP.1` = d@reductions$umap@cell.embeddings[,1],
    `UMAP.2` = d@reductions$umap@cell.embeddings[,2],
    `UMAP.3` = d@reductions$umap@cell.embeddings[,3]
    )
  # Generate df for each grouping
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
            1,
            0.95)
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
#' @return A series of plots stored as a ggplot2 object for visualizing cluster gene expression.
#' @examples
#'
#' p.umap <- sc.umap.panel(d.integrated,c("col1","col2","col3"))
#'
#' @export
sc.umap.panel.gene <- function(
    so,
    md.var,
    g.name,
    col.scheme,
    col.names,
    leg.x,
    leg.y
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
    `UMAP.1` = d@reductions$umap@cell.embeddings[,1],
    `UMAP.2` = d@reductions$umap@cell.embeddings[,2]
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

    ggplot2::geom_text_repel(data = setNames(
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
#' @return A list of plots saved as ggplot2 objects for visualizing cluster gene expression.
#' @examples
#'
#' sc.umap.panel.gene.list(list.genes,d.seurat,"seurat_clusters",col.vec,col.vec.names,0.95,0.95)
#'
#' @export
sc.umap.panel.gene.list <- function(list.g,so,md.var,g.name,col.scheme,col.names,leg.x,leg.y) {
  lg <- list.g
  d <- so
  lg <- unique(lg[lg %in% SeuratObject::Features(d)])
  lg.abs <- subset(lg, !(lg %in% SeuratObject::Features(d)))
  print(paste(lg.abs, "was not found; plots for this gene will be excluded from the final list...",sep = " "))
  # Create plots
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
          leg.y
          )
      # Save each plot
      ggplot2::ggsave(
        paste("analysis/gene/plot.umap.exp.",x,".png",sep = ""),
        pg,
        height = 24,
        width = 20,
        dpi = 700
        )
      }
    )
  }



















# fun.sc.umap.plot <- function(
#     so,
#     md.list,
#     angle.lab,
#     type.lab,
#     col.scheme
# ) {
#
#   # Format input data
#
#   d <- so
#
#   d2 <- data.frame(
#     d@meta.data,
#     `UMAP.1` = d@reductions$umap@cell.embeddings[,1],
#     `UMAP.2` = d@reductions$umap@cell.embeddings[,2],
#     `UMAP.3` = d@reductions$umap@cell.embeddings[,3]
#   )
#
#   # Generate df for each grouping
#
#   d2.list <- setNames(
#     lapply(
#       c(
#         md.list
#       ),
#       function(x)
#         d2[,c(
#             x,
#             "UMAP.1",
#             "UMAP.2",
#             "UMAP.3"
#           )]
#         ),
#     c(
#       md.list
#       )
#     )
#
#   # Generate plots
#
#   d2.plot <- lapply(
#     c(
#       md.list
#     ),
#     function(x)
#     {
#
#     p <-  ggplot(
#         d2.list[[x]],
#         aes(
#           x=`UMAP.1`,
#           y=`UMAP.2`,
#           z=`UMAP.3`,
#           color = .data[[x]],
#           label = .data[[x]]
#         )
#       ) +
#
#       scale_color_manual(
#         paste(""),
#         values = col.scheme
#       ) +
#
#       # Add 3D points, axes, and axis-labels
#
#       axes_3D(
#         linewidth = 1
#         ) +
#
#       stat_3D(
#         geom = "point",
#         shape = 16,
#         size = 1,
#         alpha = 0.6
#         ) +
#
#       stat_3D(
#         geom = type.lab,
#         data = setNames(
#           aggregate(
#             d2.list[[x]][,c(
#               "UMAP.1",
#               "UMAP.2",
#               "UMAP.3"
#               )],
#             list(
#               d2.list[[x]][[x]]
#               ),
#             FUN = median
#             ),
#           c(
#             x,
#             names(
#               d2.list[[x]][,c(
#                 "UMAP.1",
#                 "UMAP.2",
#                 "UMAP.3"
#                 )]
#               )
#             )
#           ),
#         size = 4,
#         bg.color = "white"
#         ) +
#
#       labs_3D(
#         labs = c(
#           "UMAP-1",
#           "UMAP-2",
#           "UMAP-3"
#           ),
#         vjust = c(
#           -0.2,-0.2,-0.2
#           ),
#         hjust = c(
#           0,1,1
#           ),
#         angle = c(
#           angle.lab,
#           -angle.lab,
#           90
#           )
#         ) +
#
#         # Add general multivariate plot theme and adjust axis text
#         thm.mult +
#
#         theme(
#           panel.grid.major.y = element_blank(),
#           axis.text.x = element_blank(),
#           axis.text.y = element_blank(),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           axis.ticks = element_blank(),
#           plot.margin = unit(
#             c(
#               0.1,0.1,0.1,0.1
#               ),
#             "cm"
#             ),
#           legend.position = "none"
#         )
#
#     }
#   )
#
#
#   # Combine output
#
#   if(
#     length(
#       d2.plot
#     ) >
#     1
#   ) {
#
#     d2.out <- ggarrange(
#       plotlist = d2.plot,
#       labels = c(
#         names(
#           d2.plot
#         )
#       ),
#       ncol = ifelse(
#         length(
#           d2.plot
#         ) <=3,
#         length(
#           d2.plot
#         ),
#         3
#       ),
#       nrow = ifelse(
#         length(
#           d2.plot
#         ) > 3,
#         2,
#         1
#       )
#     )
#
#   }
#
#
#   if(
#     length(
#       d2.plot
#     ) == 1
#   ) {
#
#     d2.out <- d2.plot[[1]]
#
#   }
#
#
#
#
#   return(d2.out)
#
# }



# fun.sc.umap.plot.2D <- function(
#     so,
#     md.list,
#     col.scheme
# ) {
#
#   # Format input data
#
#   d <- so
#
#   d2 <- data.frame(
#     d@meta.data,
#     `UMAP.1` = d@reductions$umap@cell.embeddings[,1],
#     `UMAP.2` = d@reductions$umap@cell.embeddings[,2],
#     `UMAP.3` = d@reductions$umap@cell.embeddings[,3]
#   )
#
#   # Generate df for each grouping
#
#   d2.list <- setNames(
#     lapply(
#       c(
#         md.list
#       ),
#       function(x)
#         d2[,c(
#           x,
#           "UMAP.1",
#           "UMAP.2",
#           "UMAP.3"
#         )]
#     ),
#     c(
#       md.list
#     )
#   )
#
#   # Generate plots
#
#   d2.plot <- lapply(
#     c(
#       md.list
#     ),
#     function(x)
#     {
#
#       p <-  ggplot(
#         d2.list[[x]],
#         aes(
#           x=`UMAP.1`,
#           y=`UMAP.2`,
#           color = .data[[x]],
#           label = .data[[x]]
#         )
#       ) +
#
#         scale_color_manual(
#           paste(""),
#           values = col.scheme
#         ) +
#
#         # Add points and labels
#
#         geom_point(
#           shape = 16,
#           size = 1,
#           alpha = 0.6
#         ) +
#
#         geom_text_repel(
#           data = setNames(
#             aggregate(
#               d2.list[[x]][,c(
#                 "UMAP.1",
#                 "UMAP.2"
#               )],
#               list(
#                 d2.list[[x]][[x]]
#               ),
#               FUN = median
#             ),
#             c(
#               x,
#               names(
#                 d2.list[[x]][,c(
#                   "UMAP.1",
#                   "UMAP.2"
#                 )]
#               )
#             )
#           ),
#           size = 4,
#           bg.color = "white"
#         ) +
#
#
#
#         # Add general multivariate plot theme and adjust axis text
#         thm.mult +
#
#         theme(
#           panel.grid.major.y = element_blank(),
#           axis.text.x = element_blank(),
#           axis.text.y = element_blank(),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           axis.ticks = element_blank(),
#           plot.margin = unit(
#             c(
#               0.1,0.1,0.1,0.1
#             ),
#             "cm"
#           ),
#           legend.position = "none"
#         )
#
#     }
#   )
#
#
#   # Combine output
#
#   if(
#     length(
#       d2.plot
#     ) >
#     1
#   ) {
#
#     d2.out <- ggarrange(
#       plotlist = d2.plot,
#       labels = c(
#         names(
#           d2.plot
#         )
#       ),
#       ncol = ifelse(
#         length(
#           d2.plot
#         ) <=3,
#         length(
#           d2.plot
#         ),
#         3
#       ),
#       nrow = ifelse(
#         length(
#           d2.plot
#         ) > 3,
#         2,
#         1
#       )
#     )
#
#   }
#
#
#   if(
#     length(
#       d2.plot
#     ) == 1
#   ) {
#
#     d2.out <- d2.plot[[1]]
#
#   }
#
#
#
#
#   return(d2.out)
#
# }














# fun.sc.umap.plot2 <- function(
#     so,
#     md.var,
#     g.list,
#     angle.lab,
#     type.lab,
#     col.scheme
# ) {
#
#   # Format input data
#
#   d <- so
#
#   d2 <- data.frame(
#     FetchData(
#       d,
#       vars = c(
#         g.list,
#         md.var
#         )
#     ),
#     `UMAP.1` = d@reductions$umap@cell.embeddings[,1],
#     `UMAP.2` = d@reductions$umap@cell.embeddings[,2],
#     `UMAP.3` = d@reductions$umap@cell.embeddings[,3]
#   )
#
#   # Generate df for each grouping
#
#   d2.list <- setNames(
#     lapply(
#       gsub(
#         "-",
#         ".",
#         g.list
#       ),
#       function(x)
#         d2[,c(
#           x,
#           md.var,
#           "UMAP.1",
#           "UMAP.2",
#           "UMAP.3"
#         )]
#     ),
#     c(
#       g.list
#     )
#   )
#
#   # Generate plots
#
#   d2.plot <- lapply(
#     c(
#       g.list
#     ),
#     function(x)
#     {
#       ## Gene expression overlay
#       p <-  ggplot(
#         d2.list[[x]],
#         aes(
#           x=`UMAP.1`,
#           y=`UMAP.2`,
#           z=`UMAP.3`,
#           color = .data[[gsub(
#             "-",
#             ".",
#             x
#             )]]
#         )
#       ) +
#
#         scale_color_gradientn(
#           name = "Relative Expression",
#           colors = col1b
#         ) +
#
#         # Add 3D points, axes, and axis-labels
#
#         axes_3D(
#           linewidth = 1
#         ) +
#
#         stat_3D(
#           geom = "point",
#           shape = 16,
#           size = 1,
#           alpha = 0.6
#         ) +
#
#         labs_3D(
#           labs = c(
#             "UMAP-1",
#             "UMAP-2",
#             "UMAP-3"
#           ),
#           vjust = c(
#             -0.2,-0.2,-0.2
#           ),
#           hjust = c(
#             0,1,1
#           ),
#           angle = c(
#             angle.lab,
#             -angle.lab,
#             90
#           )
#         ) +
#
#         ggtitle(x) +
#
#         # Add general multivariate plot theme and adjust axis text
#         thm.mult +
#
#         theme(
#           panel.grid.major.y = element_blank(),
#           axis.text.x = element_blank(),
#           axis.text.y = element_blank(),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           axis.ticks = element_blank(),
#           plot.margin = unit(
#             c(
#               0.1,0.1,0.1,0.1
#             ),
#             "cm"
#           ),
#           legend.position = c(
#             0.9,
#             0.9
#           )
#         )
#
#
#
#
#         ## Metadata overlay
#
#         p.md <-  ggplot(
#           d2.list[[x]],
#           aes(
#             x=`UMAP.1`,
#             y=`UMAP.2`,
#             z=`UMAP.3`,
#             color = .data[[md.var]],
#             label = .data[[md.var]]
#           )
#         ) +
#
#           scale_color_manual(
#             paste(""),
#             values = col.scheme
#           ) +
#
#           # Add 3D points, axes, and axis-labels
#
#           axes_3D(
#             linewidth = 1
#           ) +
#
#           stat_3D(
#             geom = "point",
#             shape = 16,
#             size = 1,
#             alpha = 0.6
#           ) +
#
#           stat_3D(
#             geom = type.lab,
#             data = setNames(
#               aggregate(
#                 d2.list[[x]][,c(
#                   "UMAP.1",
#                   "UMAP.2",
#                   "UMAP.3"
#                 )],
#                 list(
#                   d2.list[[x]][[md.var]]
#                 ),
#                 FUN = median
#               ),
#               c(
#                 md.var,
#                 names(
#                   d2.list[[x]][,c(
#                     "UMAP.1",
#                     "UMAP.2",
#                     "UMAP.3"
#                   )]
#                 )
#               )
#             ),
#             size = 4,
#             bg.color = "white"
#           ) +
#
#           labs_3D(
#             labs = c(
#               "UMAP-1",
#               "UMAP-2",
#               "UMAP-3"
#             ),
#             vjust = c(
#               -0.2,-0.2,-0.2
#             ),
#             hjust = c(
#               0,1,1
#             ),
#             angle = c(
#               angle.lab,
#               -angle.lab,
#               90
#             )
#           ) +
#
#           ggtitle(md.var) +
#
#         # Add general multivariate plot theme and adjust axis text
#         thm.mult +
#
#           theme(
#             panel.grid.major.y = element_blank(),
#             axis.text.x = element_blank(),
#             axis.text.y = element_blank(),
#             axis.title.x = element_blank(),
#             axis.title.y = element_blank(),
#             axis.ticks = element_blank(),
#             plot.margin = unit(
#               c(
#                 0.1,0.1,0.1,0.1
#               ),
#               "cm"
#             ),
#             legend.position = "none"
#           )
#
#
#         ## Violin Plot
#
#         plot.v <- ggplot(
#           d2.list[[x]],
#           aes(
#             x = .data[[md.var]],
#             y = .data[[gsub(
#               "-",
#               ".",
#               x
#             )]],
#             fill = .data[[md.var]]
#             )
#           ) +
#
#           scale_fill_manual(
#             name = md.var,
#             values = col.scheme
#             ) +
#
#           # Add violin plot and dotplot
#           geom_violin(
#             trim = T
#             ) +
#
#           geom_jitter(
#             aes(
#               alpha = 0.2
#               ),
#             shape = 16,
#             size = 0.2,
#             position = position_jitter(
#               width = 0.4
#               ),
#             show.legend = F
#             ) +
#
#           # Add Theme
#
#           thm.univ +
#           thm.leg.title.y +
#           labs(
#             y = "Relative Expression"
#             ) +
#           theme(
#             plot.margin = unit(
#               c(
#                 0.1,
#                 0.1,
#                 0.1,
#                 0.1
#                 ),
#               "cm"
#               ),
#             legend.position = "none"
#             )
#
#         return(
#           list(
#             "p.g" = p,
#             "p.d" = p.md,
#             "p.v" = plot.v
#           )
#         )
#
#
#
#     }
#   )
#
#
#   # Combine output
#
#   d2.out <- setNames(
#     lapply(
#       d2.plot,
#       function(x)
#         {
#
#         p <- ggarrange(
#           plotlist = x,
#           ncol = 3,
#           nrow = 1,
#           common.legend = F
#         )
#
#         }
#       ),
#     c(
#       g.list
#       )
#     )
#
#
#
#
#   return(d2.out)
#
# }






# fun.sc.umap.plot.single <- function(
#     so,
#     var.lab,
#     var.col,
#     angle.lab,
#     type.lab,
#     col.scheme
# ) {
#
#   # Format input data
#
#   d <- so
#
#   d2 <- data.frame(
#     d@meta.data,
#     `UMAP.1` = d@reductions$umap@cell.embeddings[,1],
#     `UMAP.2` = d@reductions$umap@cell.embeddings[,2],
#     `UMAP.3` = d@reductions$umap@cell.embeddings[,3]
#   )
#
#   d3 <- dplyr::count(
#     dplyr::group_by(
#       d2,
#       .data[[var.lab]]
#     ),
#     .data[[var.col]]
#   )
#
#   d3 <- dplyr::left_join(
#     setNames(
#       aggregate(
#         d3[["n"]],
#         list(
#           d3[[var.lab]]
#         ),
#         FUN = max
#       ),
#       c(
#         var.lab,
#         "n"
#       )
#     ),
#     d3[,
#        c(
#          var.col,
#          "n"
#        )],
#     by = "n"
#   )
#
#
#   # Generate plots
#
#   d2.plot <- ggplot(
#     d2,
#     aes(
#       x=`UMAP.1`,
#       y=`UMAP.2`,
#       z=`UMAP.3`,
#       color = .data[[var.col]],
#       label = .data[[var.lab]]
#     )
#   ) +
#
#     scale_color_manual(
#       paste(var.col),
#       values = col.scheme
#     ) +
#
#     # Add 3D points, axes, and axis-labels
#
#     axes_3D(
#       linewidth = 1
#     ) +
#
#     stat_3D(
#       geom = "point",
#       shape = 16,
#       size = 1,
#       alpha = 0.6
#     ) +
#
#     stat_3D(
#       geom = type.lab,
#       data = dplyr::left_join(
#         setNames(
#           aggregate(
#             d2[,c(
#               "UMAP.1",
#               "UMAP.2",
#               "UMAP.3"
#             )],
#             list(
#               d2[[var.lab]]
#             ),
#             FUN = median
#           ),
#           c(
#             var.lab,
#             names(
#               d2[,c(
#                 "UMAP.1",
#                 "UMAP.2",
#                 "UMAP.3"
#               )]
#             )
#           )
#         ),
#         d3[,c(
#           var.lab,
#           var.col
#         )],
#         by = var.lab
#       ),
#       size = 4,
#       bg.color = "white"
#     ) +
#
#     labs_3D(
#       labs = c(
#         "UMAP-1",
#         "UMAP-2",
#         "UMAP-3"
#       ),
#       vjust = c(
#         -0.2,-0.2,-0.2
#       ),
#       hjust = c(
#         0,1,1
#       ),
#       angle = c(
#         angle.lab,
#         -angle.lab,
#         90
#       )
#     ) +
#
#     # Add general multivariate plot theme and adjust axis text
#     thm.mult +
#
#     theme(
#       panel.grid.major.y = element_blank(),
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank(),
#       axis.ticks = element_blank(),
#       plot.margin = unit(
#         c(
#           0.1,0.1,0.1,0.1
#         ),
#         "cm"
#       )
#     )
#
#   d2.plot
#
#   }






# fun.sc.umap.plot.single.2D <- function(
#     so,
#     var.lab,
#     var.col,
#     col.scheme
# ) {
#
#   # Format input data
#
#   d <- so
#
#   d2 <- data.frame(
#     d@meta.data,
#     `UMAP.1` = d@reductions$umap@cell.embeddings[,1],
#     `UMAP.2` = d@reductions$umap@cell.embeddings[,2]
#   )
#
#   d3 <- dplyr::count(
#     dplyr::group_by(
#       d2,
#       .data[[var.lab]]
#     ),
#     .data[[var.col]]
#   )
#
#   d3 <- dplyr::left_join(
#     setNames(
#       aggregate(
#         d3[["n"]],
#         list(
#           d3[[var.lab]]
#         ),
#         FUN = max
#       ),
#       c(
#         var.lab,
#         "n"
#       )
#     ),
#     d3[,
#        c(
#          var.col,
#          "n"
#        )],
#     by = "n"
#   )
#
#
#   # Generate plots
#
#   d2.plot <- ggplot(
#     d2,
#     aes(
#       x=`UMAP.1`,
#       y=`UMAP.2`,
#       color = .data[[var.col]],
#       label = .data[[var.lab]]
#     )
#   ) +
#
#     scale_color_manual(
#       paste(var.col),
#       values = col.scheme
#     ) +
#
#     geom_point(
#       shape = 16,
#       size = 1,
#       alpha = 0.6
#     ) +
#
#     geom_text_repel(
#       data = dplyr::left_join(
#         setNames(
#           aggregate(
#             d2[,c(
#               "UMAP.1",
#               "UMAP.2"
#               )],
#             list(
#               d2[[var.lab]]
#             ),
#             FUN = median
#           ),
#           c(
#             var.lab,
#             names(
#               d2[,c(
#                 "UMAP.1",
#                 "UMAP.2"
#                 )]
#             )
#           )
#         ),
#         d3[,c(
#           var.lab,
#           var.col
#         )],
#         by = var.lab
#         ),
#       size = 4,
#       bg.color = "white"
#       ) +
#
#     # Add general multivariate plot theme and adjust axis text
#
#     thm.mult +
#
#     theme(
#       panel.grid.major.y = element_blank(),
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank(),
#       axis.ticks = element_blank(),
#       plot.margin = unit(
#         c(
#           0.1,0.1,0.1,0.1
#         ),
#         "cm"
#       )
#     )
#
#   d2.plot
#
# }











# fun.sc.umap.plot.panel2 <- function(
#     so,
#     md.var,
#     g.name,
#     angle.lab,
#     type.lab,
#     col.scheme,
#     col.names
# ) {
#
#   # Format input data
#   cols <- setNames(col.scheme,
#                    col.names)
#
#
#   d <- so
#
#   d2 <- data.frame(
#     FetchData(
#       d,
#       vars = c(
#         gsub(
#           "-",
#           "-",
#           g.name
#         ),
#         md.var
#       )
#     ),
#     `UMAP.1` = d@reductions$umap@cell.embeddings[,1],
#     `UMAP.2` = d@reductions$umap@cell.embeddings[,2],
#     `UMAP.3` = d@reductions$umap@cell.embeddings[,3]
#   )
#
#   # Generate plots
#
#   d2.plot <- ggplot(
#         d2,
#         aes(
#           x=`UMAP.1`,
#           y=`UMAP.2`,
#           z=`UMAP.3`,
#           color = .data[[gsub(
#             "-",
#             ".",
#             g.name
#           )]]
#         )
#       ) +
#
#         scale_color_gradientn(
#           name = "Relative Expression",
#           colors = col1b
#         ) +
#
#         # Add 3D points, axes, and axis-labels
#
#         axes_3D(
#           linewidth = 1
#         ) +
#
#         stat_3D(
#           geom = "point",
#           shape = 16,
#           size = 1,
#           alpha = 0.6
#         ) +
#
#         labs_3D(
#           labs = c(
#             "UMAP-1",
#             "UMAP-2",
#             "UMAP-3"
#           ),
#           vjust = c(
#             -0.2,-0.2,-0.2
#           ),
#           hjust = c(
#             0,1,1
#           ),
#           angle = c(
#             angle.lab,
#             -angle.lab,
#             90
#           )
#         ) +
#
#         ggtitle(g.name) +
#
#         # Add general multivariate plot theme and adjust axis text
#         thm.mult +
#
#         theme(
#           panel.grid.major.y = element_blank(),
#           axis.text.x = element_blank(),
#           axis.text.y = element_blank(),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           axis.ticks = element_blank(),
#           plot.margin = unit(
#             c(
#               0.1,0.1,0.1,0.1
#             ),
#             "cm"
#           ),
#           legend.position = c(
#             0.9,
#             0.9
#           )
#         )
#
#
#
#
#       ## Metadata overlay
#
#       p.md <-  ggplot(
#         d2,
#         aes(
#           x=`UMAP.1`,
#           y=`UMAP.2`,
#           z=`UMAP.3`,
#           color = .data[[md.var]],
#           label = .data[[md.var]]
#         )
#       ) +
#
#         scale_color_manual(
#           paste(""),
#           values = cols
#         ) +
#
#         # Add 3D points, axes, and axis-labels
#
#         axes_3D(
#           linewidth = 1
#         ) +
#
#         stat_3D(
#           geom = "point",
#           shape = 16,
#           size = 1,
#           alpha = 0.6
#         ) +
#
#         stat_3D(
#           geom = type.lab,
#           data = setNames(
#             aggregate(
#               d2[,c(
#                 "UMAP.1",
#                 "UMAP.2",
#                 "UMAP.3"
#               )],
#               list(
#                 d2[[md.var]]
#               ),
#               FUN = median
#             ),
#             c(
#               md.var,
#               names(
#                 d2[,c(
#                   "UMAP.1",
#                   "UMAP.2",
#                   "UMAP.3"
#                 )]
#               )
#             )
#           ),
#           size = 4,
#           bg.color = "white"
#         ) +
#
#         labs_3D(
#           labs = c(
#             "UMAP-1",
#             "UMAP-2",
#             "UMAP-3"
#           ),
#           vjust = c(
#             -0.2,-0.2,-0.2
#           ),
#           hjust = c(
#             0,1,1
#           ),
#           angle = c(
#             angle.lab,
#             -angle.lab,
#             90
#           )
#         ) +
#
#         ggtitle(md.var) +
#
#         # Add general multivariate plot theme and adjust axis text
#         thm.mult +
#
#         theme(
#           panel.grid.major.y = element_blank(),
#           axis.text.x = element_blank(),
#           axis.text.y = element_blank(),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           axis.ticks = element_blank(),
#           plot.margin = unit(
#             c(
#               0.1,0.1,0.1,0.1
#             ),
#             "cm"
#           ),
#           legend.position = "none"
#         )
#
#
#       ## Violin Plot
#
#       plot.v <- ggplot(
#         d2,
#         aes(
#           x = .data[[md.var]],
#           y = .data[[gsub(
#             "-",
#             ".",
#             g.name
#           )]],
#           fill = .data[[md.var]]
#         )
#       ) +
#
#         scale_fill_manual(
#           name = md.var,
#           values = col.scheme
#         ) +
#
#         # Add violin plot and dotplot
#         geom_violin(
#           trim = T
#         ) +
#
#         geom_jitter(
#           aes(
#             alpha = 0.2
#           ),
#           shape = 16,
#           size = 0.2,
#           position = position_jitter(
#             width = 0.4
#           ),
#           show.legend = F
#         ) +
#
#         # Add Theme
#
#         thm.univ +
#         thm.leg.title.y +
#         labs(
#           y = "Relative Expression"
#         ) +
#         theme(
#           plot.margin = unit(
#             c(
#               0.1,
#               0.1,
#               0.1,
#               0.1
#             ),
#             "cm"
#           ),
#           legend.position = "none"
#         )
#
#
#   # Combine output
#
#   d2.out <- ggarrange(
#     d2.plot,
#     p.md,
#     plot.v,
#     ncol = 3,
#     nrow = 1,
#     common.legend = F
#     )
#
#
#   return(d2.out)
#
#   }






# fun.sc.umap.plot.panel2.2D <- function(
#     so,
#     md.var,
#     g.name,
#     col.scheme,
#     col.names,
#     leg.x,
#     leg.y
# ) {
#
#   # Format input data
#   cols <- setNames(col.scheme,
#                    col.names)
#
#
#   d <- so
#
#   d2 <- data.frame(
#     FetchData(
#       d,
#       vars = c(
#         gsub(
#           "-",
#           "-",
#           g.name
#         ),
#         md.var
#       )
#     ),
#     `UMAP.1` = d@reductions$umap@cell.embeddings[,1],
#     `UMAP.2` = d@reductions$umap@cell.embeddings[,2]
#     )
#
#   # Generate plots
#
#   d2.plot <- ggplot(
#     d2,
#     aes(
#       x=`UMAP.1`,
#       y=`UMAP.2`,
#       color = .data[[gsub(
#         "-",
#         ".",
#         g.name
#       )]]
#     )
#   ) +
#
#     scale_color_gradientn(
#       name = "Relative Expression",
#       colors = col1b
#     ) +
#
#     # Add 3D points, axes, and axis-labels
#
#     geom_point(shape = 16,
#                size = 1,
#                alpha = 0.6) +
#
#     ggtitle(g.name) +
#
#     # Add general multivariate plot theme and adjust axis text
#     thm.mult +
#
#     theme(
#       panel.grid.major.y = element_blank(),
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank(),
#       axis.ticks = element_blank(),
#       plot.margin = unit(
#         c(
#           0.1,0.1,0.1,0.1
#         ),
#         "cm"
#       ),
#       legend.position = c(
#         leg.x,
#         leg.y
#       )
#     )
#
#
#
#
#   ## Metadata overlay
#
#   p.md <-  ggplot(
#     d2,
#     aes(
#       x=`UMAP.1`,
#       y=`UMAP.2`,
#       color = .data[[md.var]],
#       label = .data[[md.var]]
#     )
#   ) +
#
#     scale_color_manual(
#       paste(""),
#       values = cols
#     ) +
#
#     geom_point(shape = 16,
#                size = 1,
#                alpha = 0.6) +
#
#     geom_text_repel(data = setNames(
#       aggregate(
#         d2[,c(
#           "UMAP.1",
#           "UMAP.2"
#         )],
#         list(
#           d2[[md.var]]
#         ),
#         FUN = median
#       ),
#       c(
#         md.var,
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
#     ggtitle(md.var) +
#
#     # Add general multivariate plot theme and adjust axis text
#     thm.mult +
#
#     theme(
#       panel.grid.major.y = element_blank(),
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank(),
#       axis.ticks = element_blank(),
#       plot.margin = unit(
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
#   plot.v <- ggplot(
#     d2,
#     aes(
#       x = .data[[md.var]],
#       y = .data[[gsub(
#         "-",
#         ".",
#         g.name
#       )]],
#       fill = .data[[md.var]]
#     )
#   ) +
#
#     scale_fill_manual(
#       name = md.var,
#       values = col.scheme
#     ) +
#
#     # Add violin plot and dotplot
#     geom_violin(
#       trim = T
#     ) +
#
#     geom_jitter(
#       aes(
#         alpha = 0.2
#       ),
#       shape = 16,
#       size = 0.2,
#       position = position_jitter(
#         width = 0.4
#       ),
#       show.legend = F
#     ) +
#
#     # Add Theme
#
#     thm.univ +
#     thm.leg.title.y +
#     labs(
#       y = "Relative Expression"
#     ) +
#     theme(
#       plot.margin = unit(
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
#
#   # Combine output
#
#   d2.out <- ggarrange(
#     d2.plot,
#     p.md,
#     plot.v,
#     ncol = 3,
#     nrow = 1,
#     common.legend = F
#   )
#
#
#   return(d2.out)
#
# }














# fun.sc.umap.plot.traj.2D <- function(
#     cds.object,
#     var.lab,
#     var.col,
#     col.scheme
# ) {
#
#   # Format input data
#
#   d <- cds.object
#
#   d2 <- setNames(data.frame(
#     d@colData,
#     d@reduce_dim_aux$UMAP@listData$model$umap_model$embedding
#   ),
#   c(names(d@colData),
#     "UMAP.1",
#     "UMAP.2",
#     "UMAP.3"))
#
#   d3 <- dplyr::count(
#     dplyr::group_by(
#       d2,
#       .data[[var.lab]]
#     ),
#     .data[[var.col]]
#   )
#
#   d3 <- dplyr::left_join(
#     setNames(
#       aggregate(
#         d3[["n"]],
#         list(
#           d3[[var.lab]]
#         ),
#         FUN = max
#       ),
#       c(
#         var.lab,
#         "n"
#       )
#     ),
#     d3[,
#        c(
#          var.col,
#          "n"
#        )],
#     by = "n"
#   )
#
#
#   # Generate plots
#
#   d2.plot <- ggplot(
#     d2,
#     aes(
#       x=`UMAP.1`,
#       y=`UMAP.2`,
#       color = .data[[var.col]],
#       label = .data[[var.lab]]
#     )
#   ) +
#
#     scale_color_manual(
#       paste(var.col),
#       values = col.scheme
#     ) +
#
#     geom_point(
#       shape = 16,
#       size = 1,
#       alpha = 0.6
#     ) +
#
#     geom_text_repel(
#       data = dplyr::left_join(
#         setNames(
#           aggregate(
#             d2[,c(
#               "UMAP.1",
#               "UMAP.2"
#             )],
#             list(
#               d2[[var.lab]]
#             ),
#             FUN = median
#           ),
#           c(
#             var.lab,
#             names(
#               d2[,c(
#                 "UMAP.1",
#                 "UMAP.2"
#               )]
#             )
#           )
#         ),
#         d3[,c(
#           var.lab,
#           var.col
#         )],
#         by = var.lab
#       ),
#       size = 4,
#       bg.color = "white"
#     ) +
#
#     # Add general multivariate plot theme and adjust axis text
#
#     thm.mult +
#
#     theme(
#       panel.grid.major.y = element_blank(),
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank(),
#       axis.ticks = element_blank(),
#       plot.margin = unit(
#         c(
#           0.1,0.1,0.1,0.1
#         ),
#         "cm"
#       )
#     )
#
#   d2.plot
#
# }





# ## 3D
#
# fun.sc.umap.plot.traj <- function(
#     cds.object,
#     var.lab,
#     var.col,
#     angle.lab,
#     type.lab,
#     col.scheme
# ) {
#
#   # Format input data
#
#   d <- cds.object
#
#   d2 <- setNames(data.frame(
#     d@colData,
#     d@reduce_dim_aux$UMAP@listData$model$umap_model$embedding
#   ),
#   c(names(d@colData),
#     "UMAP.1",
#     "UMAP.2",
#     "UMAP.3"))
#
#   d3 <- dplyr::count(
#     dplyr::group_by(
#       d2,
#       .data[[var.lab]]
#     ),
#     .data[[var.col]]
#   )
#
#   d3 <- dplyr::left_join(
#     setNames(
#       aggregate(
#         d3[["n"]],
#         list(
#           d3[[var.lab]]
#         ),
#         FUN = max
#       ),
#       c(
#         var.lab,
#         "n"
#       )
#     ),
#     d3[,
#        c(
#          var.col,
#          "n"
#        )],
#     by = "n"
#   )
#
#
#   # Generate plots
#
#   d2.plot <- ggplot(
#     d2,
#     aes(
#       x=`UMAP.1`,
#       y=`UMAP.2`,
#       z=`UMAP.3`,
#       color = .data[[var.col]],
#       label = .data[[var.lab]]
#     )
#   ) +
#
#     scale_color_manual(
#       paste(var.col),
#       values = col.scheme
#     ) +
#
#     # Add 3D points, axes, and axis-labels
#
#     axes_3D(
#       linewidth = 1
#     ) +
#
#     stat_3D(
#       geom = "point",
#       shape = 16,
#       size = 1,
#       alpha = 0.6
#     ) +
#
#     stat_3D(
#       geom = type.lab,
#       data = dplyr::left_join(
#         setNames(
#           aggregate(
#             d2[,c(
#               "UMAP.1",
#               "UMAP.2",
#               "UMAP.3"
#             )],
#             list(
#               d2[[var.lab]]
#             ),
#             FUN = median
#           ),
#           c(
#             var.lab,
#             names(
#               d2[,c(
#                 "UMAP.1",
#                 "UMAP.2",
#                 "UMAP.3"
#               )]
#             )
#           )
#         ),
#         d3[,c(
#           var.lab,
#           var.col
#         )],
#         by = var.lab
#       ),
#       size = 4,
#       bg.color = "white"
#     ) +
#
#     labs_3D(
#       labs = c(
#         "UMAP-1",
#         "UMAP-2",
#         "UMAP-3"
#       ),
#       vjust = c(
#         -0.2,-0.2,-0.2
#       ),
#       hjust = c(
#         0,1,1
#       ),
#       angle = c(
#         angle.lab,
#         -angle.lab,
#         90
#       )
#     ) +
#
#     # Add general multivariate plot theme and adjust axis text
#     thm.mult +
#
#     theme(
#       panel.grid.major.y = element_blank(),
#       axis.text.x = element_blank(),
#       axis.text.y = element_blank(),
#       axis.title.x = element_blank(),
#       axis.title.y = element_blank(),
#       axis.ticks = element_blank(),
#       plot.margin = unit(
#         c(
#           0.1,0.1,0.1,0.1
#         ),
#         "cm"
#       )
#     )
#
#   d2.plot
#
# }
















