#' scRNA-Seq UMAP Plot
#'
#' Generates a panel of UMAPs given a
#' Seurat object containing PCA and UMAP results.
#'
#' @param so An object of class Seurat.
#' @param md_list A vector of character strings indicating
#' metadata columns for overlaying on a loadings plot.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A series of UMAPs with specified metadata overlays.
#' @examples
#'
#' # p_umap <- sc_umap_panel(d_integrated,c("col1","col2","col3"),"wnn.umap")
#'
#' @export
sc_umap_panel <- function(
  so,
  md_list,
  slot1
) {
  d <- so
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 3) { #nolint
    d2 <- data.frame(
      d@meta.data,
      `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
      `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2],
      `UMAP.3` = d@reductions[[slot1]]@cell.embeddings[, 3]
    )
    d2_list <- setNames(
      lapply(
        c(md_list),
        function(x) {
          d2[, c(
            x,
            "UMAP.1",
            "UMAP.2",
            "UMAP.3"
          )
          ]
        }
      ),
      c(md_list)
    )

  }
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 2) { #nolint
    d2 <- data.frame(
      d@meta.data,
      `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
      `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2]
    )
    d2_list <- setNames(
      lapply(
        c(md_list),
        function(x) {
          d2[, c(
            x,
            "UMAP.1",
            "UMAP.2"
          )
          ]
        }
      ),
      c(md_list)
    )
  }
  # Generate plots
  d2_plot <- lapply(
    c(md_list),
    function(x) {
      p <-  ggplot2::ggplot(
        d2_list[[x]],
        ggplot2::aes(
          x=`UMAP.1`, # nolint
          y=`UMAP.2`, # nolint
          color = .data[[x]], # nolint
          label = .data[[x]] # nolint
        )
      ) +
        ggplot2::scale_color_manual(
          paste(""),
          values = col_univ()
        ) +
        # Add points
        ggplot2::geom_point(
          shape = 16,
          size = 1,
          alpha = 0.6
        ) +
        sc_theme1() +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          plot.margin = ggplot2::unit(
            c(0.1, 0.1, 0.1, 0.1),
            "cm"
          ),
          legend.position = c(
            0.9,
            0.85
          )
        )
      return(p)
    }
  )
  # Combine output
  if(length(d2_plot) > 1) { # nolint
    d2_out <- ggpubr::ggarrange(
      plotlist = d2_plot,
      labels = c(
        names(
          d2_plot
        )
      ),
      ncol = ifelse(
        length(
          d2_plot
        ) <= 3,
        length(
          d2_plot
        ),
        3
      ),
      nrow = ifelse(
        length(
          d2_plot
        ) > 3,
        ceiling(
          length(
            d2_plot
          ) /
            3
        ),
        1
      )
    )
  }
  if(length(d2_plot) == 1) { # nolint
    d2_out <- d2_plot[[1]]
  }
  return(d2_out)
}

#' Visualize Individual Gene Expression
#'
#' Generates a series of plots to visualize
#' individual gene expression per cluster.
#'
#' @param so An object of class Seurat.
#' @param md_var A character string indicating
#' the clustering column for overlaying on a UMAP plot.
#' @param g_name A character string indicating the gene name.
#' @param col_scheme The color scheme to be used for
#' distinguishing between groups, provided as a vector.
#' @param col_names A vector of the same length as the
#' provided color scheme for assigning colors to each group.
#' @param leg_x A numeric value indicating the placement
#' of the figure legend on the x-axis.
#' @param leg_y A numeric value indicating the placement
#' of the figure legend on the y-axis.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A series of plots stored as a ggplot2 object
#' for visualizing cluster gene expression.
#' @examples
#'
#' # p_umap <- sc_umap_panel_gene(
#' #  d_integrated,
#' #  c("col1","col2","col3"),
#' #  "CFTR",
#' #  col_univ,
#' #  c("group1","group2"),
#' #  0.95,
#' #  0.95,
#' #  "wnn.umap"
#' # )
#'
#' @export
sc_umap_panel_gene <- function(
  so,
  md_var,
  g_name,
  col_scheme,
  col_names,
  leg_x,
  leg_y,
  slot1
) {
  # Format input data
  cols <- setNames(col_scheme,
                   col_names)
  d <- so
  d2 <- data.frame(
    Seurat::FetchData(
      d,
      vars = c(
        gsub(
          "-",
          "-",
          g_name
        ),
        md_var
      )
    ),
    `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
    `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2]
  )

  # Generate plots
  d2_plot <- ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x=`UMAP.1`, # nolint
      y=`UMAP.2`, # nolint
      color = .data[[gsub( # nolint
        "-",
        ".",
        g_name
      )]]
    )
  ) +
    ggplot2::scale_color_gradientn(
      name = "Relative Expression",
      colors = col_grad()
    ) +
    # Add 3D points, axes, and axis-labels
    ggplot2::geom_point(
      shape = 16,
      size = 1,
      alpha = 0.6
    ) +
    ggplot2::ggtitle(g_name) +
    # Add general multivariate plot theme and adjust axis text
    sc_theme1() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(
        c(
          0.1, 0.1, 0.1, 0.1
        ),
        "cm"
      ),
      legend.position = c(
        leg_x,
        leg_y
      )
    )
  ## Metadata overlay
  p_md <-  ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x=`UMAP.1`, # nolint
      y=`UMAP.2`, # nolint
      color = .data[[md_var]], # nolint
      label = .data[[md_var]] # nolint
    )
  ) +
    ggplot2::scale_color_manual(
      paste(""),
      values = cols
    ) +
    ggplot2::geom_point(
      shape = 16,
      size = 1,
      alpha = 0.6
    ) +
    ggrepel::geom_text_repel(data = setNames(
      aggregate(
        d2[, c(
          "UMAP.1",
          "UMAP.2"
        )],
        list(
          d2[[md_var]]
        ),
        FUN = median
      ),
      c(
        md_var,
        names(
          d2[, c(
            "UMAP.1",
            "UMAP.2"
          )]
        )
      )
    ),
    size = 4,
    bg.color = "white") +
    ggplot2::ggtitle(md_var) +
    # Add general multivariate plot theme and adjust axis text
    sc_theme1() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(
        c(
          0.1, 0.1, 0.1, 0.1
        ),
        "cm"
      ),
      legend.position = "none"
    )
  ## Violin Plot
  plot_v <- ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x = .data[[md_var]], # nolint
      y = .data[[gsub(
        "-",
        ".",
        g_name
      )]],
      fill = .data[[md_var]]
    )
  ) +
    ggplot2::scale_fill_manual(
      name = md_var,
      values = col_scheme
    ) +
    # Add violin plot and dotplot
    ggplot2::geom_violin(
      trim = TRUE
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
      show.legend = FALSE
    ) +
    # Add Theme
    sc_theme1() +
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
  d2_out <- ggpubr::ggarrange(
    d2_plot,
    p_md,
    plot_v,
    ncol = 3,
    nrow = 1,
    common.legend = FALSE
  )
  return(d2_out)
}

#' Visualize Gene List Expression
#'
#' Generates a series of plots to visualize
#' the expression of a list of genes per cluster.
#'
#' @param list_g A vector containing a list of genes
#' to be plotted for a given Seurat object.
#' @param so An object of class Seurat.
#' @param md_var A character string indicating the clustering
#' column for overlaying on a UMAP plot.
#' @param col_scheme The color scheme to be used for distinguishing
#' between groups, provided as a vector.
#' @param col_names A vector of the same length as the provided
#' color scheme for assigning colors to each group.
#' @param leg_x A numeric value indicating the placement
#' of the figure legend on the x-axis.
#' @param leg_y A numeric value indicating the placement
#' of the figure legend on the y-axis.
#' @param parl Logical indicating whether processing should
#' be run in parallel (Linux and WSL2 only).
#' Set to FALSE if running sequentially.
#' @param core_perc Percentage of available cores to use if
#' running in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A list of plots saved as ggplot2 objects for
#' visualizing cluster gene expression.
#' @examples
#'
#' # sc_umap_panel_gene_list(
#' #  list_genes,
#' #  d_seurat,
#' #  "seurat_clusters",
#' #  col_vec,
#' #  col_vec_names,
#' #  0.95,
#' #  0.95,
#' #  TRUE,
#' #  0.5,
#' #  "umap"
#' # )
#'
#' @export
sc_umap_panel_gene_list <- function(
  list_g,
  so,
  md_var,
  col_scheme,
  col_names,
  leg_x,
  leg_y,
  parl,
  core_perc,
  slot1
) {
  lg <- list_g
  d <- so
  lg <- unique(lg[lg %in% SeuratObject::Features(d)])
  lg_abs <- subset(lg, !(lg %in% SeuratObject::Features(d)))
  # Create plots
  if(Sys.info()[["sysname"]] != "Windows" && # nolint
      parl == TRUE
  ) {
    parallel::mclapply(
      mc.cores = ceiling(
        parallel::detectCores() *
          core_perc
      ),
      lg,
      function(x) {
        pg <- sc_umap_panel_gene(
          d,
          md_var,
          x,
          col_scheme,
          col_names,
          leg_x,
          leg_y,
          slot1
        )
        # Save each plot
        ggplot2::ggsave(
          paste(
            "analysis/gene/plot.umap.exp.",
            x,
            ".png",
            sep = ""
          ),
          pg,
          height = 12,
          width = 36,
          dpi = 700
        )
      }
    )

  }
  if(Sys.info()[["sysname"]] == "Windows") { # nolint
    lapply(
      lg,
      function(x) {
        pg <- sc_umap_panel_gene(
          d,
          md_var,
          x,
          col_scheme,
          col_names,
          leg_x,
          leg_y,
          slot1
        )
        # Save each plot
        ggplot2::ggsave(
          paste(
            "analysis/gene/plot.umap.exp.",
            x,
            ".png",
            sep = ""
          ),
          pg,
          height = 12,
          width = 36,
          dpi = 700
        )
      }
    )
  }
  if(Sys.info()[["sysname"]] != "Windows" && parl == FALSE) { # nolint
    lapply(
      lg,
      function(x) {
        pg <- sc_umap_panel_gene(
          d,
          md_var,
          x,
          col_scheme,
          col_names,
          leg_x,
          leg_y,
          slot1
        )
        # Save each plot
        ggplot2::ggsave(
          paste(
            "analysis/gene/plot.umap.exp.",
            x,
            ".png",
            sep = ""
          ),
          pg,
          height = 12,
          width = 36,
          dpi = 700
        )
      }
    )

  }
  if(length(lg_abs) > 0) { # nolint
    print(
      paste(
        lg_abs,
        "was not found; plots for this gene will be 
        excluded from the final list...",
        sep = " "
      )
    )
  }
}

#' Standard UMAP Plot
#'
#' Generates a single UMAP plot with a specific
#' metadata overlay for a Seurat object.
#'
#' @param so An object of class Seurat.
#' @param md_var A character string indicating the clustering
#' column for overlaying on a UMAP plot.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A UMAP plot with points grouped by a specific metadata column.
#' @examples
#'
#' # p_umap <- sc_umap_standard(d_integrated,"col1","wnn.umap")
#'
#' @export
sc_umap_standard <- function(
  so,
  md_var,
  slot1
) {
  # Format input data
  d <- so
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 3) { # nolint
    d2 <- data.frame(
      d@meta.data,
      `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
      `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2],
      `UMAP.3` = d@reductions[[slot1]]@cell.embeddings[, 3],
      md.var = d@meta.data[[md_var]]
    )
  }
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 2) { # nolint
    d2 <- data.frame(
      d@meta.data,
      `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
      `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2],
      md.var = d@meta.data[[md_var]]
    )
  }

  # Generate plot
  d2_plot <- ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x=`UMAP.1`, # nolint
      y=`UMAP.2`, # nolint
      color = .data[[md_var]], # nolint
      label = .data[[md_var]] # nolint
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
          d2[, c("UMAP.1", "UMAP.2")],
          list(d2[[md_var]]),
          FUN = median
        ),
        c(md_var, names(
          d2[, c("UMAP.1", "UMAP.2")]
        )
        )
      ),
      size = 4,
      bg.color = "white"
    ) +
    ggplot2::scale_color_manual(
      paste(""),
      values = col_univ()
    ) +
    sc_theme1() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(
        c(0.1, 0.1, 0.1, 0.1), "cm"
      ),
      legend.position = "none"
    )
  return(d2_plot)
}

#' scATAC-Seq Activity UMAP
#'
#' Generates a series of plots given a Seurat scATAC-Seq assay.
#'
#' @param so A Seurat object.
#' @param md_var A character string indicating the clustering
#' column for overlaying on a UMAP plot.
#' @param da_df Differential activity results data frame.
#' @param tf_name Transcription factor name, provided as a character string.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A UMAP plot with points indicating the activity of a specific
#' transcription factor motif.
#' @examples
#'
#' # p_umap_act <- sc_umap_panel_act(
#' #   d,
#' #   "CellType",
#' #   diff.activity.output,
#' #   "GHRL1",
#' #   col_univ()[1:length(levels(d@meta.data[[md_var]]))],
#' #   c(levels(d@meta.data[[md_var]]))
#' #   0.95,
#' #   0.95,
#' #   "wnn.umap.cor"
#' # )
#'
#' @export
sc_umap_panel_act <- function(
  so,
  md_var,
  da_df,
  tf_name,
  col_scheme,
  col_names,
  leg_x,
  leg_y,
  slot1
) {
  d <- so
  # Format input data
  cols <- setNames(col_scheme, col_names)

  d2 <- data.frame(
    Seurat::FetchData(
      d,
      vars = c(
        da_df[da_df[["gene"]] == tf_name, "ID"][[1]],
        md_var
      )
    ),
    `UMAP.1` = d@reductions[[slot1]]@cell.embeddings[, 1],
    `UMAP.2` = d@reductions[[slot1]]@cell.embeddings[, 2]
  )

  # Generate plots
  d2_plot <- ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x=`UMAP.1`, # nolint
      y=`UMAP.2`, # nolint
      color = .data[[da_df[da_df[["gene"]] == tf_name, "ID"][[1]]]] # nolint
    )
  ) +
    ggplot2::scale_color_gradientn(
      name = "Motif Activity",
      colors = col_grad()
    ) +
    # Add 3D points, axes, and axis-labels
    ggplot2::geom_point(
      shape = 16,
      size = 1,
      alpha = 0.6
    ) +
    ggplot2::ggtitle(tf_name) +
    # Add general multivariate plot theme and adjust axis text
    sc_theme1() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(
        c(
          0.1, 0.1, 0.1, 0.1
        ),
        "cm"
      ),
      legend.position = c(
        leg_x,
        leg_y
      )
    )
  ## Metadata overlay
  p_md <-  ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x=`UMAP.1`, # nolint
      y=`UMAP.2`, # nolint
      color = .data[[md_var]], # nolint
      label = .data[[md_var]] # nolint
    )
  ) +
    ggplot2::scale_color_manual(
      paste(""),
      values = cols
    ) +
    ggplot2::geom_point(
      shape = 16,
      size = 1,
      alpha = 0.6
    ) +

    ggrepel::geom_text_repel(data = setNames(
      aggregate(
        d2[, c(
          "UMAP.1",
          "UMAP.2"
        )],
        list(
          d2[[md_var]]
        ),
        FUN = median
      ),
      c(
        md_var,
        names(
          d2[, c(
            "UMAP.1",
            "UMAP.2"
          )]
        )
      )
    ),
    size = 4,
    bg.color = "white") +
    ggplot2::ggtitle(md_var) +
    # Add general multivariate plot theme and adjust axis text
    sc_theme1() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(
        c(
          0.1, 0.1, 0.1, 0.1
        ),
        "cm"
      ),
      legend.position = "none"
    )
  ## Violin Plot
  plot_v <- ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x = .data[[md_var]], # nolint
      y = .data[[da_df[da_df[["gene"]] == tf_name, "ID"][[1]]]],
      fill = .data[[md_var]]
    )
  ) +
    ggplot2::scale_fill_manual(
      name = md_var,
      values = col_univ()[1:length( # nolint
        levels(d@meta.data[[md_var]])
      )
      ]
    ) +
    # Add violin plot and dotplot
    ggplot2::geom_violin(
      trim = TRUE
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
      show.legend = FALSE
    ) +
    # Add Theme
    sc_theme1() +
    ggplot2::labs(
      y = "Motif Activity"
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
  d2_out <- ggpubr::ggarrange(
    d2_plot,
    p_md,
    plot_v,
    ncol = 3,
    nrow = 1,
    common.legend = FALSE
  )
  return(d2_out)
}