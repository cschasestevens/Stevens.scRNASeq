#' Phenocycler QC
#'
#' Runs quality control of segmented phenocycler data for all samples
#'
#' @param qc_type Specific QC steps to perform (
#' select from "all" or "sum_only")
#' @param ld A list of phenocycler segmentation results containing summary data
#' and segmented cells.
#' @param col_nums1 Vector of column numbers to include from segmentation
#' summary data.
#' @param md_var Vector for renaming summary data column names.
#' @param md_num Number of metadata columns in exported cell segmentation
#' results. Only used if qc_type is "all".
#' @param samp_id Sample ID column name.
#' @param g_col (optional) Add a column providing grouping information.
#' @param g_col2 (optional) Same as g_col, but adds grouping information
#' for individual segmentation results.
#' @param ch_list List of channels exported with segmentation results. Note
#' that the channel list must match the order of channels included within
#' the results.
#' @param ch_nuc Nuclei channel name.
#' @param col_nums2 If qc_type is "all", indicate the column number to
#' use as the sample ID as well as the X and Y coordinates
#' for formatting the segmentation results table.
#' @param loess_norm Normalize signal intensities by performing LOESS
#' prior to z-score scaling.
#' @return List containing QC summary statistics and plots.
#' @examples
#'
#' # pheno_data <- pc_qc(
#' #   ld = ld1,
#' #   col_nums1 = c(4, 10, 12, 13),
#' #   md_var = c("Area.um2", "Cells.raw", "Description"),
#' #   samp_id = "Code",
#' #   g_col = ifelse(grepl("KK", qc1[["Code"]]), "CF", "Norm"),
#' #   g_col2 = ifelse(grepl("KK", qc2[["Code"]]), "CF", "Norm"),
#' #   ch_list = s3,
#' #   col_nums2 = c(6, 8, 9)
#' # )
#'
#' @export
pc_qc <- function(
  qc_type = "all",
  ld,
  col_nums1,
  md_var,
  md_num = 9,
  samp_id = "Code",
  g_col = NULL,
  g_col2 = NULL,
  ch_list,
  col_nums2 = NULL,
  ch_nuc = "DAPI",
  loess_norm = FALSE
) {
  # Cast data frames and convert to individual Seurat objects
  ## Read merged df
  qc2_filt <- readRDS("analysis/qc_int_splitbychannel.rds")
  ## Split data by slide
  d1 <- setNames(
    parallel::mclapply(
      mc.cores = 4,
      seq.int(1, length(unique(qc2_filt[["Slide"]])), 1),
      function(x) {
        d2 <- qc2_filt[qc2_filt[["Slide"]] == x, ]
        return(d2)
      }
    ),
    paste("slide", unique(qc2_filt[["Slide"]]), sep = "")
  )
  lapply(d1, function(x) head(x))
  remove(qc2_filt)
  gc(reset = TRUE)
  d2 <- setNames(
    parallel::mclapply(
      mc.cores = 4,
      seq.int(1, length(d1), 1),
      function(x) {
        d1 <- d1[[x]][, -9]
        d2 <- reshape2::dcast(
          d1,
          ID + Slide + Code + Group + X + Y ~ Channel,
          value.var = "value"
        )
        return(d2)
      }
    ),
    names(d1)
  )
  lapply(d2, function(x) head(x))
  remove(d1)
  gc(reset = TRUE)
  # Convert each data frame into Seurat objects for PCA and UMAP
  d3 <- setNames(
    parallel::mclapply(
      mc.cores = 4,
      seq.int(1, length(d2), 1),
      function(x) {
        d1 <- Seurat::CreateSeuratObject(
          counts = t(as.matrix(d2[[x]][, 7:ncol(d2[[x]])])),
          meta.data = d2[[x]][, 1:6],
          assay = "PC"
        )
        return(d1)
      }
    ),
    names(d2)
  )
  # Convert Seurat objects to individual AnnData files
  parallel::mclapply(
    mc.cores = 4,
    seq.int(1, length(d3), 1),
    function(x) {
      saveRDS(d3[[x]], paste("data/d_slide", x, "_seurat.rds", sep = ""))
    }
  )
  remove(d2)
  gc(reset = TRUE)
  parallel::mclapply(
    mc.cores = 4,
    seq.int(1, length(d3), 1),
    function(x) {
      d2 <- d3[[x]]
      d1 <- Seurat::CreateAssayObject(counts = d2[["PC"]]$counts)
      d2[["PCv3"]] <- d1
      Seurat::DefaultAssay(d2) <- "PCv3"
      d2[["PC"]] <- NULL
      d2 <- SeuratObject::RenameAssays(object = d2, PCv3 = "PC")
      SeuratDisk::SaveH5Seurat(
        d2,
        filename = paste("data/pc_slide", x, ".h5Seurat", sep = "")
      )
      SeuratDisk::Convert(
        paste("data/pc_slide", x, ".h5Seurat", sep = ""),
        dest = "h5ad"
      )
      file.remove(paste("data/pc_slide", x, ".h5Seurat", sep = ""))
      return(d2)
    }
  )









  # Process all Seurat objects and run PCA/UMAP
  d3[["PC"]]
  d3 <- setNames(
    parallel::mclapply(
      mc.cores = 4,
      seq.int(1, length(d3), 1),
      function(x) {
        ## Normalize and run PCA
        options(future.globals.maxSize = 10000 * 1024^2)
        d1 <- Seurat::NormalizeData(d3[[1]])
        d1 <- Seurat::FindVariableFeatures(d1)
        d1 <- Seurat::RunPCA(
          object = Seurat::ScaleData(
            object = d1,
            verbose = TRUE
          ),
          verbose = TRUE
        )
        Seurat::VizDimLoadings(
          object = d1,
          dims = 1:2,
          reduction = "pca"
        )
        Seurat::ElbowPlot(
          d1,
          reduction = "pca",
          ndims = 50
        )
        d1 <- harmony::RunHarmony(
          d1,
          assay.use = "PC",
          group.by.vars = "Code",
          reduction.use = "pca",
          reduction.save = "PC.cor",
          project.dim = FALSE
        )
        ## UMAP
        d1 <- Seurat::RunUMAP(
          d1,
          reduction = "PC.cor",
          reduction.name = "umap.cor",
          reduction.key = "umapcor_",
          dims = 1:30,
          n.components = 3
        )
        ## Find clusters and append to dfs
        d1 <- Seurat::FindNeighbors(d1, dims = 1:30)
        d1 <- Seurat::FindClusters(
          d1, resolution = 0.5, n.iter = 10, random.seed = 1234
        )
        d1 <- Seurat::FindClusters(
          d1, resolution = 0.5, n.iter = 25, random.seed = 1234
        )
        d1 <- Seurat::AddMetaData(
          d1,
          metadata = as.factor(as.numeric(d1@meta.data$seurat_clusters)),
          col.name = "cluster_res0.5_cor"
        )
        return(d1)
      }
    ),
    names(d3)
  )

  # Seurat UMAP not technically feasible for large number of cells (>100K)

  d2 <- data.frame(
    d1@meta.data,
    `UMAP.1` = d1@reductions[["umap.cor"]]@cell.embeddings[, 1],
    `UMAP.2` = d1@reductions[["umap.cor"]]@cell.embeddings[, 2],
    `UMAP.3` = d1@reductions[["umap.cor"]]@cell.embeddings[, 3]
  )
  p2 <-  ggplot2::ggplot(
    d2,
    ggplot2::aes(
      x=`UMAP.1`, # nolint
      y=`UMAP.2`, # nolint
      color = .data[["cluster_res0.5_cor"]], # nolint
      label = .data[["cluster_res0.5_cor"]] # nolint
    )
  ) +
    ggplot2::scale_color_manual(
      paste(""),
      values = col_univ() # nolint
    ) +
    # Add points
    ggplot2::geom_point(
      shape = 16,
      size = 0.5,
      alpha = 0.6
    ) +
    pc_theme_img() + # nolint
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

  ggplot2::ggsave(
    "analysis/test_cluster_seurat_cor.png",
    ggpubr::ggarrange(
      p, p2,
      nrow = 1
    ),
    width = 12,
    height = 6,
    dpi = 300
  )

d1
rownames(d1)  
  p2
  p


  cl_mark <- Seurat::FindAllMarkers(
    d1,
    test.use = "wilcox",
    densify = TRUE,
    verbose = TRUE
  )

  msi_marker_heatmap <- function(
    so,
    h_w,
    h_h,
    fs_c,
    fs_r,
    cl_c,
    cl_r,
    rot_c,
    col1
  ) {
    d <- so
    if(!file.exists("analysis/table.marker.features.txt")) { # nolint
      print(
        "No marker feature file has been created;
        calculating marker features for each cluster..."
      )
      Seurat::DefaultAssay(d) <- "MSI"
      cl_mark <- Seurat::FindAllMarkers(
        d,
        test.use = "wilcox",
        densify = TRUE,
        verbose = TRUE
      )
      write.table(
        cl_mark,
        "analysis/table_markers_slide1_top10.txt",
        row.names = FALSE,
        sep = "\t"
      )
    }
    if(file.exists("analysis/table.marker.features.txt")) { # nolint
      cl_mark <- read.table(
        "analysis/table.marker.features.txt",
        sep = "\t",
        header = TRUE
      )
    }
    ## Marker feature input matrix (top10 per region)
    if(class(cl_mark[["cluster"]]) == "character") { # nolint
      cl_mark <- cl_mark[gtools::mixedorder(cl_mark[["cluster"]]), ]
    }
    cl_mark[["antibody"]] <- cl_mark[["gene"]]
    cl_mark <- dplyr::group_by(
      cl_mark,
      .data[["cluster"]] # nolint
    )
    ### Top 10 features per cluster (by p value then by fold change)
    cl_mark <- dplyr::slice_max(
      cl_mark[cl_mark[["avg_log2FC"]] > 0, ],
      order_by = -.data[["p_val_adj"]], # nolint
      n = 25
    )[, c(
      "antibody",
      "cluster",
      "avg_log2FC",
      "p_val_adj"
    )]
    cl_mark <- dplyr::group_by(
      cl_mark,
      .data[["cluster"]] # nolint
    )
    cl_mark <- dplyr::slice_max(
      cl_mark,
      order_by = .data[["avg_log2FC"]], # nolint
      n = 10
    )[, c(
      "antibody",
      "cluster"
    )]
    ### Subset seurat and scale
    h <- SeuratObject::FetchData(
      d1,
      vars = c(
        "cluster",
        unique(cl_mark[["antibody"]])
      )
    )
    ### Heatmap annotation (average intensity)
    h_anno <- as.data.frame(
      lapply(
        h[, 2:ncol(
          h
        )],
        function(x) {
          mean(x)
        }
      )
    )
  
    h_anno <- h_anno[, h_anno[1, ] > 0]
    ### Scale and plot average intensity per region
    h_in <- scale(
      as.matrix(
        magrittr::set_rownames(
          setNames(
            as.data.frame(
              lapply(
                h[, 2:ncol(
                  h
                )],
                function(x) {
                  dplyr::select(
                    aggregate(
                      x,
                      list(
                        h[, 1]
                      ),
                      FUN = mean
                    ),
                    c(2)
                  )
                }
              )
            ),
            names(h[, 2:ncol(h)])
          ),
          levels(h[, 1])
        )
      ),
      center = TRUE
    )
    qs <- quantile(
      h_in,
      probs = c(
        0.05,
        0.95
      ),
      na.rm = TRUE
    )
    h_in <- as.matrix(
      as.data.frame(h_in)[, unlist(
        lapply(
          seq.int(1, ncol(as.data.frame(h_in)), 1),
          function(x) {
            !anyNA(as.data.frame(h_in)[x])
          }
        )
      )
      ]
    )
    col_umap2 <- c("#2e86c1", "white", "#f5b7b1", "#e74c3c")

    fun_hm_col <- circlize::colorRamp2(
      c(
        qs[[1]],
        (qs[[1]]) / 2,
        (qs[[2]]) / 2,
        qs[[2]]
      ),
      colors = col_umap2
    )
    # Create Plot
    h_out <- ComplexHeatmap::Heatmap(
      h_in,
      col = fun_hm_col,
      name = "Scaled Intensity",
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        `Average.Intensity` = ComplexHeatmap::anno_barplot(
          as.vector(t(h_anno)),
          gp = grid::gpar(fill = col_umap2) # nolint
        ),
        annotation_name_gp = grid::gpar(
          fontsize = 10
        )
      ),
      show_column_names = TRUE,
      show_row_names = TRUE,
      heatmap_width = ggplot2::unit(12, "cm"),
      heatmap_height = ggplot2::unit(12, "cm"),
      column_title = "Top 10 Marker Antibodies",
      column_names_rot = 45,
      column_names_gp = grid::gpar(fontsize = 8),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = 8),
      cluster_columns = FALSE,
      cluster_rows = FALSE
    )
    h_out
    return(list("markers" = cl_mark, "plot" = h_out))
  }





















































  d1
  d2 <- Seurat::CreateAssayObject(counts = d1[["PC"]]$counts)
  d2
  d1[["PCv3"]] <- d2
  d1
  saveRDS(d1, "analysis/pc_data_seurat.rds")
  Seurat::DefaultAssay(d1) <- "PCv3"
  d1
  d1[["PC"]] <- NULL
  d1 <- SeuratObject::RenameAssays(object = d1, PCv3 = "PC")
  SeuratDisk::SaveH5Seurat(d1, filename = "pc_data.h5Seurat")
  SeuratDisk::Convert("pc_data.h5Seurat", dest = "h5ad")
  head(d1@meta.data)

  ## Save as RDS and as AnnData
  return(
    list(
      "summary" = qc1,
      "data" = d2
    )
  )
}

#' Phenocycler UMAP Panel
#'
#' Generates a panel of UMAPs given a
#' Seurat object containing PCA and UMAP results.
#'
#' @param so A Seurat object.
#' @param md_list A vector of character strings indicating
#' metadata columns for overlaying on a loadings plot.
#' @param slot1 A character string corresponding to the umap slot name to plot.
#' @return A series of UMAPs with specified metadata overlays.
#' @examples
#'
#' # pc_umap <- pc_umap_panel(d, c("col1","col2","col3"), "umap.cor")
#'
#' @export
pc_umap_panel <- function(
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
          values = col_univ() # nolint
        ) +
        # Add points
        ggplot2::geom_point(
          shape = 16,
          size = 0.5,
          alpha = 0.6
        ) +
        pc_theme_img() + # nolint
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

#' Top-10 Feature Heatmap
#'
#' Generates a heatmap from a Seurat Object containing
#' clustered MSI data. By default, the top-10 features per
#' cluster are returned.
#'
#' @param so A Seurat object.
#' @param h_w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param h_h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs_c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs_r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @param cl_c Cluster columns?
#' @param cl_r Cluster rows?
#' @param rot_c Rotation of column names.
#' @param col1 Gradient color scheme to use
#' (must be exactly 4 colors in length).
#' @return A ComplexHeatmap object containing a top-10 marker gene heatmap.
#' @examples
#'
#' # p_hmap_top <- msi_marker_heatmap(
#' #   so = d_seg,
#' #   h_w = 18,
#' #   h_h = 24,
#' #   fs_c = 6,
#' #   fs_r = 8,
#' #   cl_c = TRUE,
#' #   cl_r = TRUE,
#' #   rot_c = 45,
#' #   col1 = col_grad()[c(3, 6, 9, 12)]
#' # )
#'
#' @export
msi_marker_heatmap <- function(
  so,
  h_w,
  h_h,
  fs_c,
  fs_r,
  cl_c,
  cl_r,
  rot_c,
  col1
) {
  d <- so
  if(!file.exists("analysis/table.marker.features.txt")) { # nolint
    print(
      "No marker feature file has been created;
      calculating marker features for each cluster..."
    )
    Seurat::DefaultAssay(d) <- "MSI"
    cl_mark <- Seurat::FindAllMarkers(
      d,
      test.use = "wilcox",
      densify = TRUE,
      verbose = TRUE
    )
  }
  if(file.exists("analysis/table.marker.features.txt")) { # nolint
    cl_mark <- read.table(
      "analysis/table.marker.features.txt",
      sep = "\t",
      header = TRUE
    )
  }
  ## Marker feature input matrix (top10 per region)
  if(class(cl_mark[["cluster"]]) == "character") { # nolint
    cl_mark <- cl_mark[gtools::mixedorder(cl_mark[["cluster"]]), ]
  }
  cl_mark[["region"]] <- cl_mark[["cluster"]]
  cl_mark[["feature"]] <- cl_mark[["gene"]]
  cl_mark <- dplyr::group_by(
    cl_mark,
    .data[["region"]] # nolint
  )
  ### Top 10 features per cluster (by p value then by fold change)
  cl_mark <- dplyr::slice_max(
    cl_mark[cl_mark[["avg_log2FC"]] > 0, ],
    order_by = -.data[["p_val_adj"]], # nolint
    n = 25
  )[, c(
    "feature",
    "cluster",
    "avg_log2FC",
    "p_val_adj"
  )]
  cl_mark <- dplyr::group_by(
    cl_mark,
    .data[["cluster"]] # nolint
  )
  cl_mark <- dplyr::slice_max(
    cl_mark,
    order_by = .data[["avg_log2FC"]], # nolint
    n = 10
  )[, c(
    "feature",
    "cluster"
  )]
  ### Subset seurat and scale
  h <- SeuratObject::FetchData(
    d,
    vars = c(
      "cluster",
      unique(cl_mark[["feature"]])
    )
  )
  ### Heatmap annotation (average expression)
  h_anno <- as.data.frame(
    lapply(
      h[, 2:ncol(
        h
      )],
      function(x) {
        mean(x)
      }
    )
  )

  h_anno <- h_anno[, h_anno[1, ] > 0]
  ### Scale and plot average intensity per region
  h_in <- scale(
    as.matrix(
      magrittr::set_rownames(
        setNames(
          as.data.frame(
            lapply(
              h[, 2:ncol(
                h
              )],
              function(x) {
                dplyr::select(
                  aggregate(
                    x,
                    list(
                      h[, 1]
                    ),
                    FUN = mean
                  ),
                  c(2)
                )
              }
            )
          ),
          names(h[, 2:ncol(h)])
        ),
        levels(h[, 1])
      )
    ),
    center = TRUE
  )
  qs <- quantile(
    h_in,
    probs = c(
      0.05,
      0.95
    ),
    na.rm = TRUE
  )
  h_in <- as.matrix(
    as.data.frame(h_in)[, unlist(
      lapply(
        seq.int(1, ncol(as.data.frame(h_in)), 1),
        function(x) {
          !anyNA(as.data.frame(h_in)[x])
        }
      )
    )
    ]
  )
  fun_hm_col <- circlize::colorRamp2(
    c(
      qs[[1]],
      (qs[[1]]) / 2,
      (qs[[2]]) / 2,
      qs[[2]]
    ),
    colors = col1
  )
  # Create Plot
  h_out <- ComplexHeatmap::Heatmap(
    h_in,
    col = fun_hm_col,
    name = "Scaled Intensity",
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      `Average.Intensity` = ComplexHeatmap::anno_barplot(
        as.vector(t(h_anno)),
        gp = grid::gpar(fill = col1) # nolint
      ),
      annotation_name_gp = grid::gpar(
        fontsize = 10
      )
    ),
    show_column_names = TRUE,
    show_row_names = TRUE,
    heatmap_width = ggplot2::unit(h_w, "cm"),
    heatmap_height = ggplot2::unit(h_h, "cm"),
    column_title = "Top 10 Markers",
    column_names_rot = rot_c,
    column_names_gp = grid::gpar(fontsize = fs_c),
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = fs_r),
    cluster_columns = cl_c,
    cluster_rows = cl_r,
  )
  return(list("markers" = cl_mark, "plot" = h_out))
}

sc_umap_standard <- function(
  so,
  md_var,
  slot1,
  dims1 = "2D",
  col1 = col_univ(), # nolint
  pos_leg = "none"
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

  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 3 && dims1 == "3D") { # nolint
    # Generate plot
    d2_plot <- plotly::plot_ly(
      d2,
      x = ~`UMAP.1`, # nolint
      y = ~`UMAP.2`, # nolint
      z = ~`UMAP.3`, # nolint
      color = ~.data[[md_var]], # nolint
      colors = col_univ() # nolint
    ) %>% # nolint
      plotly::add_markers(marker = list(size = 3)) %>%
      plotly::layout(
        autosize = FALSE,
        width = 800,
        height = 600,
        margin = list(
          l = 50,
          r = 50,
          b = 25,
          t = 25,
          pad = 1
        )
      )
    htmlwidgets::saveWidget(
      d2_plot,
      file = paste("analysis/plot.3D.umap.", md_var, ".html", sep = "")
    )
  }

  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 2 && dims1 == "3D") { # nolint
    print(
      "Error: Only 2 components are present
      in the selected dimension reduction!"
    )
  }

  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 2 && dims1 == "2D") { # nolint
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
        size = 5.5,
        bg.color = "grey0",
        color = "grey55",
        bg.r = 0.075
      ) +
      ggplot2::scale_color_manual(
        paste(""),
        values = col1
      ) +
      sc_theme1() + # nolint
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
        legend.position = pos_leg
      )
  }
  if(ncol(d@reductions[[slot1]]@cell.embeddings) == 3 && dims1 == "2D") { # nolint
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
        size = 5.5,
        bg.color = "grey0",
        color = "grey55",
        bg.r = 0.075
      ) +
      ggplot2::scale_color_manual(
        paste(""),
        values = col1
      ) +
      sc_theme1() + # nolint
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
        legend.position = pos_leg
      )
  }
  return(d2_plot)
}

library(dplyr)
sc_umap_standard(
  so = d1,
  md_var = "cluster_res0.5_cor",
  slot1 = "umap.cor",
  dims1 = "3D",
  col1 = col_univ(),
  pos_leg = "none"
)
