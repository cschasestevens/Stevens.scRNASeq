#' Top-10 Marker Gene Heatmap (Reclustered)
#'
#' Generates a heatmap from a reclustered Seurat Object
#' and marker gene list based on the top-10 marker genes
#' for each cluster.
#'
#' @param so An object of class Seurat.
#' @param cl_var Character string containing the name of
#' the cluster variable for cell type predictions.
#' @param h_w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param h_h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs_c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs_r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @return A ComplexHeatmap object containing a reclustered
#' top-10 marker gene heatmap.
#' @examples
#'
#' # p_umap <- sc_top10_marker_heatmap(d_annotated,"seurat_clusters",18,24,6,8)
#'
#' @export
sc_top10_marker_heatmap_rc <- function(
  so,
  cl_var,
  h_w,
  h_h,
  fs_c,
  fs_r
) {
  d <- so
  if(!file.exists("analysis/recluster/table.marker.genes.txt")) { # nolint
    print(
      "No marker gene file has been created;
      calculating marker genes for each cluster..."
    )
    cl_mark <- Seurat::FindAllMarkers(d, verbose = TRUE)
    write.table(
      cl_mark,
      "analysis/recluster/table.marker.genes.txt",
      col.names = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
  }
  ## Marker gene input matrix (top10 per cell type)
  d_mark <- read.table(
    "analysis/recluster/table.marker.genes.txt",
    sep = "\t",
    header = TRUE
  )
  d_mark[["CellType.no"]] <- d_mark[["cluster"]]
  ### Top 10 genes per cluster
  d_mark <- dplyr::slice_max(
    dplyr::group_by(
      d_mark,
      .data[["CellType.no"]] # nolint
    ),
    order_by = .data[["avg_log2FC"]],
    n = 10
  )[, c(
    "gene",
    "cluster"
  )]
  #### Save table
  write.table(
    d_mark,
    "analysis/recluster/table.marker.genes.top10.txt",
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
  ### Subset seurat and scale
  SeuratObject::DefaultAssay(d) <- "RNA"
  h <- SeuratObject::FetchData(
    d,
    vars = c(
      cl_var,
      unique(
        d_mark[["gene"]]
      )
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
  ### Scale and plot average expression per cell type
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
    colors = col_grad()[c( # nolint
      1, 3,
      6, 12
    )]
  )
  # Create Plot
  h_out <- ComplexHeatmap::Heatmap(
    h_in,
    col = fun_hm_col,
    name = "Scaled Expression",
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      `Average.Expression` = ComplexHeatmap::anno_barplot(
        as.vector(t(h_anno)),
        gp = grid::gpar(fill = col_grad()) # nolint
      ),
      annotation_name_gp = grid::gpar(
        fontsize = 10
      )
    ),
    show_column_names = TRUE,
    show_row_names = TRUE,
    heatmap_width = ggplot2::unit(h_w, "cm"),
    heatmap_height = ggplot2::unit(h_h, "cm"),
    column_title = "Marker Genes (Top 10)",
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = fs_c),
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = fs_r),
    cluster_columns = FALSE,
    cluster_rows = FALSE
  )
  return(h_out)
}

#' Top-10 Marker Gene Heatmap
#'
#' Generates a heatmap from a Seurat Object and
#' marker gene list based on the top-10 marker genes for each cluster.
#'
#' @param so An object of class Seurat.
#' @param asy Assay to use (ex. "RNA").
#' @param cl_var Character string containing the name
#' of the cluster variable for cell type predictions.
#' @param h_w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param h_h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs_c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs_r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @param cl_c Cluster columns?
#' @param cl_r Cluster rows?
#' @return A ComplexHeatmap object containing a top-10 marker gene heatmap.
#' @examples
#'
#' # p_umap <- sc_top10_marker_heatmap(
#' #   d_annotated,
#' #   "RNA",
#' #   "seurat_clusters",
#' #   18,
#' #   24,
#' #   6,
#' #   8,
#' #   TRUE,
#' #   TRUE
#' # )
#'
#' @export
sc_top10_marker_heatmap <- function(
  so,
  asy,
  cl_var,
  h_w,
  h_h,
  fs_c,
  fs_r,
  cl_c,
  cl_r
) {
  d <- so
  if(!file.exists("analysis/table.marker.genes.txt") && asy == "RNA") { # nolint
    print(
      "No marker gene file has been created;
      calculating marker genes for each cluster..."
    )
    Seurat::DefaultAssay(d) <- "RNA"
    cl_mark <- Seurat::FindAllMarkers(d, verbose = TRUE)
    write.table(
      cl_mark,
      "analysis/table.marker.genes.txt",
      col.names = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
  }

  if(!file.exists("analysis/table.marker.motifs.txt") && asy == "ATAC") { # nolint
    print(
      "No marker motif file has been created;
      calculating marker motifs for each cluster..."
    )
    Seurat::DefaultAssay(d) <- "ATAC"
    cl_mark <- Seurat::FindAllMarkers(
      d,
      min.pct = 0.05,
      verbose = TRUE
    )
    names_motif <- data.frame(
      "gene" = seq.int(1, nrow(d@assays$ATAC@meta.features), 1),
      "near.gene" = paste(
        d@assays$ATAC@meta.features[["nearestGene"]],
        seq.int(1, nrow(d@assays$ATAC@meta.features), 1),
        sep = "."
      ),
      "motif" = paste(
        d@assays$ATAC@meta.features[["seqnames"]],
        paste(
          d@assays$ATAC@meta.features[["start"]],
          d@assays$ATAC@meta.features[["end"]],
          sep = "-"
        ),
        sep = ":"
      )
    )
    cl_mark <- dplyr::left_join(
      cl_mark,
      names_motif,
      by = "gene"
    )
    write.table(
      cl_mark,
      "analysis/table.marker.motifs.txt",
      col.names = TRUE,
      row.names = FALSE,
      sep = "\t"
    )
  }

  if(asy == "RNA") { # nolint
    cl_mark <- read.table(
      "analysis/table.marker.genes.txt",
      sep = "\t",
      header = TRUE
    )
  }

  if(asy == "ATAC") { # nolint
    cl_mark <- read.table(
      "analysis/table.marker.motifs.txt",
      sep = "\t",
      header = TRUE
    )
  }
  ## Marker gene input matrix (top10 per cell type)
  if(class(cl_mark[["cluster"]]) == "character") { # nolint
    cl_mark <- cl_mark[gtools::mixedorder(cl_mark[["cluster"]]), ]
  }
  cl_mark[["CellType.no"]] <- cl_mark[["cluster"]]

  cl_mark <- dplyr::group_by(
    cl_mark,
    .data[["CellType.no"]] # nolint
  )
  ### Top 10 genes per cluster (by p value then by fold change)
  cl_mark <- dplyr::slice_max(
    cl_mark[cl_mark[["avg_log2FC"]] > 0, ],
    order_by = -.data[["p_val_adj"]], # nolint
    n = 25
  )[, c(
    "gene",
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
    "gene",
    "cluster"
  )]

  #### Save table
  if(asy == "RNA") { # nolint
    write.table(
      cl_mark,
      "analysis/table.marker.genes.top10.txt",
      row.names = FALSE,
      col.names = TRUE,
      sep = "\t"
    )
    SeuratObject::DefaultAssay(d) <- "RNA"
  }
  if(asy == "ATAC") { # nolint
    write.table(
      cl_mark,
      "analysis/table.marker.motifs.top10.txt",
      row.names = FALSE,
      col.names = TRUE,
      sep = "\t"
    )
    SeuratObject::DefaultAssay(d) <- "ATAC"
  }
  ### Subset seurat and scale
  h <- SeuratObject::FetchData(
    d,
    vars = c(
      cl_var,
      unique(cl_mark[["gene"]])
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
  ### Scale and plot average expression/accessibility per cell type
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
    colors = col_grad()[c( # nolint
      1, 3,
      6, 12
    )]
  )
  # Create Plot
  h_out <- ComplexHeatmap::Heatmap(
    h_in,
    col = fun_hm_col,
    name = "Scaled Expression/Accessibility",
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      `Average.Expression` = ComplexHeatmap::anno_barplot(
        as.vector(t(h_anno)),
        gp = grid::gpar(fill = col_grad()) # nolint
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
    column_names_rot = 90,
    column_names_gp = grid::gpar(fontsize = fs_c),
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = fs_r),
    cluster_columns = cl_c,
    cluster_rows = cl_r
  )
  return(h_out)
}

#' Top-10 DEG Heatmap
#'
#' Generates a heatmap from a Seurat Object and marker gene
#' list based on the top-10 marker genes for each cluster.
#'
#' @param l_deg A list of DGEA results returned by sc.DGEA().
#' @param so An object of class Seurat.
#' @param asy Character string providing the name of the assay to use.
#' @param cl_var Character string containing the name
#' of the clustering variable.
#' @param hm_w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param hm_h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs_c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs_r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @param c_name Comparison name for plot title, provided as a character string.
#' @return A ComplexHeatmap object containing a top-10 marker gene heatmap.
#' @examples
#'
#' # p.heatmap <- sc_top10_de_da_heatmap(
#' #   dgea_output,
#' #   d_annotated,
#' #   "RNA",
#' #   "seurat.clusters",
#' #   18,
#' #   24,
#' #   6,
#' #   8,
#' #   "grp2 vs. grp1"
#' # )
#'
#' @export
sc_top10_de_da_heatmap <- function(
  l_deg,
  so,
  asy,
  cl_var,
  hm_w,
  hm_h,
  fs_c,
  fs_r,
  c_name
) {
  d1 <- so
  SeuratObject::DefaultAssay(d1) <- asy
  ld <- l_deg[[1]]

  ## Return specified number of DEGs for selected comparison and cell types
  d <- dplyr::select(
    dplyr::slice_max(
      dplyr::group_by(
        ld,
        dplyr::across(
          dplyr::all_of(
            c(
              "Comparison",
              "CellType"
            )
          )
        )
      ),
      order_by = -.data[["H.pval"]], # nolint
      n = 10
    ),
    c(
      "CellType",
      "Comparison",
      "GENE",
      "H.pval"
    )
  )
  ## Remove Mitochondrial/Ribosomal genes?
  mt_rp_remove <- TRUE
  ifelse(
    mt_rp_remove == TRUE,
    d <- d[!grepl("MT-|RP", d[["GENE"]]), ],
    d
  )
  ## Filter DGEA results based on gene list
  d2 <- dplyr::select(
    ld[
      ld[["Comparison"]] %in%
        d[["Comparison"]] &
        ld[["GENE"]] %in%
          unique(
            d[["GENE"]]
          ),
    ],
    c(
      "CellType",
      "Comparison",
      "GENE",
      "logFC"
    )
  )
  ## Format as matrix
  d3 <- reshape2::dcast(
    d2[, -c(5)],
    lazyeval::lazy_eval(
      paste(
        "CellType",
        "+ Comparison ~ GENE",
        sep = ""
      )
    ),
    value.var = "logFC"
  )
  d3[3:ncol(d3)] <- as.data.frame(
    lapply(
      d3[3:ncol(
        d3
      )],
      function(x) {
        ifelse(
          is.na(x),
          0,
          x
        )
      }
    )
  )
  ## Row annotations
  h_row_anno3 <- d3[["CellType"]]
  ## Input matrix
  d3 <- magrittr::set_rownames(
    as.matrix(d3[3:ncol(d3)]),
    d3[["Comparison"]]
  )
  ### Subset seurat and scale
  h <- SeuratObject::FetchData(
    d1,
    vars = c(
      cl_var,
      unique(
        d2[["GENE"]]
      )
    )
  )
  ### Heatmap annotation (average expression)
  h_anno <- as.data.frame(
    lapply(
      h[, 2:ncol(
        h
      )],
      function(x) mean(x)
    )
  )
  ## Color scheme
  fun_hm_col <- circlize::colorRamp2(
    c(
      min(d3), min(d3) / 2,
      0, max(d3) / 2,
      max(d3)
    ),
    colors = c(
      col_grad()[[1]], col_grad()[[3]], # nolint
      "white", col_grad()[[6]],
      col_grad()[[12]]
    )
  )
  ## Output plot
  h_out <- ComplexHeatmap::Heatmap(
    d3,
    col = fun_hm_col,
    name = "logFC",
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      `Avg.Exp` = ComplexHeatmap::anno_barplot(
        as.vector(t(h_anno)),
        gp = grid::gpar(
          fill = col_grad() # nolint
        )
      ),
      annotation_name_gp = grid::gpar(
        fontsize = 8
      ),
      annotation_name_side = "left"
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      "Cluster" = h_row_anno3,
      col = list(
        "Cluster" = setNames(
          as.vector(
            col_univ()[1:length( # nolint
              unique(h_row_anno3)
            )]
          ),
          c(unique(h_row_anno3))
        )
      ),
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontsize = 10)
    ),
    show_column_names = TRUE,
    show_row_names = FALSE,
    row_names_gp = grid::gpar(
      fontsize = fs_r
    ),
    heatmap_width = ggplot2::unit(
      hm_w, "cm"
    ),
    heatmap_height = ggplot2::unit(
      hm_h, "cm"
    ),
    column_title = paste(
      "Top 10 DEG per Cell Type",
      c_name,
      sep = ": "
    ),
    column_names_rot = 90,
    column_names_gp = grid::gpar(
      fontsize = fs_c
    ),
    cluster_columns = TRUE,
    cluster_rows = TRUE
  )
  return(h_out)
}

#' Top-10 DEG Dot Plot
#'
#' Generates a dot plot from a Seurat Object and DEG list based on
#' the top-10 DEGs for a specific cell type and comparison.
#'
#' @param p_type Should a custom gene list or DGEA results object
#' be used for plotting?
#' Type "cstm.list" for custom gene lists and "deg.list"
#' for using dgea.results objects.
#' @param l_deg A list of DGEA results returned by sc.DGEA()
#' or a vector of selected genes for plotting.
#' @param so An object of class Seurat.
#' @param ct A pattern provided as a character string for
#' matching a specific cell type or types.
#' @param cl_var Character string indicating the name of the
#' cluster/cell type variable.
#' @param split_var Logical indicating whether the cluster variable
#' should be stratified by additional group variables.
#' @param list_var (Optional) A vector of character strings
#' indicating the name(s) of up to two group variables
#' for stratifying plot points.
#' @return An input data frame and corresponding dot plot displaying
#' the expression of the top-10 DEGs for a specific cell type.
#' @examples
#'
#' # p_dotplot <- sc_top10_deg_dotplot(
#' #   # Type "deg.list" or "cstm.list" to toggle between inputs
#' #   "deg.list",
#' #   # Name of a custom gene list or dgea.results object
#' #   dgea_output,
#' #   # Seurat object
#' #   d_annotated,
#' #   # Unique character strings corresponding to cell types
#' #   "3.Se|6.Se",
#' #   # Name of clustering variable
#' #   "CellType",
#' #   TRUE,
#' #   # Vector of up to 2 variables for stratifying clustering variables
#' #   c("Knockout","Airway")
#' # )
#'
#' @export
sc_top10_deg_dotplot <- function( # nolint
  p_type,
  l_deg,
  so,
  ct,
  cl_var,
  split_var,
  list_var
) {
  # Load Seurat and change default assay to RNA
  d <- so
  SeuratObject::DefaultAssay(d) <- "RNA"
  c <- ct

  if(p_type == "deg.list") { # nolint
    ## Return specified number of DEGs for selected comparison and cell types
    ld <- l_deg[["DGEA.results"]]
    top10 <- dplyr::select(
      dplyr::slice_max(
        dplyr::group_by(
          ld,
          dplyr::across(
            dplyr::all_of(
              c(
                "Comparison",
                "CellType"
              )
            )
          )
        ),
        order_by = -.data[["H.pval"]], # nolint
        n = 10
      ),
      c(
        "CellType",
        "Comparison",
        "GENE",
        "H.pval"
      )
    )
    ## filter list for chosen cell types
    top10 <- top10[grepl(c, top10[["CellType"]]), ]
    ## return list of unique genes
    top10_pres <- unique(
      top10[["GENE"]][top10[["GENE"]] %in% SeuratObject::Features(d)]
    )
    top10_abs <- subset(
      top10[["GENE"]],
      !(top10[["GENE"]] %in% SeuratObject::Features(d))
    )
  }

  if(p_type == "cstm.list") { # nolint
    top10 <- as.vector(l_deg[[1]])
    top10_pres <- unique(top10[top10 %in% SeuratObject::Features(d)])
    top10_abs <- subset(top10, !(top10 %in% SeuratObject::Features(d)))
  }
  ## select genes from Seurat object
  if(split_var == TRUE) { # nolint
    d1 <- cbind(
      SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          list_var,
          top10_pres
        )
      )
    )
  }
  if(split_var == FALSE) { # nolint
    d1 <- cbind(
      SeuratObject::FetchData(
        d,
        vars = c(
          cl_var,
          top10_pres
        )
      )
    )
  }
  ## Subset based on cell type
  d1 <- d1[grepl(c, d1[[cl_var]]), ]
  d1[["CellType"]] <- d1[[cl_var]]
  ## Count/ratio table for creating dot plots
  d1_prc <- dplyr::bind_rows(
    setNames(
      lapply(
        top10_pres,
        function(x) {
          # Determine average expression of each gene
          ## for 2 variables:
          if(split_var == TRUE && length(list_var) == 2) { # nolint
            ### average expression
            d_avg <- setNames(
              aggregate(
                d1[[x]],
                by = list(
                  d1[,  c("CellType")],
                  d1[,  c(list_var[[1]])],
                  d1[,  c(list_var[[2]])]
                ),
                function(y) mean(y)
              ),
              c(
                "CellType", c(list_var),
                "avg.exp"
              )
            )
            d_avg[["avg.exp"]] <- round(
              d_avg[["avg.exp"]],
              digits = 2
            )
            ### percent expressed
            d_prc <- dplyr::count(
              d1, .data[["CellType"]], # nolint
              .data[[list_var[[1]]]],
              .data[[list_var[[2]]]],
              .data[[x]] > 0
            )
            d_prc <- setNames(
              dplyr::filter(
                d_prc, d_prc[4] == TRUE
              ),
              c(
                "CellType", c(list_var),
                "pres", "n"
              )
            )
            d_cnt <- setNames(
              dplyr::count(
                d1, .data[["CellType"]], # nolint
                .data[[list_var[[1]]]],
                .data[[list_var[[2]]]]
              ),
              c(
                "CellType", c(list_var),
                "n"
              )
            )
            d_comb <- dplyr::left_join(
              d_cnt,
              d_prc,
              by = c(
                "CellType", c(list_var)
              )
            )
            d_comb[is.na(d_comb)] <- 0
            d_comb <- dplyr::mutate(
              d_comb,
              "perc.exp" = round(
                d_comb[["n.y"]] /
                  d_comb[["n.x"]],
                digits = 2
              )
            )
            ### combined
            d_comb_out <- dplyr::left_join(
              d_avg,
              d_comb[,
                     c(
                       "CellType", c(list_var),
                       "perc.exp"
                     )],
              by = c(
                "CellType",
                c(list_var)
              )
            )
          }
          if(split_var == TRUE && length(list_var) < 2) { # nolint
            ### average expression
            d_avg <- setNames(
              aggregate(
                d1[[x]],
                by = list(
                  d1[, c("CellType")],
                  d1[, c(list_var[[1]])]
                ),
                function(y) mean(y)
              ),
              c(
                "CellType", c(list_var),
                "avg.exp"
              )
            )
            d_avg[["avg.exp"]] <- round(
              d_avg[["avg.exp"]],
              digits = 2
            )
            ### percent expressed
            d_prc <- dplyr::count(
              d1,.data[["CellType"]], # nolint
              .data[[list_var[[1]]]],
              .data[[x]] > 0
            )
            d_prc <- setNames(
              dplyr::filter(
                d_prc, d_prc[3] == TRUE
              ),
              c(
                "CellType", c(list_var),
                "pres", "n"
              )
            )
            d_cnt <- setNames(
              dplyr::count(
                d1, .data[["CellType"]], # nolint
                .data[[list_var[[1]]]]
              ),
              c(
                "CellType", c(list_var),
                "n"
              )
            )
            d_comb <- dplyr::left_join(
              d_cnt,
              d_prc,
              by = c(
                "CellType", c(list_var)
              )
            )
            d_comb[is.na(d_comb)] <- 0
            d_comb <- dplyr::mutate(
              d_comb,
              "perc.exp" = round(
                d_comb[["n.y"]] /
                  d_comb[["n.x"]],
                digits = 2
              )
            )
            ### combined
            d_comb_out <- dplyr::left_join(
              d_avg,
              d_comb[,
                     c(
                       "CellType", c(list_var),
                       "perc.exp"
                     )],
              by = c(
                "CellType",
                c(list_var)
              )
            )
          }
          if(split_var == FALSE) { # nolint
            ### average expression
            d_avg <- setNames(
              aggregate(
                d1[[x]],
                by = list(
                  d1[, c("CellType")]
                ),
                function(y) mean(y)
              ),
              c(
                "CellType",
                "avg.exp"
              )
            )
            d_avg[["avg.exp"]] <- round(
              d_avg[["avg.exp"]],
              digits = 2
            )
            ### percent expressed
            d_prc <- dplyr::count(
              d1,.data[["CellType"]], # nolint
              .data[[x]] > 0
            )
            d_prc <- setNames(
              dplyr::filter(
                d_prc, d_prc[2] == TRUE
              ),
              c(
                "CellType",
                "pres", "n"
              )
            )
            d_cnt <- setNames(
              dplyr::count(
                d1, .data[["CellType"]] # nolint
              ),
              c(
                "CellType",
                "n"
              )
            )
            d_comb <- dplyr::left_join(
              d_cnt,
              d_prc,
              by = c(
                "CellType"
              )
            )
            d_comb[is.na(d_comb)] <- 0
            d_comb <- dplyr::mutate(
              d_comb,
              "perc.exp" = round(
                d_comb[["n.y"]] /
                  d_comb[["n.x"]],
                digits = 2
              )
            )
            ### combined
            d_comb_out <- dplyr::left_join(
              d_avg,
              d_comb[,
                     c(
                       "CellType",
                       "perc.exp"
                     )],
              by = c(
                "CellType"
              )
            )
          }
          return(d_comb_out)
        }
      ),
      c(top10_pres)
    ),
    .id = "GENE"
  )
  ## Add row name labels and convert GENE column to factor
  if(split_var == TRUE && length(list_var) == 2) { # nolint

    d1_prc <- data.frame(
      d1_prc,
      "labs" = factor(
        paste(
          d1_prc[["CellType"]],
          d1_prc[[list_var[[1]]]],
          d1_prc[[list_var[[2]]]],
          sep = " "
        ),
        levels = gtools::mixedsort(
          unique(
            paste(
              d1_prc[["CellType"]],
              d1_prc[[list_var[[1]]]],
              d1_prc[[list_var[[2]]]],
              sep = " "
            )
          )
        )
      )
    )
    d1_prc[["GENE"]] <- factor(
      d1_prc[["GENE"]],
      levels = unique(
        d1_prc[["GENE"]]
      )
    )
  }
  ## Add row name labels and convert GENE column to factor
  if(split_var == TRUE && length(list_var) < 2) { # nolint
    d1_prc <- data.frame(
      d1_prc,
      "labs" = factor(
        paste(
          d1_prc[["CellType"]],
          d1_prc[[list_var[[1]]]],
          sep = " "
        ),
        levels = gtools::mixedsort(
          unique(
            paste(
              d1_prc[["CellType"]],
              d1_prc[[list_var[[1]]]],
              sep = " "
            )
          )
        )
      )
    )
    d1_prc[["GENE"]] <- factor(
      d1_prc[["GENE"]],
      levels = unique(
        d1_prc[["GENE"]]
      )
    )
  }
  if(split_var == FALSE) { # nolint
    d1_prc <- data.frame(
      d1_prc,
      "labs" = d1_prc[["CellType"]]
    )
    d1_prc[["GENE"]] <- factor(
      d1_prc[["GENE"]],
      levels = unique(
        d1_prc[["GENE"]]
      )
    )
  }
  ## Plot
  p_dot <- ggplot2::ggplot(
    d1_prc,
    ggplot2::aes(
      x = .data[["GENE"]], # nolint
      y = .data[["labs"]],
      fill = .data[["avg.exp"]],
      size = .data[["perc.exp"]]
    )
  ) +
    ggplot2::geom_point(
      shape = 21
    ) +
    ggplot2::scale_fill_gradientn(
      colors = col_grad() # nolint
    ) +
    ggplot2::scale_size_area(max_size = 10) +
    sc_theme1() + # nolint
    ggplot2::labs(
      fill = "Average Expression",
      size = "Percent Expressed",
      y = ""
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::unit(
        c(.2, .2,
          .2, .2),
        "cm"
      )
    )

  if(exists("top10.abs") && length(top10_abs) == 1) { # nolint
    print(
      paste(
        top10_abs,
        "was not found in the provided Seurat object;
        plots for this gene was excluded from the dot plot...",
        sep = " "
      )
    )
  }

  if(exists("top10.abs") && length(top10_abs) > 1) { # nolint
    print(
      paste(
        top10_abs,
        "were not found in the provided Seurat object;
        these genes were excluded from the dot plot...",
        sep = " "
      )
    )
  }

  return(
    list(
      "Input" = d1_prc,
      "Plot" = p_dot
    )
  )
}

#' Top Differentially Active Motif Heatmap
#'
#' Generates a heatmap from a Seurat Object and
#' differential activity analysis result.
#'
#' @param so An object of class Seurat.
#' @param mot_m Character string providing the name of the assay to use.
#' @param filt_ct Logical indicating whether the data should
#' be filtered based on cell type
#' @param ct_name Character string containing the cell type name to filter.
#' @param cl_var Character string containing the name
#' of the clustering variable.
#' @param split_var Logical indicating whether the cluster variable
#' should be stratified by additional group variables.
#' @param list_var (Optional) A vector of character strings
#' indicating the name(s) of up to two group variables
#' for stratifying plot points.
#' @param top_n Number of motifs to use (filters top-100 motifs per assay or
#' specified cell type and subsequently filters by scaled activity score).
#' @param h_w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param h_h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs_c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs_r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @return A ComplexHeatmap object containing a top motif heatmap.
#' @examples
#'
#' # tf_heatmap <- sc_top_motif_heatmap(
#' #   # Seurat object
#' #   so = d,
#' #   # Differential activity results
#' #   mot_m = diff_output_activity,
#' #   # Filter based on CellType?
#' #   filt_ct = TRUE,
#' #   # Cell type name
#' #   ct_name = "11.Secretory",
#' #   # Clustering column
#' #   cl_var = "CellType",
#' #   # Split by additional variable(s)?
#' #   split_var = TRUE,
#' #   # Additional variable(s) to split plot (if split_var = TRUE)
#' #   list_var = c("Airway"),
#' #   # Number of motifs to use
#' #   top_n = 50,
#' #   # Heatmap width
#' #   h_w = 36,
#' #   # Heatmap height
#' #   h_h = 12,
#' #   # Column font size
#' #   fs_c = 6,
#' #   # Row font size
#' #   fs_r = 8
#' # )
#'
#' @export
sc_top_motif_heatmap <- function( # nolint
  so,
  mot_m,
  filt_ct,
  ct_name,
  cl_var,
  split_var,
  list_var,
  top_n,
  h_w,
  h_h,
  fs_c,
  fs_r
) {
  d <- so
  Seurat::DefaultAssay(d) <- "chromvar"

  ## Motif input matrix (top motifs per cell type)
  d_mark <- mot_m
  if(class(d_mark[["cluster"]]) == "character") { # nolint
    d_mark <- d_mark[gtools::mixedorder(d_mark[["cluster"]]), ]
  }
  if(filt_ct == TRUE) { # nolint
    d_mark <- dplyr::bind_rows(setNames(
      lapply(
        seq.int(1, length(unique(d_mark[["gene"]])), 1),
        function(x) {
          d1 <- d_mark[
            d_mark[["gene"]] == unique(d_mark[["gene"]])[[x]] &
              d_mark[["avg_diff"]] > 0,
          ]
          d1 <- d1[order(d1[["avg_diff"]], decreasing = TRUE), ]
          d1 <- d1[1, ]
          return(d1)
        }
      ),
      c(unique(d_mark[["gene"]]))
    ))
    d_mark <- d_mark[
      d_mark[["cluster"]] == ct_name &
        !is.na(d_mark[["cluster"]]),
    ]
    d_mark[["ID"]] <- gsub(
      "^([^.]*\\.)|\\..*",
      "\\1",
      d_mark[["ID"]]
    )
    #### Save table
    write.table(
      d_mark,
      paste(
        "analysis/table.diff.active.motifs.",
        ct_name,
        ".txt",
        sep = ""
      ),
      row.names = FALSE,
      col.names = TRUE,
      sep = "\t"
    )
  }

  if(filt_ct == FALSE) { # nolint
    d_mark <- dplyr::group_by(
      d_mark,
      .data[["cluster"]] # nolint
    )

    ## Top motifs per cluster (by p value then by fold change)
    d_mark <- dplyr::slice_max(
      d_mark[d_mark[["avg_diff"]] > 0, ],
      order_by = -.data[["p_val_adj"]], # nolint
      n = 100
    )[, c(
      "gene",
      "ID",
      "cluster",
      "avg_diff",
      "p_val_adj"
    )
    ]

    d_mark <- dplyr::group_by(
      d_mark,
      .data[["cluster"]] # nolint
    )

    d_mark <- dplyr::slice_max(
      d_mark,
      order_by = .data[["avg_diff"]], # nolint
      n = top_n
    )[, c(
      "gene",
      "ID",
      "cluster"
    )
    ]
    d_mark[["ID"]] <- gsub(
      "^([^.]*\\.)|\\..*",
      "\\1",
      d_mark[["ID"]]
    )

    #### Save table
    write.table(
      d_mark,
      paste(
        "analysis/table.marker.motifs.Top",
        paste(as.character(top_n)),
        ".txt",
        sep = ""
      ),
      row.names = FALSE,
      col.names = TRUE,
      sep = "\t"
    )
  }

  ### Subset seurat and scale
  ## select genes from Seurat object
  if(split_var == TRUE) { # nolint
    h <- SeuratObject::FetchData(
      d,
      vars = c(
        cl_var,
        list_var,
        unique(
          d_mark[["ID"]]
        )
      )
    )
  }
  if(split_var == FALSE) { # nolint
    h <- SeuratObject::FetchData(
      d,
      vars = c(
        cl_var,
        unique(
          d_mark[["ID"]]
        )
      )
    )
  }
  if(split_var == TRUE && length(list_var) == 2) { # nolint
    h_in <- scale(
      as.matrix(
        magrittr::set_rownames(
          setNames(
            as.data.frame(
              lapply(
                h[, 4:ncol(
                  h
                )
                ],
                function(x) {
                  dplyr::select(
                    aggregate(
                      x,
                      list(
                        h[, 1],
                        h[, 2],
                        h[, 3]
                      ),
                      FUN = mean
                    ),
                    c(4)
                  )
                }
              )
            ),
            names(h[, 4:ncol(h)])
          ),
          paste(
            aggregate(
              h[, 4],
              list(
                h[, 1],
                h[, 2],
                h[, 3]
              ),
              FUN = mean
            )[[1]],
            aggregate(
              h[, 4],
              list(
                h[, 1],
                h[, 2],
                h[, 3]
              ),
              FUN = mean
            )[[2]],
            aggregate(
              h[, 4],
              list(
                h[, 1],
                h[, 2],
                h[, 3]
              ),
              FUN = mean
            )[[3]],
            sep = "."
          )
        )
      ),
      center = TRUE
    )
  }
  if(split_var == TRUE && length(list_var) < 2) { # nolint
    h_in <- scale(
      as.matrix(
        magrittr::set_rownames(
          setNames(
            as.data.frame(
              lapply(
                h[, 3:ncol(
                  h
                )
                ],
                function(x) {
                  dplyr::select(
                    aggregate(
                      x,
                      list(
                        h[, 1],
                        h[, 2]
                      ),
                      FUN = mean
                    ),
                    c(3)
                  )
                }
              )
            ),
            names(h[, 3:ncol(h)])
          ),
          paste(
            aggregate(
              h[, 3],
              list(
                h[, 1],
                h[, 2]
              ),
              FUN = mean
            )[[1]],
            aggregate(
              h[, 3],
              list(
                h[, 1],
                h[, 2]
              ),
              FUN = mean
            )[[2]],
            sep = "."
          )
        )
      ),
      center = TRUE
    )
  }
  if(split_var == FALSE) { # nolint
    h_in <- scale(
      as.matrix(
        magrittr::set_rownames(
          setNames(
            as.data.frame(
              lapply(
                h[, 2:ncol(
                  h
                )
                ],
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
  }


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
        function(x) !anyNA(as.data.frame(h_in)[x])
      )
    )
    ]
  )

  if(split_var == TRUE && length(list_var) == 2) { # nolint
    h_in <- h_in[
      gtools::mixedsort(
        unique(
          paste(
            h[["CellType"]],
            h[[list_var[[1]]]], # nolint
            h[[list_var[[2]]]],
            sep = "."
          )
        )
      ),
    ]
  }
  if(split_var == TRUE && length(list_var) < 2) { # nolint
    h_in <- h_in[
      gtools::mixedsort(
        unique(
          paste(
            h[["CellType"]],
            h[[list_var[[1]]]], # nolint
            sep = "."
          )
        )
      ),
    ]
  }

  tf_names <- data.frame("ID" = colnames(h_in))
  tf_names <- dplyr::left_join(
    tf_names,
    unique(d_mark[, c("ID", "gene")]),
    by = "ID"
  )
  colnames(h_in) <- tf_names[["gene"]]

  fun_hm_col <- circlize::colorRamp2(
    c(
      qs[[1]],
      (qs[[1]]) / 2,
      (qs[[2]]) / 2,
      qs[[2]]
    ),
    colors = col_grad()[c( # nolint
      1, 3,
      6, 12
    )
    ]
  )
  # Create Plot
  if(filt_ct == TRUE) { # nolint
    # Plot
    h_out <- ComplexHeatmap::Heatmap(
      h_in,
      col = fun_hm_col,
      name = "Scaled Activity Score",
      show_column_names = TRUE,
      show_row_names = TRUE,
      heatmap_width = ggplot2::unit(h_w, "cm"),
      heatmap_height = ggplot2::unit(h_h, "cm"),
      column_title = paste("Top Motifs:", ct_name),
      column_names_rot = 90,
      column_names_gp = grid::gpar(fontsize = fs_c),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = fs_r),
      cluster_columns = TRUE,
      cluster_rows = TRUE
    )
  }

  if(filt_ct == FALSE) { # nolint
    h_out <- ComplexHeatmap::Heatmap(
      h_in,
      col = fun_hm_col,
      name = "Scaled Activity Score",
      show_column_names = TRUE,
      show_row_names = TRUE,
      heatmap_width = ggplot2::unit(h_w, "cm"),
      heatmap_height = ggplot2::unit(h_h, "cm"),
      column_title = paste("Top", 10, "Motifs"),
      column_names_rot = 90,
      column_names_gp = grid::gpar(fontsize = fs_c),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = fs_r),
      cluster_columns = TRUE,
      cluster_rows = TRUE
    )
  }
  return(h_out)
}
