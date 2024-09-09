#' Process CellRanger Data Files
#'
#' Processes a single sample into a Seurat object for data integration.
#'
#' @param df_par Data frame of parameters to use for data processing.
#' @param i Numeric value for sample ID.
#' @param rho_adj Numeric value indicating a proportion to scale
#' rho values calculated during the ambient RNA removal step. 10% (0.1)
#' is generally suitable for most data sets.
#' @param m_cell Threshold for including features in a Seurat Object.
#' Run ?Seurat::CreateSeuratObject() for more details.
#' @param m_feat Threshold for including cells in a Seurat Object.
#' Run ?Seurat::CreateSeuratObject() for more details.
#' @return A processed sample file converted into a Seurat object
#' with a summary list of QC and processing details.
#' @examples
#'
#' # proc_data <- sc_process_file(
#' # # parameter list
#' # list_params,
#' # # sample ID
#' # 1,
#' # # adj.rho proportion
#' # 0.1,
#' # # minimum cells per feature
#' # 5,
#' # # minimum features per cell
#' # 200
#' # )
#'
#' @export
sc_process_file <- function(
  df_par,
  i,
  rho_adj,
  m_cell,
  m_feat
) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1234)
  # Select file from chosen input parameter df
  df_p <- df_par
  d <- df_p[i, ]
  # Estimate contamination fraction
  d <- SoupX::autoEstCont(
    SoupX::load10X(
      df_p[i, "Path"]
    ),
    doPlot = FALSE
  )
  df_p[["rho"]] <- as.vector(
    unlist(
      unique(
        d[["metaData"]][["rho"]]
      )
    )
  )
  ## Adjusted rho
  df_p[["adj.rho"]] <- df_p[["rho"]] +
    rho_adj *
      df_p[["rho"]]
  ## Cluster and cell numbers
  df_p[["Clusters"]] <- length(
    unique(
      (d)$metaData$clusters
    )
  )
  df_p[["Cell.No"]] <- nrow(d$metaData)
  # Adjust contamination fraction by selected adj.rho value
  d <- SoupX::adjustCounts(
    SoupX::setContaminationFraction(
      d,
      df_p[["adj.rho"]]
    )
  )
  # Doublet removal/Create single cell experiment
  d <- scDblFinder::scDblFinder(
    SingleCellExperiment::SingleCellExperiment(
      list(
        counts = as.matrix(d)
      )
    )
  )
  d <- SummarizedExperiment::assay(
    d[, d$scDblFinder.class == "singlet"],
    "counts"
  )
  # Create Seurat object
  d <- Seurat::CreateSeuratObject(
    counts = d,
    project = paste(df_p[i, "File.ID"]),
    min.cells = m_cell,
    min.features = m_feat
  )
  d_md <- df_p[i, 5:ncol(df_p)]
  d_md <- setNames(
    as.data.frame(
      lapply(
        seq.int(
          1,
          ncol(d_md),
          1
        ),
        function(x) {
          rep(
            d_md[1, x],
            nrow(d@meta.data)
          )
        }
      )
    ),
    names(d_md)
  )
  d <- Seurat::AddMetaData(
    object = d,
    metadata = d_md
  )
  # Add gene names
  d_g <- setNames(
    read.table(
      df_p[i,
           c("Path.feat")],
      sep = "\t",
      header = FALSE
    ),
    c("GENEID", "GENE", "TYPE")
  )
  d_g <- d_g[!duplicated(d_g[["GENE"]]), ]
  rownames(d_g) <- d_g[["GENE"]]
  d[["RNA"]] <- Seurat::AddMetaData(
    object = d[["RNA"]],
    metadata = d_g
  )
  # Subset Seurat object
  d[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    object = d,
    pattern = "^MT-"
  )
  plot_pre_qc <- Seurat::VlnPlot(
    object = d,
    features = c("nFeature_RNA",
                 "nCount_RNA",
                 "percent.mt"),
    ncol = 3,
    pt.size = 0.2
  )
  d_filt <- BiocGenerics::subset(
    d,
    subset = nFeature_RNA > 300 & # nolint
      nFeature_RNA < 7000 &  # nolint
      percent.mt < 20 # nolint
  )
  plot_pos_qc <- Seurat::VlnPlot(
    object = d_filt,
    features = c("nFeature_RNA",
                 "nCount_RNA",
                 "percent.mt"),
    ncol = 3,
    pt.size = 0.2
  )
  # Log normalize data
  d_norm <- Seurat::NormalizeData(
    object = d_filt,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  d_norm <- Seurat::FindVariableFeatures(
    object = d_norm,
    selection.method = "vst",
    nfeatures = 2000
  )
  d_norm_sum <- head(
    Seurat::VariableFeatures(d_norm),
    25
  )
  plot_var_feat <- Seurat::VariableFeaturePlot(
    d_norm,
    assay = "RNA"
  )
  plot_var_feat_out <- Seurat::LabelPoints(
    plot_var_feat,
    points = d_norm_sum,
    repel = TRUE
  )
  df_p <- df_p[i, ]

  return(
    list(
      "Seurat.obj" = d_norm,
      "Params" = df_p,
      "QC.pre" = plot_pre_qc,
      "QC.post" = plot_pos_qc,
      "var.feat.sum" = d_norm_sum,
      "var.feat.plot" = plot_var_feat_out
    )
  )

}





#' Batch Processing of CellRanger Files
#'
#' Processes a list of scRNA-Seq samples for data integration.
#'
#' @param df_par Data frame of parameters to use for data processing.
#' @param rho_adj Numeric value indicating a proportion to scale
#' rho values calculated during the ambient RNA removal step. 10% (0.1)
#' is generally suitable for most data sets.
#' @param m_cell Threshold for including features in a Seurat Object.
#' Run ?Seurat::CreateSeuratObject() for more details.
#' @param m_feat Threshold for including cells in a Seurat Object.
#' Run ?Seurat::CreateSeuratObject() for more details.
#' @param parl Logical indicating whether processing should be run in
#' parallel (Linux and WSL2 only). Set to FALSE if running sequentially.
#' @param core_perc Percentage of available cores to use if running
#' in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @return A processed list of sample files converted into Seurat
#' objects with a summary list of QC and processing details.
#' @examples
#'
#' # list_data <- sc_process_batch(
#' # # parameter list
#' # list_params,
#' # # adj.rho proportion
#' # 0.1,
#' # # minimum cells per feature
#' # 5,
#' # # minimum features per cell
#' # 200
#' # )
#'
#' @export
sc_process_batch <- function(
  df_par,
  rho_adj,
  m_cell,
  m_feat,
  parl,
  core_perc
) {
  if( # nolint
    Sys.info()[["sysname"]] != "Windows" && parl == TRUE
  ) {
    list_data <- setNames(
      parallel::mclapply(
        mc.cores = ceiling(
          parallel::detectCores() *
            core_perc
        ),
        seq.int(
          1,
          nrow(df_par),
          1
        ),
        function(x) {
          sc_process_file(
            # parameter list
            list_params, # nolint
            # sample ID
            x,
            # adj.rho proportion
            rho_adj,
            # minimum cells per feature
            m_cell,
            # minimum features per cell
            m_feat
          )
        }
      ),
      df_par[["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] != "Windows" && parl == FALSE) { # nolint
    list_data <- setNames(lapply(
      seq.int(
        1,
        nrow(df_par),
        1
      ),
      function(x) {
        sc_process_file(
          # parameter list
          list_params, # nolint
          # sample ID
          x,
          # adj.rho proportion
          rho_adj,
          # minimum cells per feature
          m_cell,
          # minimum features per cell
          m_feat
        )
      }
    ),
    df_par[["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] == "Windows" && parl == TRUE) { # nolint
    print(
      "Windows OS does not support parallel sample processing, 
      defaulting to sequential processing..."
    )
    list_data <- setNames(lapply(
      seq.int(
        1,
        nrow(df_par),
        1
      ),
      function(x) {
        sc_process_file(
          # parameter list
          list_params, # nolint
          # sample ID
          x,
          # adj.rho proportion
          rho_adj,
          # minimum cells per feature
          m_cell,
          # minimum features per cell
          m_feat
        )
      }
    ),
    df_par[["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] == "Windows" && parl == FALSE) { # nolint
    list_data <- setNames(lapply(
      seq.int(
        1,
        nrow(df_par),
        1
      ),
      function(x) {
        sc_process_file(
          # parameter list
          list_params, # nolint
          # sample ID
          x,
          # adj.rho proportion
          rho_adj,
          # minimum cells per feature
          m_cell,
          # minimum features per cell
          m_feat
        )
      }
    ),
    df_par[["File.ID"]]
    )
  }
  return(list_data)
}

#' Plot Sample Contamination Fraction
#'
#' Generates a scatter plot indicating the individual
#' and average contamination fractions for all samples.
#'
#' @param df_par Data frame of updated parameters used for data processing.
#' @return A scatter plot indicating the individual
#' and average contamination fractions for all samples.
#' @examples
#'
#' # sc_plot_rho(list_params)
#'
#' @export
sc_plot_rho <- function(
  df_par
) {
  p_rho <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = df_par,
      ggplot2::aes(
        x = .data[["File.ID"]], # nolint
        y = .data[["rho"]]
      ),
      shape = 21,
      size = 3,
      alpha = 0.25,
      fill = col_univ()[2] # nolint
    ) +
    ggplot2::geom_point(
      data = df_par,
      ggplot2::aes(
        x = .data[["File.ID"]],
        y = .data[["adj.rho"]]
      ),
      shape = 21,
      size = 3,
      alpha = 1,
      fill = col_univ()[1]
    ) +
    ggplot2::geom_hline(
      yintercept = mean(
        df_par[["rho"]]
      ),
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = mean(
        df_par[["adj.rho"]]
      ),
      linetype = "dashed",
      color = "firebrick2"
    ) +
    sc_theme1() # nolint
  ggplot2::ggsave(
    "processed/data.ambientRNA.cont.png",
    p_rho,
    width = 10,
    height = 10,
    dpi = 600
  )
}
