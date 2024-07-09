#' Process CellRanger Data Files
#'
#' Processes a single sample into a Seurat object for data integration.
#'
#' @param df.par Data frame of parameters to use for data processing.
#' @param i Numeric value for sample ID.
#' @param rho.adj Numeric value indicating a proportion to scale rho values calculated during the ambient RNA removal step. 10% (0.1)
#' is generally suitable for most data sets.
#' @param m.cell Threshold for including features in a Seurat Object. Run ?Seurat::CreateSeuratObject() for more details.
#' @param m.feat Threshold for including cells in a Seurat Object. Run ?Seurat::CreateSeuratObject() for more details.
#' @return A processed sample file converted into a Seurat object with a summary list of QC and processing details.
#' @examples
#'
#' proc.data <- sc.process.file(
#' # parameter list
#' list.params,
#' # sample ID
#' 1,
#' # adj.rho proportion
#' 0.1,
#' # minimum cells per feature
#' 5,
#' # minimum features per cell
#' 200
#' )
#'
#' @export
sc.process.file <- function(
    df.par,
    i,
    rho.adj,
    m.cell,
    m.feat
    ){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(1234)
  # Select file from chosen input parameter df
  df.p <- df.par
  d <- df.p[i,]
  # Estimate contamination fraction
  d <- SoupX::autoEstCont(
    SoupX::load10X(
      df.p[i,"Path"]
      ),
    doPlot = FALSE
    )
  df.p[["rho"]] <- as.vector(
    unlist(
      unique(
        d[["metaData"]][["rho"]]
        )
      )
    )
  ## Adjusted rho
  df.p[["adj.rho"]] <- df.p[["rho"]] +
    rho.adj*
    df.p[["rho"]]
  ## Cluster and cell numbers
  df.p[["Clusters"]] <- length(
    unique(
      (d)$metaData$clusters
      )
    )
  df.p[["Cell.No"]] <- nrow(d$metaData)
  # Adjust contamination fraction by selected adj.rho value
  d <- SoupX::adjustCounts(
    SoupX::setContaminationFraction(
      d,
      df.p[["adj.rho"]]
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
    d[,d$scDblFinder.class == "singlet"],
    "counts"
    )
  # Create Seurat object
  d <- Seurat::CreateSeuratObject(
    counts = d,
    project = paste(df.p[i,"File.ID"]),
    min.cells = m.cell,
    min.features = m.feat
    )
  d.md <- df.p[i,5:ncol(df.p)]
  d.md <- setNames(
    as.data.frame(
      lapply(
        seq.int(
          1,
          ncol(d.md),
          1
          ),
        function(x)
          rep(
            d.md[1,x],
            nrow(d@meta.data)
            )
        )
      ),
    names(d.md)
    )
  d <- Seurat::AddMetaData(
    object = d,
    metadata = d.md
    )
  # Add gene names
  d.g <- setNames(
      read.table(
      df.p[i,
           c("Path.feat")],
      sep = "\t",
      header = F
      ),
      c("GENEID","GENE","TYPE")
      )
  d.g <- d.g[!duplicated(d.g[["GENE"]]),]
  rownames(d.g) <- d.g[["GENE"]]
  d[["RNA"]] <- Seurat::AddMetaData(
    object = d[["RNA"]],
    metadata = d.g
    )
  # Subset Seurat object
  d[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    object = d,
    pattern = "^MT-"
    )
  plot.pre.qc <- Seurat::VlnPlot(
    object = d,
    features = c("nFeature_RNA",
                 "nCount_RNA",
                 "percent.mt"),
    ncol = 3,
    pt.size = 0.2
    )
  d.filt <- BiocGenerics::subset(
    d,
    subset = nFeature_RNA > 300 &
      nFeature_RNA < 7000 &
      percent.mt < 20
    )
  plot.pos.qc <- Seurat::VlnPlot(
    object = d.filt,
    features = c("nFeature_RNA",
                 "nCount_RNA",
                 "percent.mt"),
    ncol=3,
    pt.size=0.2
    )
  # Log normalize data
  d.norm <- Seurat::NormalizeData(
    object = d.filt,
    normalization.method = "LogNormalize",
    scale.factor = 10000
    )
  d.norm <- Seurat::FindVariableFeatures(
    object = d.norm,
    selection.method = "vst",
    nfeatures = 2000
    )
  d.norm.sum <- head(
    Seurat::VariableFeatures(d.norm),
    25
    )
  plot.var.feat <- Seurat::VariableFeaturePlot(
    d.norm,
    assay = "RNA"
    )
  plot.var.feat.out <- Seurat::LabelPoints(
    plot.var.feat,
    points = d.norm.sum,
    repel = T
    )
  df.p <- df.p[i,]

  return(
    list(
      "Seurat.obj" = d.norm,
      "Params" = df.p,
      "QC.pre" = plot.pre.qc,
      "QC.post" = plot.pos.qc,
      "var.feat.sum" = d.norm.sum,
      "var.feat.plot" = plot.var.feat.out
      )
    )

  }





#' Batch Processing of CellRanger Files
#'
#' Processes a list of samples into separate Seurat objects for data integration.
#'
#' @param df.par Data frame of parameters to use for data processing.
#' @param rho.adj Numeric value indicating a proportion to scale rho values calculated during the ambient RNA removal step. 10% (0.1)
#' is generally suitable for most data sets.
#' @param m.cell Threshold for including features in a Seurat Object. Run ?Seurat::CreateSeuratObject() for more details.
#' @param m.feat Threshold for including cells in a Seurat Object. Run ?Seurat::CreateSeuratObject() for more details.
#' @param parl Logical indicating whether processing should be run in parallel (Linux and WSL2 only). Set to FALSE if running sequentially.
#' @param core.perc Percentage of available cores to use if running in parallel (Linux and WSL2 only). Set to 1 if running sequentially.
#' @return A processed list of sample files converted into Seurat objects with a summary list of QC and processing details.
#' @examples
#'
#' list.data <- sc.process.batch(
#' # parameter list
#' list.params,
#' # adj.rho proportion
#' 0.1,
#' # minimum cells per feature
#' 5,
#' # minimum features per cell
#' 200
#' )
#'
#' @export
sc.process.batch <- function(
    df.par,
    rho.adj,
    m.cell,
    m.feat,
    parl,
    core.perc
    ){
  if(Sys.info()[["sysname"]] != "Windows" &
     parl == TRUE) {
  list.data <- setNames(
    parallel::mclapply(
      mc.cores = ceiling(
        parallel::detectCores()*
          core.perc
        ),
      seq.int(
        1,
        nrow(df.par),
        1
        ),
      function(x) {
        sc.process.file(
          # parameter list
          list.params,
          # sample ID
          x,
          # adj.rho proportion
          rho.adj,
          # minimum cells per feature
          m.cell,
          # minimum features per cell
          m.feat
          )
        }
      ),
    df.par[["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] != "Windows" &
     parl == FALSE) {
    list.data <- setNames(lapply(
      seq.int(
        1,
        nrow(df.par),
        1
      ),
      function(x) {
        sc.process.file(
          # parameter list
          list.params,
          # sample ID
          x,
          # adj.rho proportion
          rho.adj,
          # minimum cells per feature
          m.cell,
          # minimum features per cell
          m.feat
        )
      }
    ),
    df.par[["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] == "Windows" &
     parl == TRUE) {
    print("Windows OS does not support parallel sample processing, defaulting to sequential processing...")
    list.data <- setNames(lapply(
      seq.int(
        1,
        nrow(df.par),
        1
      ),
      function(x) {
        sc.process.file(
          # parameter list
          list.params,
          # sample ID
          x,
          # adj.rho proportion
          rho.adj,
          # minimum cells per feature
          m.cell,
          # minimum features per cell
          m.feat
        )
      }
    ),
    df.par[["File.ID"]]
    )
  }
  if(Sys.info()[["sysname"]] == "Windows" &
     parl == FALSE) {
    list.data <- setNames(lapply(
      seq.int(
        1,
        nrow(df.par),
        1
      ),
      function(x) {
        sc.process.file(
          # parameter list
          list.params,
          # sample ID
          x,
          # adj.rho proportion
          rho.adj,
          # minimum cells per feature
          m.cell,
          # minimum features per cell
          m.feat
        )
      }
    ),
    df.par[["File.ID"]]
    )
  }
  return(list.data)
}











