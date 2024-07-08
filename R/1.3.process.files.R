#' Process CellRanger Data Files
#'
#' Processes a single sample into a Seurat object for data integration.
#'
#' @param df.par Character string indicating the path to a set of CellRanger files for data processing.
#' @param i A data frame containing sample metadata variables specific to an individual study.
#' @param rho.adj Numeric value indicating a proportion to scale rho values calculated during the ambient RNA removal step. 10% (0.1)
#' is generally suitable for most data sets.
#' @param m.cell Threshold for including features in a Seurat Object. See ?Seurat::CreateSeuratObject() for more details.
#' @param m.feat Threshold for including cells in a Seurat Object. See ?Seurat::CreateSeuratObject() for more details.
#' @return Data frame containing a list of parameters to use for scRNA-Seq data processing by Seurat.
#' @examples
#'
#' # Dataset input parameters
#' list.params <- Create.proc.param(
#'   "data/",
#'   data.frame(
#'     # ID column name (splits by underscore, should be listed first in folder name as in 's01_1_KO_a')
#'     "Code" = unlist(lapply(strsplit(basename(list.files("data/")),"_",fixed = T),"[",1)),
#'     # 1st metadata column (include in CellRanger folder name)
#'     "Knockout" = as.factor(ifelse(grepl("NG",basename(list.files("data/"))),"ctrl","KO")),
#'     # 2nd metadata column
#'     "Region" = as.factor(ifelse(grepl("LAE",basename(list.files("data/"))),"1","2")),
#'     # 3rd metadata column (add/remove columns as needed)
#'     "Time" = as.factor(ifelse(grepl("D28",basename(list.files("data/"))),"a","b"))
#'   )
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
  # Select file from chosen input parameter df
  df.p <- df.par
  d <- df.p[i,]
  # Estimate contamination fraction
  d <- SoupX::autoEstCont(
    SoupX::load10X(
      df.p[i,"Path"]
      )
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
    adj.rho*
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
  d <- assay(
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
  d <- Seurat::AddMetaData(
    object = d,
    metadata = df.p[i,5:ncol(df.p)]
    )
  # Add gene names
  d.g <- setNames(
      read.table(
      df.p[count1,
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
  plot.pre.qc <- VlnPlot(
    object = d,
    features = c("nFeature_RNA",
                 "nCount_RNA",
                 "percent.mt"),
    ncol = 3,
    pt.size = 0.2
    )
  d.filt <- Seurat::subset(
    d,
    subset = nFeature_RNA > 300 &
      nFeature_RNA < 7000 &
      percent.mt < 20
    )
  plot.pos.qc <- VlnPlot(
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

  return(
    list(
      "Seurat.obj" = d.norm,
      "Updated.param.df" = df.p,
      "QC.pre" = plot.pre.qc,
      "QC.post" = plot.pos.qc,
      "var.feat.sum" = d.norm.sum,
      "var.feat.plot" = plot.var.feat.out
      )
    )

  }




