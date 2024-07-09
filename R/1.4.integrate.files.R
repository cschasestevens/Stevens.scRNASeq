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
  #### Integration ####

  fun.integration <- function(
    subset1,
    rds.path
  ) {

    ## Anchor function

    f.anchor.fun <- function(
    s1
    ) {

      ## Find integration anchors

      data.files.anchor <- FindIntegrationAnchors(
        object.list = s1,
        anchor.features = 4000,
        dims = 1:50)

      return(data.files.anchor)

    }


    ## run for chosen subset

    data.files.anchor1 <- f.anchor.fun(
      subset1
    )


    ## Integration function

    d.int.fun <- function(
    x
    ) {

      data.files.merged <- IntegrateData(
        anchorset = x,
        dims = 1:50
      )

      return(data.files.merged)

    }


    data.files.merged1 <- d.int.fun(
      data.files.anchor1
    )

    saveRDS(
      data.files.merged1,
      rds.path
    )

    return(data.files.merged1)

  }


  # Load .rds containing Seurat objects to integrate

  sc.orig <- readRDS(
    "Reference.sets/Ken_HBE10xLSAE5codes1mSAE_cls19noLQ15_seurat.rds"
  )

  sc.nkx <- readRDS(
    "Analysis/RDS/1_NKX21KO_data_annotated.rds"
  )

  sc.nkx2 <- subset(
    sc.nkx,
    subset = Time == "D28"
  )

  remove(sc.nkx)

  # Perform Integration

  d.merged <- fun.integration(
    c(
      sc.orig,
      sc.nkx2
    ),
    "Analysis/RDS/1.all.lsae.data.rds"
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





