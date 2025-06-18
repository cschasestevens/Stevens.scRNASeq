#' Install Libraries
#'
#' Installs CRAN and Bioconductor packages (if not already installed)
#'
#' @return Summary of installed packages
#' @importFrom BiocManager install
#' @examples
#'
#'  # sc_check_lib()
#'
#' @export
sc_check_lib <- function() {
  list_pkg_cran <- c(
    "ggplot2",
    "viridis",
    "ggsci",
    "Seurat",
    "SeuratObject",
    "harmony",
    "future",
    "ggpubr",
    "plotly",
    "htmlwidgets",
    "ggrepel",
    "dplyr",
    "reshape2",
    "lazyeval",
    "magrittr",
    "circlize",
    "grid",
    "gtools",
    "patchwork",
    "shadowtext",
    "stringr"
  )

  list_pkg_bioc <- c(
    "SoupX", "scDblFinder", "SingleCellExperiment",
    "SummarizedExperiment", "MAST",
    "BiocGenerics", "parallel", "ComplexHeatmap",
    "biomaRt", "topGO",
    "org.Hs.eg.db", "EnhancedVolcano", "TFBSTools",
    "JASPAR2020", "chromVAR", "BSgenome.Hsapiens.UCSC.hg38",
    "decontX", "BSgenome", "Rsamtools", "GenomeInfoDb",
    "rtracklayer"
  )

  if(length(list_pkg_cran) > 0) { # nolint
    print("Attempting to install/update the following packages...")
    print(paste("CRAN:", paste(list_pkg_cran, collapse = ", "), sep = " "))
    lapply(
      list_pkg_cran,
      function(x) {
        tryCatch(
          {install.packages(x)},
          error = function(e) {
            paste(
              "Latest version of",
              x,
              "already installed... skipping to next package.",
              sep = " "
            )
          }
        )
      }
    )
  }

  if("BiocManager" %in% installed.packages()[, "Package"] == FALSE) { # nolint
    install.packages("BiocManager")
  }

  if(length(list_pkg_bioc) > 0) { # nolint
    print("Attempting to install/update the following packages...")
    print(
      paste(
        "Bioconductor:",
        paste(list_pkg_bioc, collapse = ", "),
        sep = " "
      )
    )
    BiocManager::install(list_pkg_bioc)
  }

}