#' Install Libraries
#'
#' Installs CRAN and Bioconductor packages (if not already installed) necessary to run package
#'
#' @return Summary of installed packages
#' @examples
#'
#'  # sc.check.lib()
#'
#' @export
sc.check.lib <- function(){
  list.pkg.CRAN <- c(
    "ggplot2","dplyr","ggsci","ggrepel","gtools",
    "viridis","parallel","reshape2","ggpubr","Seurat",
    "SeuratObject","future","circlize","magrittr","lazyeval",
    "shadowtext","SeuratData"
    )
  
  list.pkg.BIOC <- c(
    "SoupX","scDblFinder","SingleCellExperiment","SummarizedExperiment","MAST",
    "BiocGenerics","parallel","ComplexHeatmap","biomaRt","topGO",
    "org.Hs.eg.db","EnhancedVolcano"
    )
  
  # mis.CRAN <- list.pkg.CRAN[!(list.pkg.CRAN %in%  as.character(installed.packages()[order(as.character(installed.packages()[,"Package"])),"Package"]))]
  
  if(length(list.pkg.CRAN) > 0){
    print("Attempting to install/update the following packages...")
    print(paste("CRAN:",paste(list.pkg.CRAN,collapse = ", "),sep = " "))
    lapply(
      list.pkg.CRAN,
      function(x) {
        tryCatch(
          {install.packages(x)},
          error = function(e) {paste("Latest version of",x,"already installed... skipping to next package.",sep = " ")})})
    }

  if("BiocManager" %in% installed.packages()[,"Package"] == FALSE){
    install.packages("BiocManager")
    }
  
  if(length(list.pkg.BIOC) > 0){
    print("Attempting to install/update the following packages...")
    print(paste("Bioconductor:",paste(list.pkg.BIOC,collapse = ", "),sep = " "))
    BiocManager::install(list.pkg.BIOC)
    }

  }



