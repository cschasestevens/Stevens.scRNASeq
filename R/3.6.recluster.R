#' Reclustering Analysis
#'
#' Performs reclustering and basic visualization of a subsetted Seurat Object. Subsets are determined based on CellGroup.
#'
#' @param so An object of class Seurat. Must contain the columns 'CellType' and 'CellGroup' in the metadata slot.
#' @param ct Character string vector containing the name(s) of up to 2 cell type groups present in the CellGroup column.
#' The provided Seurat object is then subsetted to include only clusters from the provided group names.
#' @param g.list Path to an existing text file containing genes for plotting expression.
#' @param md.list A vector of character strings indicating metadata columns for overlaying on a loadings plot.

#' @param h.h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs.c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs.r Numeric value for row fontsize (passed to ComplexHeatmap).

#' @return A ComplexHeatmap object containing a top-10 marker gene heatmap.
#' @examples
#'
#' p.umap <- sc.top10.marker.heatmap(d.annotated,"seurat.clusters",18,24,6,8)
#'
#' @export


fun.marker.gene.recluster <- function(
    so,
    ct,
    g.list,
    md.list
    ) {
  # Load gene list and Seurat object
  list.genes <- read.table(
    g.list,
    sep = "\t",
    header = T)
  d <- so

  # Subset Seurat
  if(length(ct) == 2) {
    d <- BiocGenerics::subset(
      d,
      subset = CellGroup == ct[[1]]|CellGroup == ct[[2]]
      )
    }
  if(length(ct) == 1) {
    d <- BiocGenerics::subset(
      d,
      subset = CellGroup == ct
      )
    }

  # Run PCA/UMAP
  d <- sc.pca(d)
  d.pca <- sc.pca.plot(
    d,
    c(md.list,"seurat_clusters")
    )
  ggplot2::ggsave(
    "analysis/recluster/plot.pca.panel.png",
    d.pca,
    width = 36,
    height = 12,
    dpi = 700)

  # Determine optimal number of dimensions to use, then run UMAP
  d <- SeuratObject::SetIdent(d, value = d@meta.data$CellType)
  d <- Seurat::RunUMAP(
    object = d,
    dims = 1:20,
    n.components = 3
    )
  SeuratObject::DefaultAssay(d) <- "integrated"
  d <- Seurat::FindNeighbors(
    d,
    reduction = "umap",
    dims = 1:3,
    verbose = T
    )
  d <- Seurat::FindClusters(
    d,
    cluster.name = "recluster",
    resolution = 0.3
    )
  d.umap1 <- sc.umap.panel(
    d,
    c("CellType","recluster",md.list)
    )
  ggplot2::ggsave(
    "analysis/recluster/plot.umap.panel.png",
    d.umap1,
    height = 16,
    width = 32,
    dpi = 400
    )

  # Recluster QC
  ggplot2::ggsave(
    "analysis/recluster/plot.recluster.qc.png",
    sc.integration.qc(d,"recluster"),
    width = 24,
    height = 8,
    dpi = 700)

  # Recluster UMAP
  ### Cell Type
  ggplot2::ggsave(
    "analysis/recluster/plot.umap.CellType.png",
    sc.umap.standard(
      # Seurat object
      d,
      # metadata column
      "CellType"
      ),
    height = 8,
    width = 8,
    dpi = 700
    )
  ### Cell Group
  ggplot2::ggsave(
    "analysis/recluster/plot.umap.reclustering.png",
    sc.umap.standard(
      # Seurat object
      d,
      # metadata column
      "recluster"
      ),
    height = 8,
    width = 8,
    dpi = 700
    )

  ## Create Heatmap
  p.heatmap.top10 <- sc.top10.marker.heatmap.rc(
    # Seurat object
    d,
    # Cluster column
    "recluster",
    # Heatmap width
    36,
    # Heatmap height
    12,
    # Column fontsize
    4,
    # Row fontsize
    8
    )

  ## Save
  png(
    "analysis/recluster/plot.heatmap.top10.markers.png",
    width = 40,
    height = 14,
    units = "cm",
    res = 1000
  )
  print(p.heatmap.top10)
  dev.off()













  #START HERE 7/15/24

  #---- Cluster Counts ----

  ## Source function

  source(
    "Scripts/3.3.qc.cluster.violin.R",
    local = knitr::knit_global()
  )

  ## Output and save table

  cl.counts <- fun.cluster.counts(
    list.analysis.subset,
    "recluster",
    "Code",
    "Airway",
    "Knockout",
    "secretory"
  )

  cl.counts <- fun.prop.generic(
    list.analysis.subset@meta.data,
    "recluster",
    c("Code",
      "Airway",
      "Knockout")
  )

  write.table(
    cl.counts,
    paste(
      "Analysis/Integrated/",
      "D28.integrated.prop.code.airway.knockout.secretory.only.txt",
      sep = ""
    ),
    col.names = T,
    row.names = F,
    sep = "\t"
  )


  #---- Single UMAP with Multiple Metadata Input ----

  # Load reclustered data (if not already present)

  list.analysis.subset <- readRDS("Analysis/RDS/D28.original.analysis.secretory.rds")

  ## Add custom label

  ### Re-grouping of clusters

  cl.custom <- data.frame(
    "recluster" = seq.int(
      1,18,1
    ),
    "Group" = factor(
      c(
        "Secretory.set1","Secretory.set1","SecretoryCiliated","SecretoryCiliated","Secretory.set1",
        "SecretoryCiliated","Secretory.set1","Secretory.set2","Secretory.set2","SecretoryCiliated",
        "Secretory.set2","Secretory.set2","12.Secretory","Secretory.set1","SecretoryCiliated",
        "Secretory.set1","16.SecretoryCiliated","SecretoryCiliated"),
      levels = c(
        "12.Secretory","16.SecretoryCiliated",
        "Secretory.set1","Secretory.set2","SecretoryCiliated"
      )
    )
  )

  cl.custom2 <- setNames(
    as.data.frame(
      as.numeric(
        list.analysis.subset@meta.data$recluster
      )
    ),
    c(
      "recluster"
    )
  )


  cl.custom3 <- dplyr::left_join(
    cl.custom2,
    cl.custom,
    by = "recluster"
  )

  list.analysis.subset <- AddMetaData(
    list.analysis.subset,
    metadata = cl.custom3$Group,
    col.name = "Group2"
  )



  # Save plot (either 3D or 2D)

  ggsave(
    paste(
      list.p$a.path,
      "umap.secretory.airway.2D.png",
      sep = ""
    ),
    fun.sc.umap.plot.single.2D(
      # Seurat object
      list.analysis.subset,
      # Variable for labeling clusters
      "recluster",
      # Variable for fill color
      "Group2",
      # Color scheme
      ggsci::pal_npg("nrc")(10)
    ),
    height = 8,
    width = 10,
    dpi = 700
  )



  DefaultAssay(list.analysis.subset) <- "RNA"


  ken.list <- c(
    # Secretory Cell Markers (Orig.recluster # 12)
    "SFTPB","HLA-DPB1","LTF",
    "KLK13","HP","SCGB3A2",
    "HLA-DQA1","HLA-DPA1","RNASE1",
    "HLA-DQA2",
    # Goblet Cell Markers
    "ITLN1","FOXA3","MUC5AC",
    "MUC5B","SPDEF","XBP1",
    "TF",
    # Extra markers
    "TMEM45A","ATOH8","CP",
    "CFTR","FOXI1","SCGB1A1",
    "FOXJ1","KRT5","NKX2-1","IFI27","IL1R1",
    "IL1A","IL1B","IL16",
    "ICAM1","CEACAM6","SERPINB3","BSND",
    "ASCL3","CLCNKB","ATP6V1C2",
    "POU2F3","ALOX5","IL25",
    "GRP","ASCL1","CALCA"
  )


  ggsave(
    paste(
      list.p$a.path,
      "D28.integrated.umap.secretory.subs.exp.TF.png",
      sep = ""
    ),
    fun.sc.umap.plot.panel2.2D(
      # Seurat Object
      list.analysis.subset,
      # Group variable
      "recluster",
      # Gene name
      "TF",
      # Color scheme
      col1[1:length(
        levels(
          list.analysis.subset@meta.data[["recluster"]]
        )
      )],
      # Color Names
      c(
        levels(
          list.analysis.subset@meta.data[["recluster"]]
        )
      ),
      # legend x-position
      0.95,
      # legend y-position
      0.95
    ),
    height = 8,
    width = 20,
    dpi = 700
  )

  lapply(ken.list,
         function(x)
           ggsave(
             paste(
               "Analysis/Integrated/20240603.D28.integrated.exp.umaps/",x,".png",
               sep = ""
             ),
             fun.sc.umap.plot.panel2.2D(
               # Seurat Object
               list.analysis,
               # Group variable
               "CellType",
               # Gene name
               x,
               # Color scheme
               col1[1:length(
                 levels(
                   list.analysis.subset@meta.data[["CellType"]]
                 )
               )],
               # Color Names
               c(
                 levels(
                   list.analysis.subset@meta.data[["CellType"]]
                 )
               ),
               # legend x-position
               0.95,
               # legend y-position
               0.95
             ),
             height = 10,
             width = 24,
             dpi = 300
           ))


}
























