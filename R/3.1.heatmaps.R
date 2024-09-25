#' Top-10 Marker Gene Heatmap (Reclustered)
#'
#' Generates a heatmap from a reclustered Seurat Object and marker gene list based on the top-10 marker genes for each cluster.
#'
#' @param so An object of class Seurat.
#' @param cl.var Character string containing the name of the cluster variable for cell type predictions.
#' @param h.w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param h.h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs.c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs.r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @return A ComplexHeatmap object containing a reclustered top-10 marker gene heatmap.
#' @examples
#'
#' # p.umap <- sc.top10.marker.heatmap(d.annotated,"seurat.clusters",18,24,6,8)
#'
#' @export
  sc.top10.marker.heatmap.rc <- function(
    so,
    cl.var,
    h.w,
    h.h,
    fs.c,
    fs.r
  ) {
    d <- so
    if(!file.exists("analysis/recluster/table.marker.genes.txt")) {
      print("No marker gene file has been created; calculating marker genes for each cluster...")
      cl.mark <- Seurat::FindAllMarkers(d,verbose = T)
      write.table(
        cl.mark,
        "analysis/recluster/table.marker.genes.txt",
        col.names = T,
        row.names = F,
        sep = "\t"
      )
    }
    ## Marker gene input matrix (top10 per cell type)
    d.mark <- read.table(
      "analysis/recluster/table.marker.genes.txt",
      sep = "\t",
      header = T
    )
    d.mark[["CellType.no"]] <- d.mark[["cluster"]]
    ### Top 10 genes per cluster
    d.mark <- dplyr::slice_max(
      dplyr::group_by(
        d.mark,
        .data[["CellType.no"]]),
      order_by = .data[["avg_log2FC"]],
      n = 10
    )[,c(
      "gene",
      "cluster"
    )]
    #### Save table
    write.table(
      d.mark,
      "analysis/recluster/table.marker.genes.top10.txt",
      row.names = F,
      col.names = T,
      sep = "\t"
    )
    ### Subset seurat and scale
    SeuratObject::DefaultAssay(d) <- "RNA"
    h <- SeuratObject::FetchData(
      d,
      vars = c(
        cl.var,
        unique(
          d.mark[["gene"]]
        )
      )
    )
    
    ### Heatmap annotation (average expression)
    h.anno <- as.data.frame(
      lapply(
        h[,2:ncol(
          h
        )],
        function(x)
          mean(x)
      )
    )
    
    h.anno <- h.anno[,h.anno[1,] > 0]
    ### Scale and plot average expression per cell type
    h.in <- scale(
      as.matrix(
        magrittr::set_rownames(
          setNames(
            as.data.frame(
              lapply(
                h[,2:ncol(
                  h
                )],
                function(x)
                  dplyr::select(
                    aggregate(
                      x,
                      list(
                        h[,1]
                      ),
                      FUN = mean
                    ),
                    c(2)
                  )
              )
            ),
            names(h[,2:ncol(h)])
          ),
          levels(h[,1])
        )
      ),
      center = T
    )
    qs <- quantile(
      h.in,
      probs = c(
        0.05,
        0.95
      ),
      na.rm = T
      )
    
    h.in <- as.matrix(
      as.data.frame(h.in)[,unlist(
        lapply(
          seq.int(1,ncol(as.data.frame(h.in)),1),
          function(x) 
            !anyNA(as.data.frame(h.in)[x])))])
    
    fun.hm.col <- circlize::colorRamp2(
      c(
        qs[[1]],
        (qs[[1]])/2,
        (qs[[2]])/2,
        qs[[2]]
      ),
      colors = col.grad()[c(
        1,3,
        6,12
      )]
    )
    # Create Plot
    h.out <- ComplexHeatmap::Heatmap(
      h.in,
      col = fun.hm.col,
      name = "Scaled Expression",
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        `Average.Expression` = ComplexHeatmap::anno_barplot(
          as.vector(t(h.anno)),
          gp = grid::gpar(fill = col.grad())
        ),
        annotation_name_gp = grid::gpar(
          fontsize = 10
        )
      ),
      show_column_names = T,
      show_row_names = T,
      heatmap_width = ggplot2::unit(h.w,"cm"),
      heatmap_height = ggplot2::unit(h.h,"cm"),
      column_title = "Marker Genes (Top 10)",
      column_names_rot = 90,
      column_names_gp = grid::gpar(fontsize = fs.c),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = fs.r),
      cluster_columns = F,
      cluster_rows = F
    )
    return(h.out)
  }







#' Top-10 Marker Gene Heatmap
#'
#' Generates a heatmap from a Seurat Object and marker gene list based on the top-10 marker genes for each cluster.
#'
#' @param so An object of class Seurat.
#' @param cl.var Character string containing the name of the cluster variable for cell type predictions.
#' @param h.w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param h.h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs.c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs.r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @return A ComplexHeatmap object containing a top-10 marker gene heatmap.
#' @examples
#'
#' # p.umap <- sc.top10.marker.heatmap(d.annotated,"seurat.clusters",18,24,6,8)
#'
#' @export
sc_top10_marker_heatmap <- function(
  so,
  asy,
  cl.var,
  h.w,
  h.h,
  fs.c,
  fs.r
  ) {
  d <- so
  if(!file.exists("analysis/table.marker.genes.txt") && asy == "RNA") {
    print("No marker gene file has been created; calculating marker genes for each cluster...")
    Seurat::DefaultAssay(d) <- "RNA"
    cl.mark <- Seurat::FindAllMarkers(d,verbose = T)
    write.table(
      cl.mark,
      "analysis/table.marker.genes.txt",
      col.names = T,
      row.names = F,
      sep = "\t"
      )
    }

  if(!file.exists("analysis/table.marker.motifs.txt") && asy == "ATAC") {
    print("No marker motif file has been created; calculating marker motifs for each cluster...")
    Seurat::DefaultAssay(d) <- "ATAC"
    cl.mark <- Seurat::FindAllMarkers(d,min.pct = 0.05,verbose = T)
    cl.mark <- dplyr::left_join(
      cl.mark,
      names_motif,
      by = "gene"
    )
    write.table(
      cl.mark,
      "analysis/table.marker.motifs.txt",
      col.names = T,
      row.names = F,
      sep = "\t"
      )
    }

    if(asy == "RNA") {
    cl.mark <- read.table(
      "analysis/table.marker.genes.txt",
      sep = "\t",
      header = T
      )
    }

    if(asy == "ATAC") {
    cl.mark <- read.table(
      "analysis/table.marker.motifs.txt",
      sep = "\t",
      header = T
      )
    }
    ## Marker gene input matrix (top10 per cell type)
    if(class(cl.mark[["cluster"]]) == "character"){
      d.mark <- d.mark[gtools::mixedorder(d.mark[["cluster"]]),]
      }
    cl.mark[["CellType.no"]] <- cl.mark[["cluster"]]
    
    cl.mark <- dplyr::group_by(
      cl.mark,
      .data[["CellType.no"]]
      )
    
    ### Top 10 genes per cluster (by p value then by fold change)
    cl.mark <- dplyr::slice_max(
      cl.mark[cl.mark[["avg_log2FC"]] > 0,],
      order_by = -.data[["p_val_adj"]],
      n = 25
      )[,c(
      "gene",
      "cluster",
      "avg_log2FC",
      "p_val_adj"
      )]
    
    cl.mark <- dplyr::group_by(
      cl.mark,
      .data[["cluster"]]
      )
    
    cl.mark <- dplyr::slice_max(
      cl.mark,
      order_by = .data[["avg_log2FC"]],
      n = 10
      )[,c(
      "gene",
      "cluster"
      )]

    #### Save table
    if(asy == "RNA") {
    write.table(
      cl.mark,
      "analysis/table.marker.genes.top10.txt",
      row.names = F,
      col.names = T,
      sep = "\t"
      )
    SeuratObject::DefaultAssay(d) <- "RNA"
    }
    if(asy == "ATAC") {
    write.table(
      cl.mark,
      "analysis/table.marker.motifs.top10.txt",
      row.names = F,
      col.names = T,
      sep = "\t"
      )
    SeuratObject::DefaultAssay(d) <- "ATAC"
    }
    ### Subset seurat and scale
    h <- SeuratObject::FetchData(
      d,
      vars = c(
        cl.var,
        unique(
          cl.mark[["gene"]]
          )
        )
      )
    ### Heatmap annotation (average expression)
    h.anno <- as.data.frame(
      lapply(
        h[,2:ncol(
          h
        )],
        function(x)
          mean(x)
        )
      )
    
    h.anno <- h.anno[,h.anno[1,] > 0]
    ### Scale and plot average expression/accessibility per cell type
    h.in <- scale(
      as.matrix(
        magrittr::set_rownames(
          setNames(
            as.data.frame(
              lapply(
                h[,2:ncol(
                  h
                  )],
                function(x)
                  dplyr::select(
                    aggregate(
                      x,
                      list(
                        h[,1]
                        ),
                      FUN = mean
                      ),
                    c(2)
                    )
                )
              ),
            names(h[,2:ncol(h)])
            ),
          levels(h[,1])
          )
        ),
      center = T
      )
    qs <- quantile(
      h.in,
      probs = c(
        0.05,
        0.95
        ),
      na.rm = T
      )
    
    h.in <- as.matrix(
      as.data.frame(h.in)[,unlist(
        lapply(
          seq.int(1,ncol(as.data.frame(h.in)),1),
          function(x) 
            !anyNA(as.data.frame(h.in)[x])))])
    
    fun.hm.col <- circlize::colorRamp2(
      c(
        qs[[1]],
        (qs[[1]])/2,
        (qs[[2]])/2,
        qs[[2]]
        ),
      colors = col_grad()[c(
        1,3,
        6,12
        )]
      )
    # Create Plot
    h.out <- ComplexHeatmap::Heatmap(
      h.in,
      col = fun.hm.col,
      name = "Scaled Expression/Accessibility",
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        `Average.Expression` = ComplexHeatmap::anno_barplot(
          as.vector(t(h.anno)),
          gp = grid::gpar(fill = col_grad())
          ),
        annotation_name_gp = grid::gpar(
          fontsize = 10
          )
        ),
      show_column_names = T,
      show_row_names = T,
      heatmap_width = ggplot2::unit(h.w,"cm"),
      heatmap_height = ggplot2::unit(h.h,"cm"),
      column_title = "Top 10 Markers",
      column_names_rot = 90,
      column_names_gp = grid::gpar(fontsize = fs.c),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = fs.r),
      cluster_columns = F,
      cluster_rows = F
      )
    return(h.out)
  }





#' Top-10 DEG Heatmap
#'
#' Generates a heatmap from a Seurat Object and marker gene list based on the top-10 marker genes for each cluster.
#'
#' @param l.deg A list of DGEA results returned by sc.DGEA().
#' @param so An object of class Seurat.
#' @param asy Character string providing the name of the assay to use.
#' @param cl.var Character string containing the name of the clustering variable.
#' @param hm.w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param hm.h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs.c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs.r Numeric value for row fontsize (passed to ComplexHeatmap).
#' @param c.name Comparison name for plot title, provided as a character string.
#' @return A ComplexHeatmap object containing a top-10 marker gene heatmap.
#' @examples
#'
#' # p.heatmap <- sc.top10.deg.heatmap(dgea.output,d.annotated,"seurat.clusters",18,24,6,8)
#'
#' @export
sc.top10.de.da.heatmap <- function(
  l.deg,
  so,
  asy,
  cl.var,
  hm.w,
  hm.h,
  fs.c,
  fs.r,
  c.name
  ) {
  d1 <- so
  SeuratObject::DefaultAssay(d1) <- asy
  ld <- l.deg[[1]]

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
      order_by = -.data[["H.pval"]],
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
  mt.rp.remove <- TRUE
  ifelse(
    mt.rp.remove == TRUE,
    d <- d[!grepl("MT-|RP",d[["GENE"]]),],
    d
    )
  ## Filter DGEA results based on gene list
  d2 <- dplyr::select(
    ld[ld[["Comparison"]] %in%
                  d[["Comparison"]] &
                  ld[["GENE"]] %in%
                  unique(
                    d[["GENE"]]
                  ),],
    c(
      "CellType",
      "Comparison",
      "GENE",
      "logFC"
      )
    )
  ## Format as matrix
  d3 <- reshape2::dcast(
    d2[,-c(5)],
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
      function(x)
        ifelse(
          is.na(x),
          0,
          x
          )
      )
    )
  ## Row annotations
  h.row.anno3 <- d3[["CellType"]]
  ## Input matrix
  d3 <-magrittr::set_rownames(
    as.matrix(d3[3:ncol(d3)]),
    d3[["Comparison"]]
    )
  ### Subset seurat and scale
  h <-SeuratObject::FetchData(
    d1,
    vars = c(
      cl.var,
      unique(
        d2[["GENE"]]
        )
      )
    )
  ### Heatmap annotation (average expression)
  h.anno <- as.data.frame(
    lapply(
      h[,2:ncol(
        h
      )],
      function(x)
        mean(x)
      )
    )
  ## Color scheme
  fun.hm.col <- circlize::colorRamp2(
    c(
      min(d3),min(d3)/2,
      0,max(d3)/2,
      max(d3)
    ),
    colors = c(
      col.grad()[[1]],col.grad()[[3]],
      "white",col.grad()[[6]],
      col.grad()[[12]]
      )
    )
  ## Output plot
  h.out <- ComplexHeatmap::Heatmap(
    d3,
    col = fun.hm.col,
    name = "logFC",
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      `Avg.Exp` = ComplexHeatmap::anno_barplot(
        as.vector(t(h.anno)),
        gp = grid::gpar(
          fill = col.grad()
          )
        ),
      annotation_name_gp = grid::gpar(
        fontsize = 8
        ),
      annotation_name_side = "left"
      ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      "Cluster" = h.row.anno3,
      col = list(
        "Cluster" = setNames(
          as.vector(
            col.univ()[1:length(
              unique(h.row.anno3)
              )]
            ),
          c(unique(h.row.anno3))
          )
        ),
      show_annotation_name = T,
      annotation_name_gp = grid::gpar(fontsize = 10)
      ),
    show_column_names = T,
    show_row_names = F,
    row_names_gp = grid::gpar(
      fontsize = fs.r
      ),
    heatmap_width = ggplot2::unit(
      hm.w,"cm"
      ),
    heatmap_height = ggplot2::unit(
      hm.h,"cm"
      ),
    column_title = paste("Top 10 DEG per Cell Type",c.name,sep = ": "),
    column_names_rot = 90,
    column_names_gp = grid::gpar(
      fontsize = fs.c
      ),
    cluster_columns = T,
    cluster_rows = T
    )

  return(h.out)
  }





#' Top-10 DEG Dot Plot
#'
#' Generates a dot plot from a Seurat Object and DEG list based on the top-10 DEGs for a specific cell type and comparison.
#'
#' @param p.type Should a custom gene list or DGEA results object be used for plotting? Type "cstm.list" for custom gene lists and "deg.list"
#' for using dgea.results objects.
#' @param l.deg A list of DGEA results returned by sc.DGEA() or a vector of selected genes for plotting.
#' @param so An object of class Seurat.
#' @param ct A pattern provided as a character string for matching a specific cell type or types.
#' @param cl.var Character string indicating the name of the cluster/cell type variable.
#' @param list.var A vector of character strings indicating the name(s) of up to two group variables for stratifying plot points.
#' @return An input data frame and corresponding dot plot displaying the expression of the top-10 DEGs for a specific cell type.
#' @examples
#'
#' # p.dotplot <- sc.top10.deg.dotplot(
#' # # Type "deg.list" or "cstm.list" to toggle between inputs
#' # "deg.list",
#' # # Name of a custom gene list or dgea.results object
#' # dgea.output,
#' # # Seurat object
#' # d.annotated,
#' # # Unique character strings corresponding to cell types
#' # "3.Se|6.Se",
#' # # Name of clustering variable
#' # "CellType",
#' # # Vector of up to 2 variables for stratifying clustering variables
#' # c("Knockout","Airway")
#' # )
#'
#' @export
sc.top10.deg.dotplot <- function(
  p.type,
  l.deg,
  so,
  ct,
  cl.var,
  list.var
  ) {
  # Load Seurat and change default assay to RNA
  d <- so
  SeuratObject::DefaultAssay(d) <- "RNA"
  c <- ct

  if(p.type == "deg.list"){
  ## Return specified number of DEGs for selected comparison and cell types
  ld <- l.deg[["DGEA.results"]]

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
      order_by = -.data[["H.pval"]],
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
  top10 <- top10[grepl(c,top10[["CellType"]]),]
  ## return list of unique genes
  top10.pres <- unique(top10[["GENE"]][top10[["GENE"]] %in% SeuratObject::Features(d)])
  top10.abs <- subset(top10[["GENE"]], !(top10[["GENE"]] %in% SeuratObject::Features(d)))
  }

  if(p.type == "cstm.list"){
    top10 <- as.vector(l.deg[[1]])
    top10.pres <- unique(top10[top10 %in% SeuratObject::Features(d)])
    top10.abs <- subset(top10, !(top10 %in% SeuratObject::Features(d)))
    }

  ## select genes from Seurat object
  d1 <- cbind(
    SeuratObject::FetchData(
      d,
      vars = c(
        cl.var,
        list.var,
        top10.pres
        )
      )
    )
  ## Subset based on cell type
  d1 <- d1[grepl(c,d1[[cl.var]]),]
  d1[["CellType"]] <- factor(
    as.character(d1[[cl.var]]),
    levels = gtools::mixedsort(
      unique(
        as.character(d1[[cl.var]])
        )
      )
    )
  ## Count/ratio table for creating dot plots
  d1.prc <- dplyr::bind_rows(
    setNames(
      lapply(
        top10.pres,
        function(x)
        {
          # Determine average expression of each gene per cell type and comparison
          ## for 2 variables:
          if(length(list.var) == 2) {
            ### average expression
            d.avg <- setNames(
              aggregate(
                d1[[x]],
                by = list(
                  d1[,c("CellType")],
                  d1[,c(list.var[[1]])],
                  d1[,c(list.var[[2]])]
                  ),
                function(y)
                  mean(y)
                ),
              c(
                "CellType",c(list.var),
                "avg.exp"
                )
              )
            d.avg[["avg.exp"]] <- round(
              d.avg[["avg.exp"]],
              digits = 2
              )
            ### percent expressed
            d.prc <- dplyr::count(
              d1,.data[["CellType"]],
              .data[[list.var[[1]]]],
              .data[[list.var[[2]]]],
              .data[[x]] > 0
              )
            d.prc <- setNames(
              dplyr::filter(
                d.prc,d.prc[4] == T
                ),
              c(
                "CellType",c(list.var),
                "pres","n"
                )
              )
            d.cnt <- setNames(
              dplyr::count(
                d1,.data[["CellType"]],
                .data[[list.var[[1]]]],
                .data[[list.var[[2]]]]
                ),
              c(
                "CellType",c(list.var),
                "n"
                )
              )
            d.comb <- dplyr::left_join(
              d.cnt,
              d.prc,
              by = c(
                "CellType",c(list.var)
                )
              )
            d.comb[is.na(d.comb)] <- 0
            d.comb <- dplyr::mutate(
              d.comb,
              "perc.exp" = round(
                d.comb[["n.y"]]/
                  d.comb[["n.x"]],
                digits = 2
                )
              )
            ### combined
            d.comb.out <- dplyr::left_join(
              d.avg,
              d.comb[,
                     c(
                       "CellType",c(list.var),
                       "perc.exp"
                     )],
              by = c(
                "CellType",
                c(list.var)
              )
            )
          }

          if(length(list.var) < 2) {
            ### average expression
            d.avg <- setNames(
              aggregate(
                d1[[x]],
                by = list(
                  d1[,c("CellType")],
                  d1[,c(list.var[[1]])]
                ),
                function(y)
                  mean(y)
              ),
              c(
                "CellType",c(list.var),
                "avg.exp"
              )
            )
            d.avg[["avg.exp"]] <- round(
              d.avg[["avg.exp"]],
              digits = 2
            )
            ### percent expressed
            d.prc <- dplyr::count(
              d1,.data[["CellType"]],
              .data[[list.var[[1]]]],
              .data[[x]] > 0
            )
            d.prc <- setNames(
              dplyr::filter(
                d.prc,d.prc[3] == T
              ),
              c(
                "CellType",c(list.var),
                "pres","n"
              )
            )
            d.cnt <- setNames(
              dplyr::count(
                d1,.data[["CellType"]],
                .data[[list.var[[1]]]]
              ),
              c(
                "CellType",c(list.var),
                "n"
              )
            )
            d.comb <- dplyr::left_join(
              d.cnt,
              d.prc,
              by = c(
                "CellType",c(list.var)
              )
            )
            d.comb[is.na(d.comb)] <- 0
            d.comb <- dplyr::mutate(
              d.comb,
              "perc.exp" = round(
                d.comb[["n.y"]]/
                  d.comb[["n.x"]],
                digits = 2
              )
            )
            ### combined
            d.comb.out <- dplyr::left_join(
              d.avg,
              d.comb[,
                     c(
                       "CellType",c(list.var),
                       "perc.exp"
                     )],
              by = c(
                "CellType",
                c(list.var)
              )
            )
          }



          return(d.comb.out)
        }
      ),
      c(top10.pres)
      ),
    .id = "GENE"
    )

  ## Add row name labels and convert GENE column to factor
  if(length(list.var) == 2) {

    d1.prc <- data.frame(
      d1.prc,
      "labs" = factor(
        paste(
          d1.prc[["CellType"]],
          d1.prc[[list.var[[1]]]],
          d1.prc[[list.var[[2]]]],
          sep = " "),
        levels = gtools::mixedsort(
          unique(
            paste(
              d1.prc[["CellType"]],
              d1.prc[[list.var[[1]]]],
              d1.prc[[list.var[[2]]]],
              sep = " "
              )
            )
          )
        )
      )

    d1.prc[["GENE"]] <- factor(
      d1.prc[["GENE"]],
      levels = unique(
        d1.prc[["GENE"]]
        )
      )
  }

## Add row name labels and convert GENE column to factor
if(length(list.var) < 2) {

  d1.prc <- data.frame(
    d1.prc,
    "labs" = factor(
      paste(
        d1.prc[["CellType"]],
        d1.prc[[list.var[[1]]]],
        sep = " "),
      levels = gtools::mixedsort(
        unique(
          paste(
            d1.prc[["CellType"]],
            d1.prc[[list.var[[1]]]],
            sep = " "
            )
          )
        )
      )
    )

  d1.prc[["GENE"]] <- factor(
    d1.prc[["GENE"]],
    levels = unique(
      d1.prc[["GENE"]]
      )
    )
  }



  ## Plot
  p.dot <- ggplot2::ggplot(
    d1.prc,
    ggplot2::aes(
      x = .data[["GENE"]],
      y = .data[["labs"]],
      fill = .data[["avg.exp"]],
      size = .data[["perc.exp"]]
      )
    ) +
    ggplot2::geom_point(
      shape = 21
      ) +
    ggplot2::scale_fill_gradientn(
      colors = col.grad()
      ) +
    ggplot2::scale_size_area(max_size = 10) +
    sc.theme1() +
    ggplot2::labs(
      fill = "Average Expression",
      size = "Percent Expressed",
      y = ""
      ) +
    ggplot2::theme(
      plot.margin = ggplot2::unit(
        c(.2,.2,
          .2,.2),
        "cm"
        )
      )

  if(exists("top10.abs") & length(top10.abs) == 1) {
    print(
      paste(
        top10.abs,
        "was not found in the provided Seurat object; plots for this gene was excluded from the dot plot...",
        sep = " "
      )
    )
  }

  if(exists("top10.abs") & length(top10.abs) > 1) {
    print(
      paste(
        top10.abs,
        "were not found in the provided Seurat object; these genes were excluded from the dot plot...",
        sep = " "
      )
    )
  }

  return(
    list(
      "Input" = d1.prc,
      "Plot" = p.dot
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
#' @param top_n Character string containing the name of the clustering variable.
#' @param hm.w Numeric value for heatmap width (passed to ComplexHeatmap).
#' @param hm.h Numeric value for heatmap height (passed to ComplexHeatmap).
#' @param fs.c Numeric value for column fontsize (passed to ComplexHeatmap).
#' @param fs.r Numeric value for row fontsize (passed to ComplexHeatmap).
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
sc_top_motif_heatmap <- function(
  so,
  mot_m,
  filt_ct,
  ct_name,
  cl_var,
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
  d_mark[["CellType.no"]] <- d_mark[["cluster"]]

  d_mark <- dplyr::group_by(
    d_mark,
    .data[["CellType.no"]] # nolint
  )

  ## Top motifs per cluster (by p value then by fold change)
  d_mark <- dplyr::slice_max(
    d_mark[d_mark[["avg_diff"]] > 0, ],
    order_by = -.data[["p_val_adj"]], # nolint
    n = top_n
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

  if(filt_ct == TRUE) { # nolint
    d_mark <- d_mark[d_mark[["cluster"]] == ct_name, ]

    #### Save table
    write.table(
      d_mark,
      paste(
        "analysis/table.marker.motifs.Top",
        paste(as.character(top_n)),
        ".",
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
  h <- SeuratObject::FetchData(
    d,
    vars = c(
      cl_var,
      unique(
        d_mark[["ID"]]
      )
    )
  )

  ### Scale and plot average expression/accessibility per cell type
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
    colors = col_grad()[c(
      1, 3,
      6, 12
    )
    ]
  )
  # Create Plot
  if(filt_ct == TRUE) { # nolint
    h_out <- ComplexHeatmap::Heatmap(
      h_in,
      col = fun_hm_col,
      name = "Scaled Activity Score",
      show_column_names = TRUE,
      show_row_names = TRUE,
      heatmap_width = ggplot2::unit(h_w, "cm"),
      heatmap_height = ggplot2::unit(h_h, "cm"),
      column_title = paste("Top", top_n, "Motifs:", ct_name),
      column_names_rot = 90,
      column_names_gp = grid::gpar(fontsize = fs_c),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = fs_r),
      cluster_columns = TRUE,
      cluster_rows = FALSE
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
      column_title = paste("Top", top_n, "Motifs"),
      column_names_rot = 90,
      column_names_gp = grid::gpar(fontsize = fs_c),
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = fs_r),
      cluster_columns = FALSE,
      cluster_rows = FALSE
    )
  }
  return(h_out)
}
