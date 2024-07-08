#### Heatmaps ####

#---- Top 10 marker genes based on Seurat object ----

  fun.hm.marker.gene <- function(
    lo,
    h.w,
    h.h,
    fs.c,
    fs.r
    ) {
  
    ## Marker gene input matrix (top10 per cell type)
    
    attach(lo)
    
    if(
      exists(
        "Marker Genes"
      )
    ) {
      
      d.mark <- lo[["Marker Genes"]]
      
    }
    
    
    if(
      !exists(
        "Marker Genes"
      )
    ) {
      
      d.mark <- read.table(
        paste(
          list.p$a.path,
          "marker.genes.top10.txt",
          sep = ""
        ),
        sep = "\t",
        header = T
      )
      
    }
    

    d.mark[["CellType.no"]] <- ifelse(
      class(
        d.mark[["cluster"]]
        ) == "factor",
      d.mark[["cluster"]],
      factor(
        d.mark[["cluster"]] +
          1,
        levels = c(
          seq.int(
            1,
            length(
              unique(
                d.mark[["cluster"]] +
                  1
                )
              ),
              1
            )
          )
        )
      )
    
    
    ### Add cell types (requires existing cluster assignments object)
    
    d.mark <- dplyr::left_join(
      d.mark,
      data.frame(
        setNames(
          unique(
            lo[["cluster.assignments"]]
          ),
          c(
            "CellType",
            "CellGroup"
          )
        ),
        "CellType.no" = factor(
          gsub(
            "\\..*",
            "",
            as.character(
              unique(
                lo[["cluster.assignments"]][[1]]
              )
            )
          ),
          levels = c(
            seq.int(
              1,
              length(
                as.character(
                  unique(
                    lo[["cluster.assignments"]][[1]]
                  )
                )
              )
              ,
              1
            )
          )
        )
      ),
      by = "CellType.no")
    
    ### Top 10 genes per cluster
    
    d.mark <- dplyr::slice_max(
      dplyr::group_by(
        d.mark,
        .data[["CellType.no"]]),
      order_by = .data[["avg_log2FC"]],
      n = 10
    )[,c(
      "gene",
      "CellType",
      "CellGroup"
    )]
    
    #### Save table
    
    write.table(
      d.mark,
      paste(
        list.p$a.path,
        "marker.genes.top10.txt",
        sep = ""
      ),
      row.names = F,
      col.names = T,
      sep = "\t"
    )
    
    
    ### Subset seurat and scale
    
    h <- FetchData(
      lo[["Predicted Clusters"]],
      vars = c(
        "CellType",
        "CellGroup",
        unique(
          d.mark[["gene"]]
        )
      )
    )
    
    
    ### Heatmap annotation (average expression)
    
    h.anno <- as.data.frame(
      lapply(
        h[,3:ncol(
          h
        )],
        function(x) 
          mean(x)
      )
    )
    
    
    ### Scale and plot average expression per cell type
    
    h.in <- scale(
      as.matrix(
        magrittr::set_rownames(
          setNames(
            as.data.frame(
              lapply(
                h[,3:ncol(
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
                    c(
                      2
                    )
                  )
              )
            ),
            names(
              h[,3:ncol(
                h
              )]
            )
          ),
          levels(
            h[,1]
          )
        )
      ),
      center = T
    )
    
    
    qs <- quantile(
      h.in,
      probs = c(
        0.05,
        0.95
      )
    )
    
    
    fun.hm.col <- colorRamp2(
      c(
        qs[[1]],
        (qs[[1]])/2,
        (qs[[2]])/2,
        qs[[2]]
      ),
      colors = col2a[c(
        1,3,
        6,12
      )]
    )
    
    
    
    
    
    list(
      "CellGroup" = setNames(
        as.vector(
          col3a
        ),
        c(
          gsub(
            "^.*\\.",
            "",
            row.names(
              h.in
            )
          )
        )
      ))
    
    
    
    
    
    
    
    h.out <- ComplexHeatmap::Heatmap(
      h.in,
      col = fun.hm.col,
      name = "Scaled Expression",
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        `Average.Expression` = ComplexHeatmap::anno_barplot(
          as.vector(
            t(
              h.anno
            )
          ),
          gp = gpar(
            fill = col1b
          )
        ),
        annotation_name_gp = gpar(
          fontsize = 10
        )
      ),
      left_annotation = ComplexHeatmap::rowAnnotation(
        "CellGroup" = gsub(
          "^.*\\.",
          "",
          row.names(
            h.in
          )
        ),
        col = list(
          "CellGroup" = setNames(
            as.vector(
              col3a[1:length(
                row.names(
                  h.in
                )
              )]
            ),
            c(
              gsub(
                "^.*\\.",
                "",
                row.names(
                  h.in
                )
              )
            )
          )
        ),
        show_annotation_name = T,
        annotation_name_gp = gpar(
          fontsize = 10
        )),
      show_column_names = T,
      show_row_names = T,
      heatmap_width = unit(
        h.w,"cm"
      ),
      heatmap_height = unit(
        h.h,"cm"
      ),
      column_title = "Cluster Marker Genes (Top 10)",
      column_names_rot = 90,
      column_names_gp = gpar(
        fontsize = fs.c
      ),
      row_names_side = "left",
      row_names_gp = gpar(
        fontsize = fs.r
      ),
      cluster_columns = F,
      cluster_rows = F
    )
    
    return(
      h.out
    )
    
  }
  




  #---- Top 10 marker genes based on Seurat object (Reclustered) ----
  
  fun.hm.marker.gene.recluster <- function(
    lo,
    h.w,
    h.h,
    fs.c,
    fs.r
  ) {
    
    ## Marker gene input matrix (top10 per cell type)
    
    attach(lo)
    
    if(
      exists(
        "Marker Genes"
      )
    ) {
      
      d.mark <- lo[["Marker Genes"]]
      
    }
    
    
    if(
      !exists(
        "Marker Genes"
      )
    ) {
      
      d.mark <- read.table(
        paste(
          list.p$a.path,
          "marker.genes.top10.txt",
          sep = ""
        ),
        sep = "\t",
        header = T
      )
      
    }
    
    
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
      paste(
        list.p$a.path,
        "marker.genes.top10.txt",
        sep = ""
      ),
      row.names = F,
      col.names = T,
      sep = "\t"
    )
    
    
    ### Subset seurat and scale
    
    DefaultAssay(lo[["Predicted Clusters"]]) <- "RNA"
    
    h <- FetchData(
      lo[["Predicted Clusters"]],
      vars = c(
        "recluster",
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
                    c(
                      2
                    )
                  )
              )
            ),
            names(
              h[,2:ncol(
                h
              )]
            )
          ),
          levels(
            h[,1]
          )
        )
      ),
      center = T
    )
    
    
    qs <- quantile(
      h.in,
      probs = c(
        0.05,
        0.95
      )
    )
    
    
    fun.hm.col <- colorRamp2(
      c(
        qs[[1]],
        (qs[[1]])/2,
        (qs[[2]])/2,
        qs[[2]]
      ),
      colors = col2a[c(
        1,3,
        6,12
      )]
    )
    
    
    
    
    h.out <- ComplexHeatmap::Heatmap(
      h.in,
      col = fun.hm.col,
      name = "Scaled Expression",
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        `Average.Expression` = ComplexHeatmap::anno_barplot(
          as.vector(
            t(
              h.anno
            )
          ),
          gp = gpar(
            fill = col1b
          )
        ),
        annotation_name_gp = gpar(
          fontsize = 10
        )
      ),
      show_column_names = T,
      show_row_names = T,
      heatmap_width = unit(
        h.w,"cm"
      ),
      heatmap_height = unit(
        h.h,"cm"
      ),
      column_title = "Recluster Marker Genes (Top 10)",
      column_names_rot = 90,
      column_names_gp = gpar(
        fontsize = fs.c
      ),
      row_names_side = "left",
      row_names_gp = gpar(
        fontsize = fs.r
      ),
      cluster_columns = F,
      cluster_rows = F
    )
    
    return(
      h.out
    )
    
  }
  
  
  
  
  
## Alternate without cluster.assignments object
  
  fun.hm.marker.gene.alt <- function(
    so,
    mark.gene,
    h.w,
    h.h,
    fs.c,
    fs.r
  ) {
    
    ## Marker gene input matrix (top10 per cell type)

      d.mark <- read.table(
        paste(
          list.p$a.path,
          "marker.genes.txt",
          sep = ""
        ),
        sep = "\t",
        header = T
      )
    
    
      ifelse(
        class(d.mark[["cluster"]]) == "factor",
        d.mark[["CellType.no"]] <- d.mark[["cluster"]],
        d.mark[["CellType.no"]] <- factor(
          d.mark[["cluster"]] +
            1,
          levels = c(
              seq.int(
                1,
                length(
                  unique(
                    d.mark[["cluster"]] + 1
                    )
                  ),
                1
                )
              )
          )
        )
  
        
    ### Add cell types (requires existing cluster assignments object)
    
    d.mark <- dplyr::left_join(
      d.mark,
      data.frame(
          unique(
            so@meta.data[,c(
              "CellType",
              "CellGroup"
              )]
            ),
        "CellType.no" = factor(
          gsub(
            "\\..*",
            "",
            as.character(
              unique(
                so@meta.data[,c(
                  "CellType"
                )]
              )
            )
          ),
          levels = c(
            seq.int(
              1,
              length(
                as.character(
                  unique(
                    so@meta.data[,c(
                      "CellType"
                    )]
                  )
                )
              )
              ,
              1
            )
          )
        )
      ),
      by = "CellType.no")
    
    
    ### Top 10 genes per cluster
    
    d.mark <- dplyr::slice_max(
      dplyr::group_by(
        d.mark,
        .data[["CellType.no"]]),
      order_by = .data[["avg_log2FC"]],
      n = 10
    )[,c(
      "gene",
      "CellType",
      "CellGroup"
    )]
    
    #### Save table
    
    write.table(
      d.mark,
      paste(
        list.p$a.path,
        "marker.genes.top10.txt",
        sep = ""
      ),
      row.names = F,
      col.names = T,
      sep = "\t"
    )
    
    
    ### Subset seurat and scale
    
    h <- FetchData(
      so,
      vars = c(
        "CellType",
        "CellGroup",
        unique(
          d.mark[["gene"]]
        )
      )
    )
    
    
    ### Heatmap annotation (average expression)
    
    h.anno <- as.data.frame(
      lapply(
        h[,3:ncol(
          h
        )],
        function(x) 
          mean(x)
      )
    )
    
    
    ### Scale and plot average expression per cell type
    
    h.in <- scale(
      as.matrix(
        magrittr::set_rownames(
          setNames(
            as.data.frame(
              lapply(
                h[,3:ncol(
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
                    c(
                      2
                    )
                  )
              )
            ),
            names(
              h[,3:ncol(
                h
              )]
            )
          ),
          levels(
            h[,1]
          )
        )
      ),
      center = T
    )
    
    
    qs <- quantile(
      h.in,
      probs = c(
        0.05,
        0.95
      )
    )
    
    
    fun.hm.col <- colorRamp2(
      c(
        qs[[1]],
        (qs[[1]])/2,
        (qs[[2]])/2,
        qs[[2]]
      ),
      colors = col2a[c(
        1,3,
        6,12
      )]
    )
    
    
    
    h.out <- ComplexHeatmap::Heatmap(
      h.in,
      col = fun.hm.col,
      name = "Scaled Expression",
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        `Average.Expression` = ComplexHeatmap::anno_barplot(
          as.vector(
            t(
              h.anno
            )
          ),
          gp = gpar(
            fill = col1b
          )
        ),
        annotation_name_gp = gpar(
          fontsize = 10
        )
      ),
      left_annotation = ComplexHeatmap::rowAnnotation(
        "CellGroup" = gsub(
          "^.*\\.",
          "",
          row.names(
            h.in
          )
        ),
        col = list(
          "CellGroup" = setNames(
            as.vector(
              col3a[1:length(
                row.names(
                  h.in
                )
              )]
            ),
            c(
              gsub(
                "^.*\\.",
                "",
                row.names(
                  h.in
                )
              )
            )
          )
        ),
        show_annotation_name = T,
        annotation_name_gp = gpar(
          fontsize = 10
        )),
      show_column_names = T,
      show_row_names = T,
      heatmap_width = unit(
        h.w,"cm"
      ),
      heatmap_height = unit(
        h.h,"cm"
      ),
      column_title = "Cluster Marker Genes (Top 10)",
      column_names_rot = 90,
      column_names_gp = gpar(
        fontsize = fs.c
      ),
      row_names_side = "left",
      row_names_gp = gpar(
        fontsize = fs.r
      ),
      cluster_columns = F,
      cluster_rows = F
    )
    
    return(
      h.out
    )
    
  }
  
  
  
  
  
  #---- Top 10 marker genes based on Seurat object (Reclustered) ----
  
  fun.hm.marker.gene.recluster <- function(
    lo,
    mark.genes,
    h.w,
    h.h,
    fs.c,
    fs.r
  ) {
    
    ## Marker gene input matrix (top10 per cell type)
      
      d.mark <- read.table(
        paste(
          "Analysis/Integrated/",
          mark.genes,
          sep = ""
        ),
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
      paste(
        "Analysis/Integrated/",
        "marker.genes.top10.txt",
        sep = ""
      ),
      row.names = F,
      col.names = T,
      sep = "\t"
    )
    
    
    ### Subset seurat and scale
    
    DefaultAssay(lo[["Predicted Clusters"]]) <- "RNA"
    
    h <- FetchData(
      lo[["Predicted Clusters"]],
      vars = c(
        "recluster",
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
                    c(
                      2
                    )
                  )
              )
            ),
            names(
              h[,2:ncol(
                h
              )]
            )
          ),
          levels(
            h[,1]
          )
        )
      ),
      center = T
    )
    
    
    qs <- quantile(
      h.in,
      probs = c(
        0.05,
        0.95
      )
    )
    
    
    fun.hm.col <- colorRamp2(
      c(
        qs[[1]],
        (qs[[1]])/2,
        (qs[[2]])/2,
        qs[[2]]
      ),
      colors = col2a[c(
        1,3,
        6,12
      )]
    )
    
    
    
    
    h.out <- ComplexHeatmap::Heatmap(
      h.in,
      col = fun.hm.col,
      name = "Scaled Expression",
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        `Average.Expression` = ComplexHeatmap::anno_barplot(
          as.vector(
            t(
              h.anno
            )
          ),
          gp = gpar(
            fill = col1b
          )
        ),
        annotation_name_gp = gpar(
          fontsize = 10
        )
      ),
      show_column_names = T,
      show_row_names = T,
      heatmap_width = unit(
        h.w,"cm"
      ),
      heatmap_height = unit(
        h.h,"cm"
      ),
      column_title = "Recluster Marker Genes (Top 10)",
      column_names_rot = 90,
      column_names_gp = gpar(
        fontsize = fs.c
      ),
      row_names_side = "left",
      row_names_gp = gpar(
        fontsize = fs.r
      ),
      cluster_columns = F,
      cluster_rows = F
    )
    
    return(
      h.out
    )
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ## Marker gene input matrix (top10 per cell type)
  
  attach(list.analysis.subset2)
  
  if(
    exists(
      "Marker Genes"
    )
  ) {
    
    d.mark <- list.analysis.subset2[["Marker Genes"]]
    
  }
  
  
  if(
    !exists(
      "Marker Genes"
    )
  ) {
    
    d.mark <- read.table(
      paste(
        list.p$a.path,
        "marker.genes.top10.txt",
        sep = ""
      ),
      sep = "\t",
      header = T
    )
    
  }
  
  
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
    paste(
      list.p$a.path,
      "marker.genes.top10.txt",
      sep = ""
    ),
    row.names = F,
    col.names = T,
    sep = "\t"
  )
  
  
  ### Subset seurat and scale
  
  DefaultAssay(list.analysis.subset2$`Predicted Clusters`) <- "RNA"
  
  h <- FetchData(
    list.analysis.subset2[["Predicted Clusters"]],
    vars = c(
      "recluster",
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
                  c(
                    2
                  )
                )
            )
          ),
          names(
            h[,2:ncol(
              h
            )]
          )
        ),
        levels(
          h[,1]
        )
      )
    ),
    center = T
  )
  
  
  qs <- quantile(
    h.in,
    probs = c(
      0.05,
      0.95
    )
  )
  
  
  fun.hm.col <- colorRamp2(
    c(
      qs[[1]],
      (qs[[1]])/2,
      (qs[[2]])/2,
      qs[[2]]
    ),
    colors = col2a[c(
      1,3,
      6,12
    )]
  )
  
  
  
  
  
  list(
    "CellGroup" = setNames(
      as.vector(
        col3a
      ),
      c(
        gsub(
          "^.*\\.",
          "",
          row.names(
            h.in
          )
        )
      )
    ))
  
  
  
  
  
  
  
  
  
  #---- Top 10 marker genes based on Seurat object (Integrated) ----
  
  fun.hm.top10.mark <- function(
    so,
    mark.genes,
    var.clus,
    h.w,
    h.h,
    fs.c,
    fs.r
  ) {
    
    ## Marker gene input matrix (top10 per cell type)
    
    d.mark <- read.table(
      paste(
        list.p$a.path,
        mark.genes,
        sep = ""
      ),
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
      paste(
        list.p$a.path,
        "marker.genes.top10.txt",
        sep = ""
      ),
      row.names = F,
      col.names = T,
      sep = "\t"
    )
    
    
    ### Subset seurat and scale
    
    DefaultAssay(so) <- "RNA"
    
    h <- FetchData(
      so,
      vars = c(
        var.clus,
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
                    c(
                      2
                    )
                  )
              )
            ),
            names(
              h[,2:ncol(
                h
              )]
            )
          ),
          levels(
            h[,1]
          )
        )
      ),
      center = T
    )
    
    
    qs <- quantile(
      h.in,
      probs = c(
        0.05,
        0.95
      )
    )
    
    
    fun.hm.col <- colorRamp2(
      c(
        qs[[1]],
        (qs[[1]])/2,
        (qs[[2]])/2,
        qs[[2]]
      ),
      colors = col2a[c(
        1,3,
        6,12
      )]
    )
    
    
    
    
    h.out <- ComplexHeatmap::Heatmap(
      h.in,
      col = fun.hm.col,
      name = "Scaled Expression",
      top_annotation = ComplexHeatmap::HeatmapAnnotation(
        `Average.Expression` = ComplexHeatmap::anno_barplot(
          as.vector(
            t(
              h.anno
            )
          ),
          gp = gpar(
            fill = col1b
          )
        ),
        annotation_name_gp = gpar(
          fontsize = 10
        )
      ),
      show_column_names = T,
      show_row_names = T,
      heatmap_width = unit(
        h.w,"cm"
      ),
      heatmap_height = unit(
        h.h,"cm"
      ),
      column_title = "Marker Genes (Top 10)",
      column_names_rot = 90,
      column_names_gp = gpar(
        fontsize = fs.c
      ),
      row_names_side = "left",
      row_names_gp = gpar(
        fontsize = fs.r
      ),
      cluster_columns = F,
      cluster_rows = F
    )
    
    return(
      h.out
    )
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
