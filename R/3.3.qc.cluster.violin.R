#### General Violin Plot Function ####


#---- Plot ----

fun.qc.vio <- function(
    df,
    v1,
    var.cl
) {
  
  v <- ggplot(
    df,
    aes(
      x = .data[[var.cl]],
      y = .data[[v1]],
      fill = .data[[var.cl]]
    )
  ) +
    
    scale_fill_manual(
      name = "Cluster",
      values = col1
    ) +
    
    # Add violin plot and dotplot
    
    geom_violin(
      trim = F
    ) +
    
    geom_jitter(
      aes(
        alpha = 0.2
      ),
      shape = 16,
      size = 0.1,
      position = position_jitter(
        width = 0.4
      ),
      show.legend = F
    ) +
    
    # Add Theme
    
    thm.univ +
    thm.leg.title.y +
    labs(
      y = v1
    ) +
    theme(
      plot.margin = unit(
        c(
          0.1,
          0.1,
          0.1,
          0.1
        ),
        "cm"
      )
    )
  
  return(v)
  
}





#---- Count proportion of cells in each cluster by group and output as table ----

fun.cluster.counts <- function(
    so,
    var.clus,
    var1,
    title1
) {
  
  d <- so  
  
  cl.counts <- dplyr::count(
    dplyr::group_by(
      d@meta.data[,c(
        var.clus,
        var1
      )],
      .data[[var.clus]]
    ),
    .data[[var1]]
  )
  
  cl.counts <- dplyr::left_join(
    cl.counts,
    setNames(
      aggregate(
        cl.counts$n,
        list(
          cl.counts[[var.clus]]
        ),
        FUN = sum
      ),
      c(
        var.clus,
        "Total Cells"
      )
    ),
    by = var.clus
  ) %>%
    dplyr::mutate(
      "Proportion" = round(
        .data[["n"]]/
          .data[["Total Cells"]],
        digits = 2
      )
    )
  
  # Save as table
  
  write.table(
    cl.counts,
    paste(
      "Analysis/Integrated/",
      "umap.",
      title1,
      ".",
      var1,
      ".counts.txt",
      sep = ""
    ),
    sep = "\t",
    col.names = T,
    row.names = F
  )
  
  return(cl.counts)
  
}






