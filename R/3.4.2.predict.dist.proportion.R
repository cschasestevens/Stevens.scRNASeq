#### Cell Type Prediction Score Distribution ####

fun.dist.score <- function(x) 
{
  
  p.score.dist <- ggplot(
    x,
    aes(
      x = .data[["prediction.score.max"]]
      )
    ) +
    
    geom_density(
      color = col1[[4]],
      fill = col1[[3]]
      ) +
    
    labs(
      x = "Prediction Score",
      y = "Density",
      title = "Prediction Score Distribution"
      ) +
    
    thm.univ
  
}






#### Predicted Cell Type Proportions for Each Cluster ####

# Calculate proportion of predicted cell types within each cluster

fun.predict.prop <- function(
    x,
    md
    ) {
  
  data.pred.prop2 <- setNames(
    dplyr::count(
      x,
      .data[["seurat.clusters"]]
      ),
    c(
      "seurat.clusters",
      "Total Cells"
      )
    )
  
  data.pred.prop <- dplyr::count(
    x,
    x[,c(
      md,
      "seurat.clusters",
      "predicted.id"
    )]
    )
  
  data.pred.prop <- dplyr::left_join(
    data.pred.prop,
    data.pred.prop2,
    by = "seurat.clusters"
    )
  
  data.pred.prop[["Proportion"]] <- round(
    data.pred.prop$n/
      data.pred.prop$`Total Cells`,
    digits = 2
    )
  
  return(data.pred.prop)
  
}





# Alternate counting function

fun.predict.prop.alt <- function(
    x,
    c,
    md
) {
  
  data.pred.prop2 <- setNames(
    dplyr::count(
      x,
      .data[[c]]
    ),
    c(
      c,
      "Total Cells"
    )
  )
  
  data.pred.prop <- dplyr::count(
    x,
    x[,c(
      md,
      c
    )]
  )
  
  data.pred.prop <- dplyr::left_join(
    data.pred.prop,
    data.pred.prop2,
    by = c
  )
  
  data.pred.prop[["Proportion"]] <- round(
    data.pred.prop$n/
      data.pred.prop$`Total Cells`,
    digits = 2
  )
  
  return(data.pred.prop)
  
}

fun.prop.generic <- function(
    x,
    c,
    md
    ) {
  
  data.pred.prop2 <- setNames(
    dplyr::count(
      x,
      .data[[c]]
      ),
    c(
      c,
      "Total Cells"
      )
    )
  
  data.pred.prop <- dplyr::count(
    x,
    x[,c(
      md,
      c
      )]
    )
  
  data.pred.prop <- dplyr::left_join(
    data.pred.prop,
    data.pred.prop2,
    by = c
    )
  
  data.pred.prop[["Proportion"]] <- round(
    data.pred.prop$n/
      data.pred.prop$`Total Cells`,
    digits = 3
    )
  
  return(data.pred.prop)
  
  }









