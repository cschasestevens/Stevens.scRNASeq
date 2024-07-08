#### Function for Creating GO Enrichment Plots ####

# GO Enrichment Plot

fun.GO.plot <- function(df,
                        p.title) {
  
  # Input
  
  go.in <- slice_min(
    df,
    order_by = .data[["P-Value"]],
    n = 10
    )
  
  go.in[["FDR"]] <- ifelse(
    is.na(
      go.in[["FDR"]]
      ),
    1e-30,
    go.in[["FDR"]]
    )
  
  go.in[["sig.rat"]] <- go.in[["Significant"]]/
    go.in[["Annotated"]]
  
  
  # Plot
  
  cr.plot <- ggplot(
    go.in,
    aes(
      x = Term,
      y = -log10(FDR),
      fill = sig.rat)
    ) +
    
    ## Plot points
    geom_point(
      shape = 21,
      col = "black",
      aes(
        size = .data[["Annotated"]]
        ),
      alpha = 0.7
      ) +
    
    ## Significance line
    geom_hline(
      yintercept = 1.3,
      linetype = "dashed"
      ) +
    
    ## Data scale (only if needed)
    scale_y_continuous(
      limits = c(
        0,10
        ),
      breaks = c(
        0,2,4,
        6,8,10)
      ) +
    
    ## Theme
    thm.mult +
    theme(
      axis.text.x = element_text(
        size = 12
        ),
      plot.margin = unit(
        c(
          0.5,0.5,
          0.5,7),
        "cm"
        ),
      legend.title = element_text(
        size = 16
        ),
      legend.text = element_text(
        size = 12
        ),
      legend.key.size = unit(
        0.4,
        "cm")
      ) +
    
    ## Axis labels
    xlab("") +
    ylab("-log10(p)") +
    ggtitle(p.title) +
    
    ## Gradient and size scaling
    scale_fill_gradientn(
      name = expression(
        bold(
          paste(
            " Proportion \n Sig. in Term"
            )
          )
        ),
      colors = col2a[c(1,3,6,9,12)],
      breaks = c(
        0,.25,.5,
        .75,1
        ),
      labels = c(
        "0","0.25","0.5",
        "0.75","1"
        )
      ) +
    scale_size_area(
      name = "Genes in Term",
      max_size = 16
      )
  
  
  # Output
  
  return(cr.plot)
  
}




