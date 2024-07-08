#### Single Volcano Plot ####

  # Plot function
  fun.vol.s <- function (
    df,
    group_fc, 
    group_p,
    name.lab,
    fc.lim,
    title
    ) {
    # Plot function
    v <- EnhancedVolcano(
      df, 
      lab = df[[name.lab]],
      title = element_blank(),
      subtitle = element_blank(), 
      caption = element_blank(), 
      x= group_fc, 
      y= group_p,
      pCutoff = 0.005, 
      FCcutoff = 0.25,
      cutoffLineType = 'twodash',
      legendLabels = c('NS','Fold Change',
                       'p-value','FC+P'), 
      legendLabSize = 12,
      labFace = 'bold', 
      col = ggsci::pal_npg("nrc")(10)[c(4,3,5,8)],
      colAlpha = 0.7,
      legendIconSize = 4,
      pointSize = 2,
      border = 'full',
      borderWidth = 1.5,
      legendPosition = 'right',
      labSize = 3,
      drawConnectors = T,
      typeConnectors = "open",
      min.segment.length = unit(1,
                                "mm")
      ) +
      
      thm.univ +
      labs(color='Key') +
      ggtitle(title) +
      theme(plot.margin = unit(c(.2,.2,.2,.2),"cm")) +
      ggplot2::coord_cartesian(xlim=c(-fc.lim,fc.lim),
                               ylim = c(0,350)) +
      ggplot2::scale_x_continuous(breaks=seq(-5,5,1))
    
    return(v)
    
    }

  
  #### Single Volcano Plot ####
  
  # Plot function
  # fun.vol.s <- function (
  #   df,
  #   group_fc, 
  #   group_p,
  #   name.lab,
  #   fc.lim,
  #   title
  # ) {
  #   # Plot function
  #   v <- EnhancedVolcano(
  #     df, 
  #     lab = df[[name.lab]],
  #     title = element_blank(),
  #     subtitle = element_blank(), 
  #     caption = element_blank(), 
  #     x= group_fc, 
  #     y= group_p,
  #     pCutoff = 0.005, 
  #     FCcutoff = 0.25,
  #     cutoffLineType = 'twodash',
  #     legendLabels = c('NS','Fold Change',
  #                      'p-value','FC+P'), 
  #     legendLabSize = 12,
  #     labFace = 'bold', 
  #     col = ggsci::pal_npg("nrc")(10)[c(4,3,5,8)],
  #     colAlpha = 0.7,
  #     legendIconSize = 4,
  #     pointSize = 2,
  #     border = 'full',
  #     borderWidth = 1.5,
  #     legendPosition = 'right',
  #     labSize = 3,
  #     drawConnectors = T,
  #     typeConnectors = "open",
  #     min.segment.length = unit(1,
  #                               "mm")
  #   ) +
  #     
  #     thm.univ +
  #     labs(color='Key') +
  #     ggtitle(title) +
  #     theme(plot.margin = unit(c(.2,.2,.2,.2),"cm")) +
  #     ggplot2::coord_cartesian(xlim=c(-fc.lim,fc.lim),
  #                              ylim = c(0,350)) +
  #     ggplot2::scale_x_continuous(breaks=seq(-5,5,1))
  #   
  #   return(v)
  #   
  # }  
  
  




