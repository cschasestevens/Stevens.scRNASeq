#' Volcano Plot
#'
#' Generates a volcano plot from a DEG results object for a chosen cell type and comparison.
#'
#' @param l.deg A list of DGEA results returned by sc.DGEA().
#' @param ct A pattern provided as a character string for matching a specific cell type.
#' @param comp.name comparison name specified by the user and returned by sc.DGEA, provided as a character string.
#' @param p.cut Numeric value for p-value cutoff to indicate significance.
#' @param f.cut Numeric value for fold change cutoff to indicate a high effect size.
#' @param f.lim Numeric value for fold change limits on the x-axis (given in multiples of 6).
#' @return A volcano plot for the chosen cell type and treatment comparison.
#' @examples
#'
#' p.vol <- sc.volcano(dgea.output,"1.Secretory","KO vs. Control",0.005,0.25,6)
#'
#' @export
sc.volcano <- function (
    l.deg,
    ct,
    comp.name,
    p.cut,
    f.cut,
    f.lim
    ) {
    # Subset Input data
    ld <- l.deg[[1]]
    ld <- ld[ld[["CellType"]] == ct & ld[["Comparison"]] == comp.name,]

    # Plot function
    v <- EnhancedVolcano::EnhancedVolcano(
      ld,
      lab = ld[["GENE"]],
      title = ggplot2::element_blank(),
      subtitle = ggplot2::element_blank(),
      caption = ggplot2::element_blank(),
      x= "log2FC",
      y= "H.qval",
      pCutoff = p.cut,
      FCcutoff = f.cut,
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
      min.segment.length = ggplot2::unit(1,
                                "mm")
      ) +
      sc.theme1() +
      ggplot2::labs(color='Key') +
      ggplot2::ggtitle(paste(comp.name,ct,sep = " ")) +
      ggplot2::theme(plot.margin = ggplot2::unit(c(.2,.2,.2,.2),"cm")) +
      ggplot2::coord_cartesian(xlim=c(-f.lim,f.lim),
                               ylim = c(0,350)) +
      ggplot2::scale_x_continuous(breaks=seq(-f.lim,f.lim,(f.lim/6)))

    return(v)

    }




