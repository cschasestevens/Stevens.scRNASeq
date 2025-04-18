#' Universal Color Palette
#'
#' Combines the npg, aaas, and lancet ggsci palettes for use with datasets
#' containing up to 36 groups.
#'
#' @return Vector of colors to replace default discrete color scale.
#' @examples
#'
#'  # col_univ()
#'
#' @export
col_univ <- function() {
  c(
    ggsci::pal_npg("nrc")(10),
    ggsci::pal_aaas("default")(10),
    ggsci::pal_lancet("lanonc")(9),
    ggsci::pal_frontiers("default")(7),
    ggsci::pal_nejm("default")(8),
    ggsci::pal_jco("default")(10)
  )
}

#' Gradient Color Palette
#'
#' Returns 12 colors from the viridis color palette.
#'
#' @return Vector of colors to replace default gradient color scale.
#' @examples
#'
#' col_grad()
#'
#' @export
col_grad <- function() {
  viridis::viridis(12)
}


#' Generic Plot Theme
#'
#' General plotting theme.
#'
#' @return ggplot2 theme parameters to replace default plot theme.
#' @examples
#'
#' # sc_theme1()
#'
#' @export
sc_theme1 <- function() {
  thm_gen <- ggplot2::theme(
    # Plot Title
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    # Panel
    panel.border = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(colour = "grey85"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    # Axes
    axis.ticks.y = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      face = "bold",
      size = 14,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 14
    ),
    axis.title.x = ggplot2::element_text(
      face = "bold",
      size = 14
    ),
    axis.title.y = ggplot2::element_text(
      face = "bold",
      size = 14
    ),
    # Strip
    strip.background = ggplot2::element_rect(
      fill = "slategray2"
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    # Margins
    plot.margin = ggplot2::unit(
      c(0.5, 0.25, 0.5, 0.25),
      "cm"
    )
  )

  thm_leg_main <- ggplot2::theme(
    legend.title = ggplot2::element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(size = 10),
    legend.key.size = ggplot2::unit(0.2, "cm"),
    legend.key = ggplot2::element_blank(),
    legend.position.inside = c(0.95, 0.95)
  )
  thm_mult <- thm_gen +
    thm_leg_main
  return(thm_mult)
}
