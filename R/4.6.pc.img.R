#' Phenocycler Cell Map
#'
#' Generates phenocycler image scatterplot from Seurat object; Must have
#' X and Y coordinates available in the input dataframe.
#'
#' @param so Input Seurat object containing coordinates and meta data.
#' @param type Scatterplot type ("int" or "var").
#' @param md_var If type is "int," specify channel name.
#' @param col1 Gradient color column to use.
#' @return Cell map scatter plot for selected variable.
#' @examples
#'
#' # psc <- pc_scatter(
#' #   so = d1,
#' #   md_var = "cluster"
#' # )
#'
#' @export
pc_scatter <- function(
  so,
  type = "var",
  md_var,
  col1 = col_grad(scm = 4)
) {
  # Extract data frame
  d <- so
  d2 <- SeuratObject::FetchData(
    d,
    vars = c(
      "X",
      "Y",
      md_var
    )
  )
  if(type == "var") { # nolint
    p1 <- ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x = X, # nolint
        y = Y # nolint
      )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(
          color = .data[[md_var]] # nolint
        ),
        shape = 16,
        size = 0.5,
        alpha = 0.7
      ) +
      pc_theme_img() + # nolint
      ggplot2::labs(fill = md_var) +
      ggplot2::scale_y_reverse() +
      ggplot2::scale_color_manual(
        values = col1 # nolint
      )
  }
  if(type == "int") { # nolint
    p1 <- ggplot2::ggplot(
      d2,
      ggplot2::aes(
        x = X, # nolint
        y = Y # nolint
      )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(
          color = .data[[md_var]] # nolint
        ),
        shape = 16,
        size = 0.5
      ) +
      pc_theme_img() + # nolint
      ggplot2::labs(fill = md_var) +
      ggplot2::scale_y_reverse() +
      ggplot2::scale_color_gradientn(
        colors = col1, # nolint
        limits = c(
          0,
          quantile(
            d2[[md_var]],
            0.99
          )
        ),
        na.value = col1[[1]]
      )
  }
  return(p1)
}
