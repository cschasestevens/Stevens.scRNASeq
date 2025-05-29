#' Phenocycler QC
#'
#' Runs quality control of segmented phenocycler data for all samples
#'
#' @param qc_type Specific QC steps to perform (
#' select from "all" or "sum_only")
#' @param ld A list of phenocycler segmentation results containing summary data
#' and segmented cells.
#' @param col_nums1 Vector of column numbers to include from segmentation
#' summary data.
#' @param md_var Vector for renaming summary data column names.
#' @param md_num Number of metadata columns in exported cell segmentation
#' results. Only used if qc_type is "all".
#' @param samp_id Sample ID column name.
#' @param g_col (optional) Add a column providing grouping information.
#' @param g_col2 (optional) Add a column providing grouping information.
#' @param ch_list List of channels exported with segmentation results. Note
#' that the channel list must match the order of channels included within
#' the results.
#' @param ch_nuc Nuclei channel name.
#' @param col_nums2 If qc_type is "all", indicate the column number to
#' use as the sample ID as well as the X and Y coordinates
#' for formatting the segmentation results table.
#' @param norm_data Perform data normalization if TRUE.
#' @param loess_norm Normalize signal intensities by performing LOESS
#' prior to z-score scaling.
#' @return List containing QC summary statistics and plots.
#' @examples
#'
#' # pheno_data <- pc_qc(
#' #   ld = ld1,
#' #   col_nums1 = c(4, 10, 12, 13),
#' #   md_var = c("Area.um2", "Cells.raw", "Description"),
#' #   samp_id = "Name",
#' #   g_col = "Name",
#' #   g_col2 = "Parent",
#' #   ch_list = s3,
#' #   col_nums2 = c(6, 8, 9)
#' # )
#'
#' @export
pc_qc <- function(
  qc_type = "all",
  ld,
  col_nums1,
  md_var,
  md_num = 9,
  samp_id = "Code",
  g_col = NULL,
  g_col2 = NULL,
  ch_list,
  col_nums2 = NULL,
  ch_nuc = "DAPI",
  norm_data = FALSE,
  loess_norm = FALSE
) {
  # Define parameters
  d1 <- ld
  cn1 <- col_nums1
  cn2 <- md_var
  cn1b <- col_nums2
  sid <- samp_id
  chl <- ch_list
  mdn <- md_num
  chn <- ch_nuc
  # Perform QC
  if(qc_type == "all") { # nolint
    # Summary statistics
    ## Count section dimensions and detected cells per section
    qc1 <- setNames(dplyr::bind_rows(lapply(
      seq.int(1, length(d1), 1),
      function(x) {
        data.frame(
          "Slide" = x,
          d1[[x]][[1]][, sid],
          d1[[x]][[1]][, cn1]
        )
      }
    )), c("Slide", "Code", cn2))
    if(is.null(g_col) == FALSE) { # nolint
      qc1 <- dplyr::select(
        dplyr::mutate(
          qc1,
          "Group" = ifelse(grepl("KK", qc1[["Code"]]), "CF", "Norm")
        ),
        c("Slide", "Code", "Group"), everything() # nolint
      )
    }
    # Intensity-based filtering (by nuclei and total signal)
    ## Format segmentation results (Cellpose results contain summary stats!)
    ### Formatting for exported Cellpose measurements
    qc2 <- dplyr::bind_rows(lapply(
      seq.int(1, length(d1), 1),
      function(x) {
        setNames(
          data.frame(
            "Slide" = rep(x, nrow(d1[[x]][[2]])),
            d1[[x]][[2]]
          ),
          c(
            "Slide",
            names(d1[[x]][[2]][1:mdn]),
            unlist(
              lapply(
                seq.int(1, nrow(chl), 1),
                function(x) {
                  st1 <- c("mean", "median", "min", "max", "sd")
                  st2 <- unlist(
                    lapply(
                      seq.int(1, length(st1), 1),
                      function(y) {
                        paste(
                          chl[x, 1],
                          st1[[y]],
                          sep = "_"
                        )
                      }
                    )
                  )
                  return(st2) # nolint
                }
              )
            )
          )
        )
      }
    ))
    qc2 <- setNames(data.frame(
      "ID" = seq.int(1, nrow(qc2), 1),
      qc2
    ), c("ID", names(qc2)))
    ## remove all columns except mean intensity per cell
    qc2 <- qc2[
      ,
      c(
        names(qc2[, 1:(mdn + 2)]),
        names(qc2[, grepl("mean", names(qc2))])
      )
    ]
    qc2 <- setNames(
      reshape2::melt(
        qc2,
        id.vars = 1:(mdn + 2),
        variable.name = "column",
        value.name = "value"
      )[, c(
        1, 2,
        (cn1b[[1]] + 2),
        (cn1b[[2]] + 2),
        (cn1b[[3]] + 2),
        (mdn + 3),
        (mdn + 4)
      )],
      c("ID", "Slide", "Code", "X", "Y", "column", "value")
    )
    if(is.null(g_col2) == FALSE) { # nolint
      qc2 <- dplyr::select(
        dplyr::mutate(
          qc2,
          "Group" = ifelse(grepl("KK", qc2[["Code"]]), "CF", "Norm"),
          "Channel" = gsub("\\_.*", "", qc2[["column"]])
        ),
        c("ID", "Slide", "Code", "Group", "X", "Y", "Channel", "value")
      )
    }
    if(is.null(g_col2) == TRUE) { # nolint
      qc2 <- dplyr::select(
        dplyr::mutate(
          qc2,
          "Channel" = gsub("\\_.*", "", qc2[["column"]])
        ),
        c("ID", "Slide", "Code", "X", "Y", "Channel", "value")
      )
    }
    # Calculate sum intensity per cell
    if(Sys.info()[["sysname"]] != "Windows") { # nolint
      qc_sum <- dplyr::bind_rows(
        parallel::mclapply(
          mc.cores = 2,
          seq.int(1, length(unique(qc2[["Code"]])), 1),
          function(x) {
            md1 <- unique(qc2[["Code"]])
            d1 <- qc2[qc2[["Code"]] == md1[[x]], ]
            d1 <- setNames(
              aggregate(
                d1[["value"]],
                list(d1[["ID"]]),
                function(y) sum(y)
              ),
              c("ID", "sum.int.raw")
            )
            d1 <- dplyr::select(
              dplyr::mutate(
                d1,
                "Code" = rep(md1[[x]], nrow(d1))
              ),
              "Code",
              everything() # nolint
            )
            return(d1) # nolint
          }
        )
      )
    }
    if(Sys.info()[["sysname"]] == "Windows") { # nolint
      qc_sum <- dplyr::bind_rows(
        lapply(
          seq.int(1, length(unique(qc2[["Code"]])), 1),
          function(x) {
            md1 <- unique(qc2[["Code"]])
            d1 <- qc2[qc2[["Code"]] == md1[[x]], ]
            d1 <- setNames(
              aggregate(
                d1[["value"]],
                list(d1[["ID"]]),
                function(y) sum(y)
              ),
              c("ID", "sum.int.raw")
            )
            d1 <- dplyr::select(
              dplyr::mutate(
                d1,
                "Code" = rep(md1[[x]], nrow(d1))
              ),
              "Code",
              everything() # nolint
            )
            return(d1) # nolint
          }
        )
      )
    }
    ### Flag cells based on Total intensity
    filt1 <- qc_sum[
      qc_sum[["sum.int.raw"]] < quantile(qc_sum[["sum.int.raw"]], 0.01) |
        qc_sum[["sum.int.raw"]] > quantile(qc_sum[["sum.int.raw"]], 0.99),
      "ID"
    ]
    ### Flag cells based on DAPI intensity
    filt2 <- qc2[qc2[["Channel"]] == chn, ]
    filt2 <- filt2[
      filt2[["value"]] < quantile(filt2[["value"]], 0.01) |
        filt2[["value"]] > quantile(filt2[["value"]], 0.99),
      "ID"
    ]
    ### Remove flagged cells from individual and aggregated data
    qc2_filt <- qc2[
      qc2[["ID"]] %in% filt1 == FALSE &
        qc2[["ID"]] %in% filt2 == FALSE,
    ]
    qc_sum_filt <- qc_sum[
      qc_sum[["ID"]] %in% filt1 == FALSE &
        qc_sum[["ID"]] %in% filt2 == FALSE,
    ]
    fprp <- round((1 - (nrow(qc2_filt) / nrow(qc2))), digits = 2) * 100
    print(paste(fprp, "% ", "of cells removed after filtering...", sep = ""))
    qc_sum_cnt <- setNames(dplyr::count(
      qc_sum_filt,
      .data[["Code"]] # nolint
    ), c("Code", "Cells.filt"))
    qc1 <- dplyr::left_join(
      qc1,
      qc_sum_cnt,
      by = "Code"
    )
    ### Merge coordinates with sum intensity data frame
    qc_sum_filt <- dplyr::left_join(
      qc_sum_filt,
      qc2_filt[!duplicated(qc2_filt[["ID"]]), c("ID", "X", "Y")],
      by = "ID"
    )
    # Normalizations
    if(norm_data == TRUE) { # nolint
      print("Normalizing data...")
      ## Z-score only
      #### formula (x - mean)/sd # nolint
      qc2_filt[["value.z"]] <- (
        qc2_filt[["value"]] - mean(qc2_filt[["value"]])
      ) /
        sd(qc2_filt[["value"]])
      ## LOESS + Z-score (Experimental) # nolint
      if(loess_norm == TRUE) { # nolint
      }
      # Append z-normalized intensities to data frame
      if(Sys.info()[["sysname"]] != "Windows") { # nolint
        qc_sum_filt2 <- dplyr::bind_rows(
          parallel::mclapply(
            mc.cores = 2,
            seq.int(1, length(unique(qc2_filt[["Code"]])), 1),
            function(x) {
              md1 <- unique(qc2_filt[["Code"]])
              d1 <- qc2_filt[qc2_filt[["Code"]] == md1[[x]], ]
              d1 <- setNames(
                aggregate(
                  d1[["value.z"]],
                  list(d1[["ID"]]),
                  function(y) sum(y)
                ),
                c("ID", "sum.int.znorm")
              )
              d1 <- dplyr::select(
                dplyr::mutate(
                  d1,
                  "Code" = rep(md1[[x]], nrow(d1))
                ),
                "Code",
                everything() # nolint
              )
              return(d1) # nolint
            }
          )
        )
      }
      if(Sys.info()[["sysname"]] == "Windows") { # nolint
        qc_sum_filt2 <- dplyr::bind_rows(
          lapply(
            seq.int(1, length(unique(qc2_filt[["Code"]])), 1),
            function(x) {
              md1 <- unique(qc2_filt[["Code"]])
              d1 <- qc2_filt[qc2_filt[["Code"]] == md1[[x]], ]
              d1 <- setNames(
                aggregate(
                  d1[["value.z"]],
                  list(d1[["ID"]]),
                  function(y) sum(y)
                ),
                c("ID", "sum.int.znorm")
              )
              d1 <- dplyr::select(
                dplyr::mutate(
                  d1,
                  "Code" = rep(md1[[x]], nrow(d1))
                ),
                "Code",
                everything() # nolint
              )
              return(d1) # nolint
            }
          )
        )
      }
      qc_sum_filt <- dplyr::select(dplyr::left_join(
        qc_sum_filt,
        qc_sum_filt2[, c("ID", "sum.int.znorm")],
        by = "ID"
      ), c("Code", "ID", "sum.int.raw", "sum.int.znorm"), everything()) # nolint
      # remove redundant objects
      remove(qc_sum_filt2, qc2, qc_sum, qc_sum_cnt, d1, ld)
      gc(reset = TRUE)
    }
    if(norm_data == FALSE) { # nolint
      # Return summary stats
      qc_norm <- data.frame(
        "mean" = c(
          round(mean(qc2_filt[["value"]]), digits = 2)
        ),
        "std.dev" = c(
          round(sd(qc2_filt[["value"]]), digits = 2)
        )
      )
      print("Returning raw values...")
      qc_norm
      # Remove objects
      remove(qc_sum_filt2, qc2, qc_sum, qc_sum_cnt, d1, ld)
      gc(reset = TRUE)
    }
  }
  return(
    list(
      "summary" = qc1,
      "data_ind" = qc2_filt,
      "data_total" = qc_sum_filt
    )
  )
}

#' Phenocycler QC Boxplot
#'
#' Generates boxplot of signal intensity
#'
#' @param df Input data frame containing signal intensities.
#' @param ch_list Channel list for specifying plot order.
#' @param ch_col Channel column for spliting individual boxplots. Used to
#' specify the x-axis and can alternatively plot a different variable if
#' individual channel intensities are not available
#' @param int_col Intensity column to use.
#' @return Boxplot of individual channel or summary of signal intensities.
#' @examples
#'
#' # pheno_qc_box <- pc_qc_box(
#' #   df = ld1,
#' #   ch_list = s3[[1]]
#' #   int_col = "value.z"
#' # )
#'
#' @export
pc_qc_box <- function(
  df,
  ch_list = NULL,
  ch_col = "Channel",
  int_col = NULL
) {
  d1 <- df
  if(is.null(ch_list) == FALSE) { # nolint
    p1 <- ggplot2::ggplot(
      d1,
      ggplot2::aes(
        x = factor(.data[[ch_col]], levels = ch_list), # nolint
        y = .data[[int_col]]
      )
    ) +
      ggplot2::geom_boxplot() +
      ggplot2::coord_flip() +
      # ggplot2::geom_hline(
      #   color = "black", # nolint
      #   linetype = "dashed", # nolint
      #   yintercept = quantile(qc2_filt[["value"]], 0.01) # nolint
      # ) +
      ggplot2::labs(
        y = "Intensity",
        x = "Channel"
      ) +
      sc_theme1() # nolint
  }
  if(is.null(ch_list) == TRUE) { # nolint
    p1 <- ggplot2::ggplot(
      d1,
      ggplot2::aes(
        x = as.factor(.data[[ch_col]]), # nolint
        y = .data[[int_col]]
      )
    ) +
      ggplot2::geom_boxplot() +
      ggplot2::coord_flip() +
      # ggplot2::geom_hline(
      #   color = "black", # nolint
      #   linetype = "dashed", # nolint
      #   yintercept = quantile(qc2_filt[["value"]], 0.01) # nolint
      # ) +
      ggplot2::labs(
        y = "Intensity",
        x = ch_col
      ) +
      sc_theme1() # nolint
  }
  return(p1)
}

#' Phenocycler QC Scatterplot
#'
#' Generates scatterplot of signal intensity; Must have
#' X and Y coordinates available in the input dataframe.
#'
#' @param df Input data frame containing signal intensities.
#' @param ch_col Channel column; Used to
#' specify the x-axis and can alternatively plot a different variable if
#' individual channel intensities are not available.
#' @param ch_id If plotting an individual channel, which channel should
#' be plotted?
#' @param int_col Intensity column to use.
#' @param id_col ID column containing sample information.
#' @param id_samp Specific sample ID to plot.
#' @return Boxplot of individual channel or summary of signal intensities.
#' @examples
#'
#' # pheno_qc_scatter <- pc_qc_scatter(
#' #   df = qc2_filt,
#' #   ch_id = "DAPI",
#' #   int_col = "value.z",
#' #   id_samp = "H10687"
#' # )
#'
#' @export
pc_qc_scatter <- function(
  df,
  ch_col = "Channel",
  ch_id = NULL,
  int_col,
  id_col = "Code",
  id_samp
) {
  d1 <- df
  if(is.null(ch_id) == TRUE) { # nolint
    p1 <- ggplot2::ggplot(
      d1[d1[[id_col]] == id_samp, ],
      ggplot2::aes(
        x = X, # nolint
        y = Y # nolint
      )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(
          color = .data[[int_col]] # nolint
        ),
        shape = 16,
        size = 0.5
      ) +
      pc_theme_img() + # nolint
      ggplot2::labs(fill = "Intensity") +
      ggplot2::scale_y_reverse() +
      ggplot2::scale_color_gradientn(
        colors = col_grad(scm = 4), # nolint
        limits = c(
          0,
          quantile(
            d2[[int_col]],
            0.99
          )
        ),
        na.value = col_grad(scm = 4)[[1]]
      )
  }
  if(is.null(ch_id) == FALSE) { # nolint
    d2 <- d1[d1[[ch_col]] == ch_id, ]
    p1 <- ggplot2::ggplot(
      d2[d2[[id_col]] == id_samp, ],
      ggplot2::aes(
        x = X, # nolint
        y = Y # nolint
      )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(
          color = .data[[int_col]] # nolint
        ),
        shape = 16,
        size = 0.5
      ) +
      pc_theme_img() + # nolint
      ggplot2::labs(fill = paste("Intensity: ", ch_id, sep = "")) +
      ggplot2::scale_y_reverse() +
      ggplot2::scale_color_gradientn(
        colors = col_grad(scm = 4), # nolint
        limits = c(
          0,
          quantile(
            d2[[int_col]],
            0.99
          )
        ),
        na.value = col_grad(scm = 4)[[1]]
      )
  }
  return(p1)
}
