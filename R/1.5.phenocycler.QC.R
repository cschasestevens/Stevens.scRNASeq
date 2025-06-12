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
#' @return List containing QC summary statistics and plots.
#' @examples
#'
#' # pheno_qc <- pc_qc(
#' #   ld = ld1,
#' #   samp_id = "Name",
#' #   g_col = "Name",
#' #   g_col2 = "Parent",
#' #   ch_list = s3
#' # )
#'
#' @export
pc_qc <- function(
  qc_type = "all",
  ld,
  col_nums1 = c(10, 12, 13),
  md_var = c("Area.um2", "Cells.raw", "Description"),
  md_num = 9,
  samp_id = "Code",
  g_col = NULL,
  g_col2 = NULL,
  ch_list,
  col_nums2 = c(6, 8, 9),
  ch_nuc = "DAPI"
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
  ### Flag cells based on Total intensity (remove bottom 5% and top 5%)
  # Determine intensity quantiles per slide
  thr1 <- data.frame(
    "Percentile" = paste(
      c(
        0, 5, 10, 20, 30, 40,
        50, 60, 70, 80, 90,
        95, 100
      ),
      "%",
      sep = ""
    ),
    dplyr::bind_cols(setNames(lapply(
      unique(qc_sum[["Code"]]),
      function(x) {
        qcd <- qc_sum[qc_sum[["Code"]] == x, ]
        qcd <- round(quantile(
          qcd[["sum.int.raw"]],
          c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
        ), digits = 0)
        qcd <- setNames(
          data.frame("val" = qcd),
          c(x)
        )
        return(qcd) # nolint
      }
    ), unique(qc_sum[["Code"]])))
  )
  print("Quantiles for raw total intensities are as follows:")
  print(thr1)
  head(qc2)
  thr2 <- data.frame(
    "Percentile" = paste(
      c(
        0, 5, 10, 20, 30, 40,
        50, 60, 70, 80, 90,
        95, 100
      ),
      "%",
      sep = ""
    ),
    dplyr::bind_cols(setNames(lapply(
      unique(qc_sum[["Code"]]),
      function(x) {
        qcd <- qc2[qc2[["Code"]] == x & qc2[["Channel"]] == chn, ]
        qcd <- round(quantile(
          qcd[["value"]],
          c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
        ), digits = 0)
        qcd <- setNames(
          data.frame("val" = qcd),
          c(x)
        )
        return(qcd) # nolint
      }
    ), unique(qc2[["Code"]])))
  )
  print("Quantiles for raw DAPI intensities are as follows:")
  print(thr2)
  # Visualize raw, unfiltered intensities
  qcp1 <- ggplot2::ggplot(
    data = qc_sum,
    ggplot2::aes(
      x = ID, # nolint
      y = sum.int.raw, # nolint
      color = Code # nolint
    )
  ) +
    ggplot2::geom_point(
      shape = 16,
      size = 0.5
    ) +
    ggplot2::scale_color_manual(values = col_univ()) +
    sc_theme1()
  qcp2 <- ggplot2::ggplot(
    data = qc2[qc2[["Channel"]] == chn, ],
    ggplot2::aes(
      x = ID, # nolint
      y = value, # nolint
      color = Code # nolint
    )
  ) +
    ggplot2::geom_point(
      shape = 16,
      size = 0.5
    ) +
    ggplot2::scale_color_manual(values = col_univ()) +
    sc_theme1()
  qcp3 <- ggpubr::ggarrange(
    qcp1, qcp2, labels = c("Total Raw", "DAPI Raw"),
    nrow = 1,
    ncol = 2
  )
  print(qcp3)
  # Flag based on total and DAPI intensity
  filt1 <- unlist(lapply(
    names(thr1[2:ncol(thr1)]),
    function(x) {
      # Total
      qcd <- qc_sum[qc_sum[["Code"]] == x, ]
      qcd1 <- nrow(qcd)
      print(paste(qcd1, " cells present in ", x, ".", sep = ""))
      qcd <- qcd[
        qcd[["sum.int.raw"]] < thr1[thr1[["Percentile"]] == "5%", x] |
          qcd[["sum.int.raw"]] > thr1[thr1[["Percentile"]] == "95%", x], # nolint
      ]
      qcd2 <- nrow(qcd)
      print(paste(
        qcd2, " cells in ", x,
        " flagged based on total intensity filter.", sep = ""
      ))
      # DAPI
      qcd3 <- qc2[qc2[["Code"]] == x & qc2[["Channel"]] == chn, ]
      qcd3 <- qcd3[
        qcd3[["value"]] < thr2[thr2[["Percentile"]] == "5%", x] |
          qcd3[["value"]] > thr2[thr2[["Percentile"]] == "95%", x], # nolint
      ]
      qcd5 <- nrow(qcd3[qcd[["ID"]] %in% qcd3[["ID"]] == FALSE, ])
      print(paste(
        qcd5, " additional cells in ", x,
        " flagged based on DAPI intensity filter.", sep = ""
      ))
      qcd <- unique(c(qcd[["ID"]], qcd3[["ID"]]))
      print(paste(
        (round(((length(qcd) / qcd1)), digits = 2) * 100),
        "% of cells from ", x,
        " will be removed after filtering.", sep = ""
      ))
      return(qcd) # nolint
    }
  ))
  filt1[1:20]
  ### Remove flagged cells from individual and aggregated data
  qc2_filt <- qc2[qc2[["ID"]] %in% filt1 == FALSE, ]
  qc_sum_filt <- qc_sum[qc_sum[["ID"]] %in% filt1 == FALSE, ]
  fprp <- round((1 - (nrow(qc2_filt) / nrow(qc2))), digits = 2) * 100
  print(paste(
    fprp, "% ", "of all cells removed after filtering...", sep = ""
  ))
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
  # Return summary stats
  qc_norm <- data.frame(
    "mean" = c(
      round(mean(qc2_filt[["value"]]), digits = 2)
    ),
    "std.dev" = c(
      round(sd(qc2_filt[["value"]]), digits = 2)
    )
  )
  print("Returning raw mean and standard deviation...")
  print(qc_norm)
  # Remove objects
  remove(qc2, qc_sum, qc_sum_cnt, d1, ld)
  gc(reset = TRUE)
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

#' Phenocycler Channel QC
#'
#' Generates scatterplots of segmented cell signal intensity; Must have
#' X and Y coordinates available in the input dataframe.
#'
#' @param type Either "img" or "sum."
#' @param so Input Seurat object.
#' @param ch_id Channel name.
#' @param samp_id Sample ID column.
#' @param col1 Color scheme to use.
#' @return Panel containing cell map and axis intensity scatter plots
#' for a specified channel.
#' @examples
#'
#' # d <- pc_qc_ch(
#' #   ld_ser[[1]],
#' #   ch_id = "DAPI"
#' # )
#'
#' @export
pc_qc_ch <- function(
  type1 = "img",
  so,
  samp_id = "Code",
  ch_id,
  col1 = col_grad(scm = 4),
  loess_smooth = FALSE
) {
  # Extract data
  d <- so
  sid <- samp_id
  chid <- ch_id
  cols <- col1
  d2 <- SeuratObject::FetchData(
    d,
    vars = c(
      "X",
      "Y",
      sid,
      chid
    )
  )
  ## LOESS overlays (uses 5% of pixel number for smoothing)
  if(type1 == "img") { # nolint
    if(loess_smooth == TRUE) { # nolint
      d2x_fit <- as.data.frame(
        lowess(
          x = d2[["X"]],
          y = d2[[chid]],
          # use 10% of the average pixel number per sample for span
          f = 0.05,
          iter = 1,
          delta = 0
        )
      )
      d2y_fit <- as.data.frame(
        lowess(
          x = d2[["Y"]],
          y = d2[[chid]],
          # use 10% of the average pixel number per sample for span
          f = 0.05,
          iter = 1,
          delta = 0
        )
      )
      # Plot
      ## Cell map
      p1 <- ggplot2::ggplot(
        d2,
        ggplot2::aes(
          x = X, # nolint
          y = Y # nolint
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(
            color = .data[[chid]] # nolint
          ),
          shape = 16,
          size = 0.5
        ) +
        pc_theme_img() + # nolint
        ggplot2::labs(fill = chid) +
        ggplot2::scale_y_reverse() +
        ggplot2::scale_color_gradientn(
          colors = cols, # nolint
          limits = c(
            0,
            quantile(
              d2[[chid]],
              0.99
            )
          ),
          na.value = cols[[1]]
        )
      p2 <- ggplot2::ggplot(
        d2,
        ggplot2::aes(
          x = X, # nolint
          y = .data[[chid]] # nolint
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(
            color = .data[[chid]] # nolint
          ),
          shape = 16,
          size = 0.5
        ) +
        ggplot2::geom_line(
          data = d2x_fit,
          ggplot2::aes(
            x = x, # nolint
            y = y  # nolint
          )
        ) +
        sc_theme1() + # nolint
        ggplot2::labs(fill = chid, x = "X", y = paste(chid, "Intensity")) +
        ggplot2::scale_color_gradientn(
          colors = cols, # nolint
          limits = c(
            0,
            max(d2[[chid]])
          ),
          na.value = cols[[1]]
        )
      p3 <- ggplot2::ggplot(
        d2,
        ggplot2::aes(
          x = Y, # nolint
          y = .data[[chid]] # nolint
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(
            color = .data[[chid]] # nolint
          ),
          shape = 16,
          size = 0.5
        ) +
        ggplot2::geom_line(
          data = d2y_fit,
          ggplot2::aes(
            x = x, # nolint
            y = y  # nolint
          )
        ) +
        sc_theme1() + # nolint
        ggplot2::labs(fill = chid, x = "Y", y = paste(chid, "Intensity")) +
        ggplot2::scale_color_gradientn(
          colors = cols, # nolint
          limits = c(
            0,
            max(d2[[chid]])
          ),
          na.value = cols[[1]]
        )
    }
    if(loess_smooth == FALSE) { # nolint
      # Plot
      ## Cell map
      p1 <- ggplot2::ggplot(
        d2,
        ggplot2::aes(
          x = X, # nolint
          y = Y # nolint
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(
            color = .data[[chid]] # nolint
          ),
          shape = 16,
          size = 0.5
        ) +
        pc_theme_img() + # nolint
        ggplot2::labs(fill = chid) +
        ggplot2::scale_y_reverse() +
        ggplot2::scale_color_gradientn(
          colors = cols, # nolint
          limits = c(
            0,
            quantile(
              d2[[chid]],
              0.99
            )
          ),
          na.value = cols[[1]]
        )
      p2 <- ggplot2::ggplot(
        d2,
        ggplot2::aes(
          x = X, # nolint
          y = .data[[chid]] # nolint
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(
            color = .data[[chid]] # nolint
          ),
          shape = 16,
          size = 0.5
        ) +
        ggplot2::geom_smooth() +
        sc_theme1() + # nolint
        ggplot2::labs(fill = chid, x = "X", y = paste(chid, "Intensity")) +
        ggplot2::scale_color_gradientn(
          colors = cols, # nolint
          limits = c(
            0,
            max(d2[[chid]])
          ),
          na.value = cols[[1]]
        )
      p3 <- ggplot2::ggplot(
        d2,
        ggplot2::aes(
          x = Y, # nolint
          y = .data[[chid]] # nolint
        )
      ) +
        ggplot2::geom_point(
          ggplot2::aes(
            color = .data[[chid]] # nolint
          ),
          shape = 16,
          size = 0.5
        ) +
        ggplot2::geom_smooth() +
        sc_theme1() + # nolint
        ggplot2::labs(fill = chid, x = "Y", y = paste(chid, "Intensity")) +
        ggplot2::scale_color_gradientn(
          colors = cols, # nolint
          limits = c(
            0,
            max(d2[[chid]])
          ),
          na.value = cols[[1]]
        )
    }
    gc(reset = TRUE)
    d2 <- ggpubr::ggarrange(
      p1, p2, p3, common.legend = TRUE, nrow = 1, ncol = 3
    )
  }
  if(type1 == "sum") { # nolint
    gc(reset = TRUE)
    ## Channel mean, median, sd, min, and max
    d2 <- data.frame(
      "Code" = unique(d2[[sid]]),
      "Channel" = chid,
      "mean" = round(mean(d2[[chid]]), digits = 2),
      "median" = round(median(d2[[chid]]), digits = 2),
      "sd" = round(sd(d2[[chid]]), digits = 2),
      "min" = round(min(d2[[chid]]), digits = 2),
      "max" = round(max(d2[[chid]]), digits = 2)
    )
  }
  return(d2) # nolint
}

#' Phenocycler Normalization
#'
#' Z-score or sparse implementation of LOESS normalization for correcting
#' gradient artifacts in phenocycler data. Uses method described
#' by Stevens et. al. 2025 (doi: 10.1038/s41467-025-58135-4).
#'
#' @param so Input Seurat object.
#' @param mdn Number of metadata columns in the input Seurat object.
#' @param mtd Normalization method to use (either "zsc" or "slo").
#' @param list_ch Vector names of channels to normalize.
#' @param span1 LOESS span to use (expressed as a proportion of the total cells
#' present in the dataset). Use a proportion of 0.1 or lower for highly
#' heterogenous tissues such as lung.
#' @param mcc Number of cores to use for sLOESS calculation.
#' @param its Number of iterations for fine tuning sLOESS model.
#' @return A normalized Seurat Object for subsequent clustering
#' and dimension reduction.
#' @examples
#'
#' # d <- pc_norm(
#' #   so = ld_ser[[1]],
#' #   mdn = 9,
#' #   list_ch = ch_filt
#' # )
#'
#' @export
pc_norm <- function(
  so,
  mdn,
  mtd = "slo",
  list_ch,
  span1 = 0.05,
  mcc = 2,
  its = 3
) {
  # Input data
  d <- so
  d <- Seurat::AddMetaData(
    d,
    metadata = seq.int(1, nrow(d@meta.data), 1),
    col.name = "ID"
  )
  d2 <- SeuratObject::FetchData(
    d,
    vars = c(
      "X",
      "Y",
      "ID",
      list_ch
    )
  )
  sp1 <- span1
  it1 <- its
  mc1 <- mcc
  # z-score Normalization
  if(mtd == "zsc") { # nolint
    print("Normalizing data by z-score...")
    ## Z-score only
    #### formula (x - mean)/sd # nolint
    qc2_filt[["value.z"]] <- (
      qc2_filt[["value"]] - mean(qc2_filt[["value"]])
    ) /
      sd(qc2_filt[["value"]])
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
  }
  if(mtd == "slo") { # nolint
    # sLOESS Normalization
    d2 <- setNames(data.frame(
      d2[, c(1:3)],
      dplyr::bind_cols(
        parallel::mclapply(
          mc.cores = mc1,
          seq.int(4, ncol(d2), 1),
          function(x) {
            dl <- d2[, c(1:3, x)]
            # Filter 0 values to ignore in LOESS calculation
            dl2 <- dl[dl[[4]] > 0, ]
            # Split data along both axes
            dl2x <- dl2[, c(1, 3:4)]
            dl2x <- dl2x[order(dl2x[[1]]), ]
            dl2y <- dl2[, c(2, 3:4)]
            dl2y <- dl2y[order(dl2y[[1]]), ]
            # Fit LOESS to X and Y axes
            fitx <- as.data.frame(
              lowess(
                x = dl2x[[1]],
                y = dl2x[[3]],
                # span
                f = sp1,
                iter = it1,
                delta = 0
              )
            )
            fitx <- setNames(cbind(
              dl2x[, 2],
              fitx
            ), c(names(dl2x)[c(2, 1)], "int"))
            fity <- as.data.frame(
              lowess(
                x = dl2y[[1]],
                y = dl2y[[3]],
                # span
                f = sp1,
                iter = it1,
                delta = 0
              )
            )
            fity <- setNames(cbind(
              dl2y[, 2],
              fity
            ), c(names(dl2y)[c(2, 1)], "int"))
            # Calculate fitted average between axes
            fitc <- dplyr::left_join(
              fitx,
              fity,
              by = "ID"
            )
            fitc[["int.avg"]] <- (
              fitc[["int.x"]] +
                fitc[["int.y"]]
            ) / 2
            fitc <- fitc[order(fitc[["ID"]]), ][, -c(2:5)]
            # Normalize to fitted model
            dl2[["norm"]] <- (dl2[[4]] / fitc[["int.avg"]]) *
              mean(dl2[[4]])
            dl2 <- dplyr::left_join(
              dl,
              dl2,
              by = "ID"
            )[, "norm"]
            dl2[is.na(dl2)] <- 0
            dl2 <- setNames(as.data.frame(dl2), paste("X", 1, sep = "."))
            return(dl2) # nolint
          }
        )
      )
    ), names(d2))
    # Create Seurat object
    ## Merge meta.data
    d3 <- dplyr::left_join(
      d@meta.data,
      d2,
      by = "ID"
    )
    d3 <- setNames(dplyr::select(
      d3,
      c("Slide", "Code", "Group", "X.x", "Y.x", "ID"),
      (mdn + 3):ncol(d3)
    ),
    c("Slide", "Code", "Group", "X", "Y", "ID", names(d3)[(mdn + 3):ncol(d3)]))
    d2 <- Seurat::CreateSeuratObject(
      counts = t(as.matrix(d3[, 7:ncol(d3)])),
      meta.data = d3[, c(1:6)],
      assay = "PC.norm"
    )
  }
  return(d2) # nolint
}

#' Phenocycler Normalization Summary
#'
#' Helper function for generating a short summary and images for
#' a normalized Phenocycler dataset.
#'
#' @param type1 Data type; provide a character string indicating
#' if the data have been normalized.
#' @param so Input Seurat object.
#' @param samp_name Sample name.
#' @param mcc Number of cores to use.
#' @return A summary of normalization performance for the
#' selected Seurat Object.
#' @examples
#'
#' # d <- pc_qc_sum(
#' #   type1 = "raw",
#' #   so = dser,
#' #   samp_name = "sample1"
#' # )
#'
#' @export
pc_qc_sum <- function(
  type1 = "norm",
  so,
  samp_name,
  mcc = 2
) {
  parallel::mclapply(
    mc.cores = mcc,
    seq.int(1, length(rownames(so)), 1),
    function(x) {
      ifelse(
        dir.exists(paste("analysis/qc/", samp_name, "/", sep = "")) == FALSE,
        dir.create(paste("analysis/qc/", samp_name, "/", sep = "")),
        print(
          paste("Saving QC cellmap in analysis/qc/", samp_name, "/", sep = "")
        )
      )
      ggplot2::ggsave(
        paste(
          "analysis/qc/", samp_name, "/cellmap_intensity_",
          gsub("/| |-", "_", rownames(so)[[x]]), "_", type1, ".png", sep = ""
        ),
        pc_qc_ch(
          so = so,
          ch_id = rownames(so)[[x]]
        ),
        height = 6,
        width = 18,
        dpi = 300
      )
    }
  )
  ld_sum <- dplyr::bind_rows(
    parallel::mclapply(
      mc.cores = mcc,
      seq.int(1, length(rownames(so)), 1),
      function(x) {
        pc_qc_ch(
          type = "sum",
          so = so,
          ch_id = rownames(so)[[x]]
        )
      }
    )
  )
  write.table(
    ld_sum,
    paste("analysis/qc/", samp_name, "/qc_summary_", type1, ".txt", sep = ""),
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t"
  )
  return(ld_sum) # nolint
}
